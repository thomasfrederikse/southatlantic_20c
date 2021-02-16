# ---------------------------------------------------------------------------
# Compute vertical land motion from GIA and contemporary mass redistribution
# from GRACE JPL mascons and Caron et al. (2018) GIA ensemble
# ---------------------------------------------------------------------------
import numpy as np
import pySLE3
import os
from netCDF4 import Dataset
import multiprocessing as mp
import ctypes as ct
import datetime as dt
import sys

def main():
    global settings, dt_start
    set_settings()
    read_grace()
    read_mask()
    read_love()
    dt_start = dt.datetime.now().replace(microsecond=0)
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(func=process_ensemble_member, iterable=settings['ens_range'])
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']    = os.getenv('HOME')+'/Data/'
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['dir_save_ens'] = os.getenv('HOME')+'/Data/SA/rad_grd_ens/'
    settings['fn_grace']    = settings['dir_data'] + 'GRACE/JPL_mascon/JPL_mascon_RL06v02_CRI_noGIA_noEQ_noseas.nc'
    settings['fn_slm']      = settings['dir_data'] + 'GRACE/JPL_mascon/mask.npy'
    settings['fn_mscn_coords'] = settings['dir_data'] + 'GRACE/JPL_mascon/mascon_coords.npy'
    settings['fn_gia_ens_rad'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rad_ens_05.nc'
    settings['fn_gia_ens_ewh'] = settings['dir_data'] + 'GIA/Caron/Ensemble/ewh_ens_05.nc'
    settings['fn_love']        = settings['dir_data'] + 'Budget_20c/grd_prep/love.npy'
    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
        settings['dir_scratch'] = os.getenv('HOME') + '/Scratch/'
        settings['ens_range'] = np.arange(100)
    else:
        settings['nproc'] = 4
        settings['dir_scratch'] = '/tmp/'
        settings['ens_range'] = np.arange(int(sys.argv[1]),int(sys.argv[2]))
    return

def read_grace():
    # Read GRACE data and interpolate on the common monthly grid
    print(' Reading GRACE data...')
    global settings
    global lat, lon, time, ewh_nogia, ewh_ste
    # Read GRACE data
    file_handle = Dataset(settings['fn_grace'],'r')
    file_handle.set_auto_mask(False)
    lon = mp_filled_float(file_handle.variables['lon'][:])
    lat = mp_filled_float(file_handle.variables['lat'][:])
    time  = file_handle.variables['time'][:]
    ewh_nogia   = file_handle.variables['ewh_noseas'][:]
    ewh_ste_grid = file_handle.variables['ewh_ste'][:]
    file_handle.close()
    ewh_ste = mp_filled_float(grid2mascon_3d(lat,lon,ewh_ste_grid,settings))
    return

def read_mask():
    global settings,slm
    slm = mp_filled_bool(~np.load(settings['fn_slm'],allow_pickle=True).all()['land'])
    return

def read_love():
    global love, settings
    love = np.load(settings['fn_love'],allow_pickle=True).all()
    return

def process_ensemble_member(ens):
    global settings,dt_start
    dt_now = dt.datetime.now().replace(microsecond=0)
    print('   Ensemble '+str(ens+1)+'/'+str(len(settings['ens_range']))+' Elapsed time:',(dt_now - dt_start))
    print('    GIA '+str(ens+1))
    ewh_gia, rad_gia = read_gia_ensemble_member(ens)
    print('    EWH '+str(ens+1))
    ewh_ptb  = perturb_ewh_estimate(ewh_gia)
    print('    SLE '+str(ens+1))
    rad_grd = solve_sle(ewh_ptb) + rad_gia * (time-time.mean())[:,np.newaxis,np.newaxis]
    print('    save '+str(ens+1))
    save_ens(ens,rad_grd)
    print('   Ensemble '+str(ens+1)+'/'+str(len(settings['ens_range']))+' done. Took:',(dt.datetime.now().replace(microsecond=0) - dt_now))
    return

def read_gia_ensemble_member(ens):
    global settings,lat,lon
    file_handle = Dataset(settings['fn_gia_ens_rad'],'r')
    file_handle.set_auto_mask(False)
    rad_gia   = file_handle.variables['rad'][ens,...].astype(np.float32)
    file_handle.close()
    file_handle = Dataset(settings['fn_gia_ens_ewh'],'r')
    file_handle.set_auto_mask(False)
    ewh_gia   = file_handle.variables['ewh'][ens,...].astype(np.float32)
    file_handle.close()
    ewh_gia = (masconize_2d(lat,lon,ewh_gia,settings)).astype(np.float32)
    return(ewh_gia, rad_gia)

def perturb_ewh_estimate(ewh_gia):
    # Correct GRACE estimate with GIA ensemble member and perturb with measurement noise
    global lat, lon, time, ewh_nogia, ewh_ste, slm,settings
    np.random.seed() # Reset random number generator to avoid clustering in MP
    rnd_mscn = np.random.normal(0,1,ewh_ste.shape)*ewh_ste    # Random pertubations from GRACE uncertainty
    ewh_grd = ewh_nogia + mascon2grid_3d(lat,lon,rnd_mscn, settings) - (ewh_gia*(time-time.mean())[:,np.newaxis,np.newaxis])
    ewh_grd[:,slm] = 0
    return(ewh_grd)

def solve_sle(ewh_ptb):
    global slm, lat,lon,time,love, settings
    # Solve the sea-level equation for each perturbed GRACE field
    # return RSL and radial deformation
    # D/O 90, rotation, CM
    rad_grd = 1000*pySLE3.solve(lat=lat, lon=lon, time=time, load=ewh_ptb, slm=slm,love=love,lmax=179,rotational_feedback=True,geoid_out=False,rad_out=True,rsl_out=False,barystatic_out=False,lod_out=False,verbose=False)['rad']
    return(rad_grd)

def save_ens(ens,rad_grd):
    global settings
    global lat,lon,time
    fname = settings['dir_save_ens'] + 'rad_grd_ens_'+str(ens).zfill(4)+'.nc'
    file_handle = Dataset(fname, 'w')
    file_handle.createDimension('lon', len(lon))
    file_handle.createDimension('lat', len(lat))
    file_handle.createDimension('time', len(time))

    file_handle.createVariable('lon', 'f4', ('lon',), zlib=True)[:] = lon
    file_handle.createVariable('lat', 'f4', ('lat',), zlib=True)[:] = lat
    file_handle.createVariable('time', 'f4', ('time',), zlib=True)[:] = time

    # 3D fiels
    file_handle.createVariable('rad_grd', 'i2', ('time', 'lat', 'lon',), zlib=True, complevel=4)[:] = (np.rint(20*rad_grd)).astype(int)
    file_handle.variables['rad_grd'].setncattr('scale_factor', 0.05)
    file_handle.close()
    return

def grid2mascon_3d(lat,lon,grid,settings):
    print('   Transforming field into mascons...')
    coords = np.load(settings['fn_mscn_coords'])
    mscn = np.zeros([coords.shape[0],grid.shape[0]],dtype=np.float32)
    for k in range(len(coords)):
        lat_acc = np.where((lat >= coords[k,0]) & (lat < coords[k,1]))[0]
        lon_acc = np.where((lon >= coords[k,2]) & (lon < coords[k,3]))[0]
        weight = np.cos(np.deg2rad(lat[lat_acc])) / np.mean(np.cos(np.deg2rad(lat[lat_acc])))  # Weight by cos lat
        mscn[k,:] = np.nanmean(weight[:,np.newaxis]*grid[:,lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1] + 1],axis=(1,2))
    return(mscn)

def mascon2grid_3d(lat,lon,mscn,settings):
    grid = np.zeros([mscn.shape[1],len(lat),len(lon)],dtype=np.float32)
    coords = np.load(settings['fn_mscn_coords'])
    for k in range(len(coords)):
        lat_acc = np.where((lat >= coords[k,0]) & (lat < coords[k,1]))[0]
        lon_acc = np.where((lon >= coords[k,2]) & (lon < coords[k,3]))[0]
        grid[:,lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1] + 1] = mscn[k,:][:,np.newaxis,np.newaxis]
    return(grid)


def masconize_2d(lat,lon,field,settings):
    coords = np.load(settings['fn_mscn_coords'])
    field_mscn = np.zeros([len(lat), len(lon)],dtype=np.float32)
    for k in range(len(coords)):
        lat_acc = np.where((lat >= coords[k,0]) & (lat < coords[k,1]))[0]
        lon_acc = np.where((lon >= coords[k,2]) & (lon < coords[k,3]))[0]
        weight = np.cos(np.deg2rad(lat[lat_acc])) / np.mean(np.cos(np.deg2rad(lat[lat_acc])))  # Weight by cos lat
        field_mscn[lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1]+1] = np.nanmean(weight[:,np.newaxis] * field[lat_acc[0]:lat_acc[-1]+1,lon_acc[0]:lon_acc[-1]+1])
    return(field_mscn)


def mp_empty_float(shape):
    shared_array_base = mp.RawArray(ct.c_float, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_empty_int(shape):
    shared_array_base = mp.RawArray(ct.c_int, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_empty_bool(shape):
    shared_array_base = mp.RawArray(ct.c_bool, int(np.prod(shape)))
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_float(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_float, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_int(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_int, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

def mp_filled_bool(input_array):
    shape = input_array.shape
    shared_array_base = mp.RawArray(ct.c_bool, input_array.flatten())
    shared_array = np.ctypeslib.as_array(shared_array_base).reshape(*shape)
    return shared_array

if __name__ == '__main__':
    main()
