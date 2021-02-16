# --------------------------------------------------
# Compute trends in gridded sea-level observations:
# EN4
# CZ16
# WOA
# I17
# Compute SA trend and global trend 1957-2018
# --------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import mod_gentools as gentools
import os
from scipy.interpolate import interp2d

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'

    settings['fn_mask']      = settings['dir_data'] + 'GRACE/JPL_mascon/mask.npy'
    settings['fn_CZ16'] = settings['dir_data'] + 'Steric/Cheng/Cheng_1940_2019.nc'
    settings['fn_I17'] = settings['dir_data'] + 'Steric/I17/I17_1955_2019.nc'
    settings['fn_WOA'] = settings['dir_data'] + 'Steric/Levitus/Levitus_1957_2019.nc'
    settings['fn_EN4'] = settings['dir_data'] + 'Steric/EN4/EN4_L09_1950_2019.nc'
    settings['prod_list'] = ['WOA','I17','CZ16','EN4']
    settings['fn_steric_trends'] = settings['dir_project'] + 'GMT/steric_trends.nc'

    steric_trend = {}
    for prod in settings['prod_list']:
        steric_trend[prod] = process_indiv(settings['fn_'+str(prod)],settings)

    # Save data
    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask_SA = (mask['basin']==4)*1.0
    mask_SA = np.hstack([mask_SA[:, 360:], mask_SA[:, :360]])
    mask['lon'] = np.hstack([mask['lon'][360:], mask['lon'][:360]])
    mask['lon'][mask['lon']>180.1] = mask['lon'][mask['lon']>180.1]-360

    file_handle = Dataset(settings['fn_steric_trends'], 'w')
    file_handle.createDimension('lon', len(mask['lon']))
    file_handle.createDimension('lat', len(mask['lat']))
    file_handle.createVariable('lon', 'f4', ('lon',), zlib=True)[:] = mask['lon']
    file_handle.createVariable('lat', 'f4', ('lat',), zlib=True)[:] = mask['lat']
    file_handle.createVariable('EN4', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = steric_trend['EN4']['trend']
    file_handle.createVariable('I17', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = steric_trend['I17']['trend']
    file_handle.createVariable('CZ16', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = steric_trend['CZ16']['trend']
    file_handle.createVariable('WOA', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = steric_trend['WOA']['trend']
    file_handle.createVariable('mask', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = mask_SA
    file_handle.close()

    # Save trends
    for prod in settings['prod_list']:
        fn_SA = settings['dir_project'] + 'GMT/steric/'+prod+'_SA.txt'
        fn_glb = settings['dir_project'] + 'GMT/steric/'+prod+'_glb.txt'
        np.savetxt(fn_SA,np.array([steric_trend[prod]['SA']]),fmt='%4.2f')
        np.savetxt(fn_glb,np.array([steric_trend[prod]['glb']]),fmt='%4.2f')


    return

def process_indiv(fn,settings):
    result = {}
    fh = Dataset(fn,'r')
    fh.set_auto_mask(False)
    lat = fh['lat'][:]
    lon = fh['lon'][:]
    time = fh['time'][:]
    slm = fh['slm'][:]


    time_acc = (time>1957) & (time<2019)
    steric = fh['totalsteric_2d'][time_acc,...]
    time = time[time_acc]
    steric = np.dstack([steric[:,:,180:], steric[:,:, :180]])
    slm = np.hstack([slm[:, 180:], slm[:, :180]])
    lon = np.hstack([lon[180:], lon[:180]])
    lon[lon>180.1] = lon[lon>180.1]-360

    steric_trend_grid = gentools.field_trend(time,steric)

    mask = np.load(settings['fn_mask'],allow_pickle=True).all()
    mask_SA = (mask['basin']==4)*1.0
    mask_SA = np.hstack([mask_SA[:, 360:], mask_SA[:, :360]])
    mask_lon = mask['lon']
    mask_lon = np.hstack([mask_lon[360:], mask_lon[:360]])
    mask_lon[mask_lon>180.1] = mask_lon[mask_lon>180.1]-360

    msk_interp = np.rint(interp2d(mask_lon, mask['lat'], mask_SA, kind='linear', bounds_error=False, fill_value=0)(lon, lat)) * slm

    area = gentools.grid_area(lat,lon)
    result['glb'] = np.nansum(slm*steric_trend_grid*area)/np.sum(slm*area)
    result['SA'] = np.nansum(msk_interp*steric_trend_grid*area)/np.sum(msk_interp*area)

    steric_trend_grid[np.isnan(steric_trend_grid)] = 0
    result['trend'] = interp2d(lon,lat, steric_trend_grid, kind='linear', bounds_error=False)(mask_lon, mask['lat'])
    return(result)






