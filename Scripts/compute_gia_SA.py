# Read Carons ensemble and derive mean and standard deviation
import numpy as np
from netCDF4 import Dataset
import mod_gentools as gentools
from scipy.interpolate import interp2d
import os

def main():
    rsl = read_Caron()
    area = gentools.grid_area(rsl['lat'], rsl['lon'])

    mask = np.load(os.getenv('HOME') + '/Data/GRACE/JPL_mascon/mask.npy',allow_pickle=True).all()
    mask_SA = (mask['basin']==4)*1.0
    mask_SA = np.rint(interp2d(mask['lon'], mask['lat'], mask_SA, kind='linear', bounds_error=False, fill_value=0)(rsl['lon'], rsl['lat']))

    # Compute mean over SA
    GIA_rsl_SA_ens   = np.sum(mask_SA * area * rsl['rsl'],axis = (1,2)) / np.sum(mask_SA * area)
    GIA_rsl_SA_mean  = np.sum(GIA_rsl_SA_ens*rsl['probability'])
    GIA_rsl_SA_sterr = np.sqrt(np.sum((GIA_rsl_SA_ens-GIA_rsl_SA_mean)**2*rsl['probability']))

    # Mean and standard error of field
    GIA_rsl_mean = np.sum(rsl['probability'][:,np.newaxis,np.newaxis] * rsl['rsl'],axis=0 )
    GIA_rsl_sterr = np.sqrt(np.sum(rsl['probability'][:,np.newaxis,np.newaxis] * (rsl['rsl'] - GIA_rsl_mean[np.newaxis,:,:])**2,axis=0))


    # Fields to -180,180
    GIA_rsl_mean = np.hstack((GIA_rsl_mean[:,180:],GIA_rsl_mean[:,:180]))
    GIA_rsl_sterr = np.hstack((GIA_rsl_sterr[:,180:],GIA_rsl_sterr[:,:180]))

    rsl['lon'] =  np.hstack((rsl['lon'][180:],rsl['lon'][:180]))
    rsl['lon'][rsl['lon']>=180] = rsl['lon'][rsl['lon']>=180]-360

    # Save data
    fn_save = '/Users/tfrederi/Projects/2020_SA/GMT/GIA/GIA_mean_sterr.nc'
    file_handle = Dataset(fn_save, 'w')
    file_handle.createDimension('lon', len(rsl['lon']))
    file_handle.createDimension('lat', len(rsl['lat']))
    file_handle.createVariable('lon', 'f4', ('lon',), zlib=True, complevel=6)[:] = rsl['lon']
    file_handle.createVariable('lat', 'f4', ('lat',), zlib=True, complevel=6)[:] = rsl['lat']
    file_handle.createVariable('mean', 'f4', ('lat', 'lon',), zlib=True, complevel=6)[:] = GIA_rsl_mean
    file_handle.createVariable('sterr', 'f4', ('lat', 'lon',), zlib=True, complevel=6)[:] = GIA_rsl_sterr
    file_handle.close()
    return


def read_Caron():
    # Read Caron data
    fn_rsl  = '/Users/tfrederi/Data/GIA/Caron/Stats/indiv_rsl.nc'
    rsl = {}
    rsl['lat'] = Dataset(fn_rsl, 'r').variables['y'][:]._get_data()
    rsl['lon']   = Dataset(fn_rsl, 'r').variables['x'][:]._get_data()
    rsl['probability']   = Dataset(fn_rsl, 'r').variables['probability'][:]._get_data()
    rsl['rsl']   = Dataset(fn_rsl, 'r').variables['rsl'][:]._get_data()
    return(rsl)




