# ------------------------------------------------------
# Merge ERA20c and ERA5 by adjusting the mean over 1979.
# Average the results into monthly means
# Compute wind stress from wind fields
# ------------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
import mod_gentools as gentools
def main():
    settings = {}
    settings['dir_data']   = os.getenv('HOME') + '/Data/'
    settings['fn_ERA5']    = settings['dir_data'] + 'Reanalyses/ERA5/ERA5.nc'
    settings['fn_ERA20c']  = settings['dir_data'] + 'Reanalyses/ERA20c/ERA20c_1900_1980_sfc.nc'
    settings['fn_save']    = settings['dir_data'] + 'SA/ERA_monthly.nc'
    settings['time'] = np.arange(1900+1/24,2019+1/24,1/12)

    ERA5   = read_dataset(settings['fn_ERA5'])
    ERA20c = read_dataset(settings['fn_ERA20c'])
    ERA    = merge_datasets(ERA5, ERA20c, settings)
    save_data(ERA,settings)
    return

def read_dataset(fn):
    Data = {}
    file_handle = Dataset(fn)
    file_handle.set_auto_mask(False)
    Data['time'] = 1900 + file_handle.variables['time'][:]/365/24
    Data['lon'] = file_handle.variables['longitude'][:]
    Data['lat'] = np.flipud(file_handle.variables['latitude'][:])
    Data['uwind'] = np.fliplr(file_handle.variables['u10'][:])
    Data['vwind'] = np.fliplr(file_handle.variables['v10'][:])
    Data['mslp'] = np.fliplr(file_handle.variables['msl'][:])
    file_handle.close()
    return(Data)

def merge_datasets(ERA5,ERA20c,settings):
    # overlap indices
    ERA5_ovl = (ERA5['time'] >=1979) & (ERA5['time'] < 1980)
    ERA20c_ovl = (ERA20c['time'] >=1979) & (ERA20c['time'] < 1980)

    # Remove bias in 1979
    ERA20c['mslp'] = ERA20c['mslp'] - ERA20c['mslp'][ERA20c_ovl,:,:].mean(axis=0)[np.newaxis,:,:] + ERA5['mslp'][ERA5_ovl,:,:].mean(axis=0)[np.newaxis,:,:]
    ERA20c['uwind'] = ERA20c['uwind'] - ERA20c['uwind'][ERA20c_ovl,:,:].mean(axis=0)[np.newaxis,:,:] + ERA5['uwind'][ERA5_ovl,:,:].mean(axis=0)[np.newaxis,:,:]
    ERA20c['vwind'] = ERA20c['vwind'] - ERA20c['vwind'][ERA20c_ovl,:,:].mean(axis=0)[np.newaxis,:,:] + ERA5['vwind'][ERA5_ovl,:,:].mean(axis=0)[np.newaxis,:,:]

    total_time = gentools.monthly_time(1900,2018)
    mslp  = np.vstack([ERA20c['mslp'][:-12,:,:],ERA5['mslp']])
    uwind = np.vstack([ERA20c['uwind'][:-12,:,:],ERA5['uwind']])
    vwind = np.vstack([ERA20c['vwind'][:-12,:,:],ERA5['vwind']])

    # From wind speed to wind stress
    uws = (0.8 + 0.065 * np.sqrt(uwind ** 2 + vwind ** 2)) * uwind * np.sqrt(uwind ** 2 + vwind ** 2)
    vws = (0.8 + 0.065 * np.sqrt(uwind ** 2 + vwind ** 2)) * vwind * np.sqrt(uwind ** 2 + vwind ** 2)

    ERA = {}
    ERA['lat']  = ERA5['lat']
    ERA['lon']  = ERA5['lon']
    ERA['time'] = settings['time']
    ERA['mslp'] = mslp
    ERA['uws']  = uws
    ERA['vws']  = vws
    return(ERA)

def save_data(ERA,settings):
    file_handle = Dataset(settings['fn_save'], 'w')
    file_handle.createDimension('lon', len(ERA['lon']))
    file_handle.createDimension('lat', len(ERA['lat']))
    file_handle.createDimension('time', len(ERA['time']))
    file_handle.createVariable('lon', 'f4', ('lon',),zlib=True)[:] = ERA['lon']
    file_handle.createVariable('lat', 'f4', ('lat',),zlib=True)[:] = ERA['lat']
    file_handle.createVariable('time', 'i4', ('time',),zlib=True)[:] = ERA['time']
    file_handle.createVariable('mslp', 'f4', ('time', 'lat', 'lon',),zlib=True)[:] = ERA['mslp']
    file_handle.createVariable('uws', 'f4', ('time', 'lat', 'lon',),zlib=True)[:] = ERA['uws']
    file_handle.createVariable('vws', 'f4', ('time', 'lat', 'lon',),zlib=True)[:] = ERA['vws']
    file_handle.close()
    return
