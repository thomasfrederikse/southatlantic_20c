# ------------------------------------
# Read 20c trends from F2020 and A2018
# Compute mean and standard errors
# Save for GMT
# ------------------------------------
import numpy as np
from netCDF4 import Dataset
import mod_gentools as gentools
from scipy.interpolate import interp2d
import os

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['fn_grd_trends'] = settings['dir_project'] + 'GMT/PD_20c/grd_trends.nc'
    settings['fn_mask']      = settings['dir_data'] + 'GRACE/JPL_mascon/mask.npy'

    mask = np.load(settings['fn_mask'],allow_pickle=True).all()

    mask_SA = (mask['basin']==4)*1.0
    mask_glb = (~mask['land'])*1.0

    adhikari   = read_adhikari(settings)
    frederikse = read_frederikse(settings)

    area = gentools.grid_area(mask['lat'], mask['lon'])

    # Compute mean over SA
    adhikari_SA   = np.sum(mask_SA * area * adhikari) / np.sum(mask_SA * area)
    frederikse_SA   = np.sum(mask_SA * area * frederikse) / np.sum(mask_SA * area)

    adhikari_glb   = np.sum(mask_glb * area * adhikari) / np.sum(mask_glb * area)
    frederikse_glb   = np.sum(mask_glb * area * frederikse) / np.sum(mask_glb * area)

    file_handle = Dataset(settings['fn_grd_trends'], 'w')
    file_handle.createDimension('lon', len(mask['lon']))
    file_handle.createDimension('lat', len(mask['lat']))
    file_handle.createVariable('lon', 'f4', ('lon',), zlib=True)[:] = mask['lon']
    file_handle.createVariable('lat', 'f4', ('lat',), zlib=True)[:] = mask['lat']
    file_handle.createVariable('adhikari', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = adhikari
    file_handle.createVariable('frederikse', 'f4', ('lat', 'lon',), zlib=True, complevel=4)[:] = frederikse
    file_handle.close()


    return



def read_adhikari(settings):
    fh = Dataset(settings['dir_data']+'sle/Surendra_20c/fp_20c.nc','r')
    fh.set_auto_mask(False)
    adhikari = fh['rsl_mean'][:]
    return(adhikari)

def read_frederikse(settings):
    fh = Dataset(settings['dir_data']+'Budget_20c/results/spatial_trends.nc','r')
    fh.set_auto_mask(False)
    frederikse = fh['grd_1900_2018'][:]
    return(frederikse)


