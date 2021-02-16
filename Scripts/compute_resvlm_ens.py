# --------------------------------------------
# Compute ensemble of residual VLM trends
# from UNR GPS database and GRD from GIA/GRACE
# following Frederikse et al. 2019
# --------------------------------------------
import numpy as np
import os
from netCDF4 import Dataset
import multiprocessing as mp
import ctypes as ct
import datetime as dt
import urllib.request
import mod_midas_py_ens as midas_py_ens
import mod_midas_py_single as midas_py_single
from scipy.interpolate import interp1d

def main():
    set_settings()
    read_station_list()
    download_gps_tseries()
    estimate_resvlm_trends()
    save_resvlm()
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']    = os.getenv('HOME')+'/Data/'
    settings['dir_scratch'] = os.getenv('HOME')+'/Scratch/'
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['dir_grd_ens'] = os.getenv('HOME')+'/Data/SA/rad_grd_ens/'
    settings['dir_gps'] = settings['dir_project'] + 'Data/GPS/'
    settings['dir_ens'] = settings['dir_data']+'SA/rad_grd_ens/'

    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
        settings['dir_scratch'] = os.getenv('HOME') + '/Scratch/'
    else:
        settings['nproc'] = 200
        settings['dir_scratch'] = '/tmp/'
    settings['ens_range'] = np.arange(5000)
    settings['fn_station_list'] = settings['dir_project'] + 'Data/station_list.npy'
    settings['fn_residual_vlm'] = settings['dir_project'] + 'Data/residual_vlm.npy'
    return

def read_station_list():
    global settings, station_list
    station_list = np.load(settings['fn_station_list'],allow_pickle=True)
    return

def download_gps_tseries():
    global settings, station_list
    print(' Downloading data holdings...')
    fn_dataholdings = settings['dir_scratch']+'DataHoldings.txt'
    null = urllib.request.urlretrieve('http://geodesy.unr.edu/NGLStationPages/DataHoldings.txt',fn_dataholdings)

    # Read filelist
    statlist = np.loadtxt(fn_dataholdings,skiprows=1,usecols=(0,1,2),dtype='object')
    statlist[:,1:3] = statlist[:,1:3].astype(float)
    os.remove(fn_dataholdings)

    # Steps database
    print(' Downloading steps...')

    fn_steplist = settings['dir_scratch']+'steps.txt'
    fn_decyear = settings['dir_scratch']+'decyear.txt'
    void = urllib.request.urlretrieve('http://geodesy.unr.edu/NGLStationPages/steps.txt', fn_steplist)
    void = urllib.request.urlretrieve('http://geodesy.unr.edu/NGLStationPages/decyr.txt', fn_decyear)
    steplist = np.loadtxt(fn_steplist, delimiter='  ', dtype=object, usecols=(0,1,2))
    steplist[:,2] = steplist[:,2].astype(int)
    decyr_lookup = np.loadtxt(fn_decyear, delimiter=' ', dtype=object, usecols=(0,1))
    decyr_dict = dict(zip(decyr_lookup[:,0],decyr_lookup[:,1].astype(float)))
    step_stats = steplist[:,0]
    step_times = steplist[:,1]
    step_decyears = np.zeros(len(step_times))
    for i in range(len(step_times)):
        step_decyears[i] = decyr_dict[step_times[i]]
    steplist = np.vstack([step_stats,step_decyears]).T
    os.remove(fn_steplist)
    os.remove(fn_decyear)

    # List of downloads
    gps_list = []
    for stat in station_list:
        for gps in stat['gps_list']:
            gps_list.append(gps)
    print(' Downloading GPS stations...')
    for gps in gps_list:
        print('   Processing '+gps+'...')
        gps_data = {}
        idx = np.where(statlist[:, 0] == gps)[0][0]
        gps_data['lat'] = statlist[idx,1]
        gps_data['lon'] = statlist[idx,2]
        url = 'http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/' + gps + '.tenv3'
        fn_gps = settings['dir_scratch']+gps+'.txt'
        null = urllib.request.urlretrieve(url,fn_gps)  # Download
        data_raw = np.loadtxt(fn_gps, usecols=(2, 12), skiprows=1)
        time_acc = (data_raw[:, 0] > 2002.25) & (data_raw[:, 0] < 2020.6)
        gps_data['time']   = data_raw[time_acc, 0]
        gps_data['height'] = 1000*data_raw[time_acc, 1]
        os.remove(fn_gps)
        # Find steps
        idx = np.where(steplist[:, 0] == gps)[0]
        if len(idx) > 0:
            gps_data['steplist'] = steplist[idx,1]
        np.save(settings['dir_gps']+gps+'.npy',gps_data)
    return

def estimate_resvlm_trends():
    global settings, station_list
    global gps_coords,gps_rad_array
    # Get GPS coordinates
    print(' Obtaining GPS coordinates...')
    gps_list = []
    for stat in station_list:
        for gps in stat['gps_list']:
            gps_list.append(gps)
    lat = np.arange(-89.75,90.25,0.5)
    lon = np.arange(0.25,360.25,0.5)
    gps_coords = mp_empty_int([len(gps_list),2])
    for idx,gps in enumerate(gps_list):
        gps_data = np.load(settings['dir_gps']+gps+'.npy',allow_pickle=True).all()
        gps_coords[idx,0] = np.argmin(np.abs(gps_data['lat']-lat))
        gps_coords[idx,1] = np.argmin(np.abs(gps_data['lon']-lon))

    # Fill rad array
    print(' Reading GRD time series...')
    gps_rad_array = mp_empty_float([187,len(gps_list),len(settings['ens_range'])])
    pool = mp.Pool(settings['nproc'])
    out  = pool.map(func=read_grd_ens, iterable=settings['ens_range'])

    # Trends
    grace_time = Dataset(settings['dir_ens']+'rad_grd_ens_0000.nc','r')["time"][:]._get_data()
    global station_resvlm_trends
    station_resvlm_trends = mp_empty_float([len(station_list),len(settings['ens_range'])])
    print(' Computing resvlm trends...')
    for idx_stat, stat in enumerate(station_list):
        print('   Station '+stat['name']+'...')
        if stat['name'] == 'Kerguelen':
            indiv_resvlm_trend = np.zeros([len(stat['gps_list'])+1,len(settings['ens_range'])])
            indiv_resvlm_sterr = np.zeros([len(stat['gps_list'])+1,len(settings['ens_range'])])
        else:
            indiv_resvlm_trend = np.zeros([len(stat['gps_list']),len(settings['ens_range'])])
            indiv_resvlm_sterr = np.zeros([len(stat['gps_list']),len(settings['ens_range'])])
        for idx_gps, gps in enumerate(stat['gps_list']):
            gps_data = np.load(settings['dir_gps']+gps+'.npy',allow_pickle=True).all()
            gps_index = gps_list.index(gps)
            resvlm_ens = gps_data['height'][:,np.newaxis] - interp1d(grace_time,gps_rad_array[:,gps_index,:],axis=0,kind='linear',fill_value='extrapolate')(gps_data['time'])
            time_pad = np.zeros(9999)
            time_pad[:len(gps_data['time'])] = gps_data['time']
            ts_pad = np.zeros([9999,5000])
            ts_pad[:len(gps_data['time']),:len(settings['ens_range'])] = resvlm_ens
            if 'steplist' in gps_data:
                steplist = np.zeros(100)
                steplist[:len(gps_data['steplist'])] = gps_data['steplist']
                resvlm_ens_trend, resvlm_ens_sterr = midas_py_ens.midas_py(len(gps_data['time']), time_pad, ts_pad, len(settings['ens_range']), len(gps_data['steplist']), steplist)
            else:
                resvlm_ens_trend, resvlm_ens_sterr = midas_py_ens.midas_py(len(gps_data['time']), time_pad, ts_pad, len(settings['ens_range']), 0, np.zeros(100))
            indiv_resvlm_trend[idx_gps,:] = resvlm_ens_trend[:len(settings['ens_range'])]
            indiv_resvlm_sterr[idx_gps,:] = resvlm_ens_sterr[:len(settings['ens_range'])]

        # Add DORIS for Kerguelen
        if stat['name'] == 'Kerguelen':
            doris_tstart = 2007.3142
            doris_tstop  = 2014.9227
            doris_vlm_trend = -0.1
            doris_vlm_sterr =  0.2
            grace_time_acc = (grace_time>doris_tstart) & (grace_time<doris_tstop)
            grd_trend = np.zeros(len(settings['ens_range']))
            amat = np.ones([grace_time_acc.sum(), 2])
            amat[:, 1] = grace_time[grace_time_acc] - grace_time[grace_time_acc].mean()
            amat_T = amat.T
            amat_sq = np.linalg.inv(np.dot(amat_T, amat))
            for i in settings['ens_range']: grd_trend[i] = np.dot(amat_sq, np.dot(amat_T, gps_rad_array[grace_time_acc,gps_index,i]))[1]
            indiv_resvlm_trend[-1,:] = doris_vlm_trend - grd_trend
            indiv_resvlm_sterr[-1, :]= doris_vlm_sterr

        # Weigh inverse of
        weight = 1/(indiv_resvlm_sterr**2) / (1/(indiv_resvlm_sterr**2)).sum(axis=0)
        rnd_ptb = np.random.normal(loc=0,scale=1,size=indiv_resvlm_trend.shape)
        station_resvlm_trends[idx_stat,:] = (weight*(indiv_resvlm_trend+rnd_ptb*indiv_resvlm_sterr)).sum(axis=0)
    return

def read_grd_ens(ens):
    print('   GRD ensemble '+str(ens)+'...')
    global gps_coords,gps_rad_array
    global settings
    fn = settings['dir_ens'] + 'rad_grd_ens_'+str(ens).zfill(4)+'.nc'
    fh = Dataset(fn,'r')
    fh.set_auto_mask(False)
    gps_rad_array[:,:,ens] = fh['rad_grd'][:][:,gps_coords[:,0],gps_coords[:,1]]
    fh.close()
    return

def save_resvlm():
    print(' Saving...')
    global settings, station_list, station_resvlm_trends
    resvlm = np.zeros(len(station_list),dtype=dict)
    for idx, stat in enumerate(station_list):
        resvlm[idx] = {}
        resvlm[idx]['name'] = stat['name']
        resvlm[idx]['gps_list'] = stat['gps_list']
        resvlm[idx]['resvlm'] = station_resvlm_trends[idx,:]
    np.save(settings['fn_residual_vlm'],resvlm)
    return

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





