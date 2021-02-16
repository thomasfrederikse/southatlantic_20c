# -----------------------------------------------------
# Compute an ensemble of rsl estimates for each station
# Correct for nodal cycle, meteo
# Merge individual tide gauge records per station
# Correct for GIA, GRD and residual VLM
# Generate ensembles
# -----------------------------------------------------
import numpy as np
from netCDF4 import Dataset
import os
from scipy.interpolate import interp1d
import urllib.request
import multiprocessing as mp
import ctypes as ct
import mod_gentools as gentools
import mod_hector as hector

def main():
    set_settings()
    proc_mask()
    read_station_list()
    compute_station_ensembles()
    compute_station_corrected_ens()
    compute_virtual_station()
    save_ensembles()
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']    = os.getenv('HOME')+'/Data/'
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['dir_grd_ens'] = os.getenv('HOME')+'/Data/Budget_20c/grd/'
    settings['dir_sl_ensembles'] = settings['dir_project']+'/Data/sl_ensembles/'
    settings['dir_sl_plw'] = settings['dir_project']+'Data/PW_tidegauge/'

    settings['fn_mask']      = settings['dir_data'] + 'GRACE/JPL_mascon/mask.npy'
    settings['fn_gia_ens_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'

    if os.uname().nodename == 'MT-110180':
        settings['nproc'] = 4
        settings['dir_scratch'] = os.getenv('HOME') + '/Scratch/'
    elif os.uname().nodename == 'debian':
        settings['nproc'] = 4
        settings['dir_scratch'] = os.getenv('HOME') + '/Scratch/'
    else:
        settings['nproc'] = 48
        settings['dir_scratch'] = '/tmp/'

    settings['ens_range'] = np.arange(5000)
    settings['falklands_test'] = False
    settings['orig_test'] = False
    settings['fn_station_list'] = settings['dir_project'] + 'Data/station_list.npy'
    settings['fn_residual_vlm'] = settings['dir_project'] + 'Data/residual_vlm.npy'

    settings['fn_statlist_PSMSL'] = settings['dir_project'] + 'Data/sealevel/statlist_PSMSL.txt'
    settings['fn_ERA_monthly'] = settings['dir_data'] + 'SA/ERA_monthly.nc'
    settings['fn_nodal_amp'] = settings['dir_data'] + 'Nodal/Nodal.nc'
    settings['fn_dakar']       = settings['dir_data'] + 'TideGauges/Dakar/Dakar_monthly_completedwithPSMSL_Thomas.dat'
    settings['fn_falklands']   = settings['dir_project'] + 'Data/sealevel/Falklands_2020_12.txt'

    if settings['falklands_test']:
        settings['fn_sl_ensembles_stat'] = settings['dir_sl_ensembles'] + 'sl_ensembles_stat_falklands_test.npz'
        settings['fn_sl_ensembles_SA'] = settings['dir_sl_ensembles'] + 'sl_ensembles_SA_falklands_test.npz'
    elif settings['orig_test']:
        settings['fn_sl_ensembles_stat'] = settings['dir_sl_ensembles'] + 'sl_ensembles_stat_orig_test.npz'
        settings['fn_sl_ensembles_SA'] = settings['dir_sl_ensembles'] + 'sl_ensembles_SA_orig_test.npz'
    else:
        settings['fn_sl_ensembles_stat'] = settings['dir_sl_ensembles'] + 'sl_ensembles_stat.npz'
        settings['fn_sl_ensembles_SA'] = settings['dir_sl_ensembles'] + 'sl_ensembles_SA.npz'

    settings['time'] = np.arange(1900+1/24,2019+1/24,1/12)
    return


def read_station_list():
    global settings, station_list
    station_list = np.load(settings['fn_station_list'],allow_pickle=True)
    if settings['orig_test']: station_list = station_list[[0,1,2,3,5,6]]
    return

def proc_mask():
    global settings, mask_SA, area
    mask_SA = mp_filled_bool(np.load(settings['fn_mask'],allow_pickle=True).all()['basin']==4)
    lat = np.arange(-89.75, 90.25, 0.5)
    lon = np.arange(0.25, 360.25, 0.5)
    area = mp_filled_float(gentools.grid_area(lat,lon))

    return

def compute_station_corrected_ens():
    print(' Compute ensembles of sea-level data at each station...')
    global station_ensembles, station_ensembles_gia, station_ensembles_gia_grd, station_ensembles_resvlm, station_ensembles_full, settings, station_list, stat_coords
    station_ensembles_gia     = mp_filled_float(station_ensembles)
    station_ensembles_gia_grd = mp_filled_float(station_ensembles)
    station_ensembles_resvlm  = mp_filled_float(station_ensembles)
    station_ensembles_full    = mp_filled_float(station_ensembles)

    # Residual vlm
    residual_vlm = np.load(settings['fn_residual_vlm'],allow_pickle=True)
    if settings['orig_test']: residual_vlm = residual_vlm[[0,1,2,3,5,6]]

    for idx, stat in enumerate(station_list):
        assert(stat['name'] == residual_vlm[idx]['name'])
        station_ensembles_resvlm[:,:,idx]+=residual_vlm[idx]['resvlm'][np.newaxis,settings['ens_range']] * (settings['time'] - settings['time'].mean())[:,np.newaxis]
        station_ensembles_full[:,:,idx]+=residual_vlm[idx]['resvlm'][np.newaxis,settings['ens_range']] * (settings['time'] - settings['time'].mean())[:,np.newaxis]
    # GIA and GRD
    # Coords
    lat = np.arange(-89.75, 90.25, 0.5)
    lon = np.arange(0.25, 360.25, 0.5)
    stat_coords = mp_empty_int([len(station_list), 2])
    for idx, stat in enumerate(station_list):
        stat_coords[idx, 0] = np.argmin(np.abs(stat['lat'] - lat))
        stat_coords[idx, 1] = np.argmin(np.abs(stat['lon'] - lon))

    pool = mp.Pool(settings['nproc'])
    out  = pool.map(func=corr_gia_grd_ens, iterable=settings['ens_range'])
    return

def compute_virtual_station():
    print(' Compute virtual-station solutions...')
    global station_ensembles, station_ensembles_gia, station_ensembles_gia_grd, station_ensembles_resvlm, station_ensembles_full, settings, station_list
    global SA_ensembles, SA_ensembles_gia, SA_ensembles_gia_grd, SA_ensembles_resvlm, SA_ensembles_full
    station_order = determine_station_order()
    SA_ensembles = compute_virtual_station_indiv(station_order, station_ensembles)
    SA_ensembles_gia = compute_virtual_station_indiv(station_order, station_ensembles_gia)
    SA_ensembles_gia_grd = compute_virtual_station_indiv(station_order, station_ensembles_gia_grd)
    SA_ensembles_resvlm = compute_virtual_station_indiv(station_order, station_ensembles_resvlm)
    SA_ensembles_full = compute_virtual_station_indiv(station_order, station_ensembles_full)
    return

def compute_virtual_station_indiv(station_order,ensemble):
    global settings
    for merge_step, merge_stats in enumerate(station_order):
        combine_height = np.zeros([2, len(settings['time']),len(settings['ens_range'])])
        combine_height[0, ...] = ensemble[:, :, merge_stats[0]]
        combine_height[1, ...] = ensemble[:, :, merge_stats[1]]
        num_obs = np.sum(np.isfinite(combine_height[:,:,0]), axis=0).astype(float)
        num_obs[num_obs == 0] = np.nan
        combine_height[0,...] -= np.mean(combine_height[0,num_obs == 2,:], axis=0)[np.newaxis,:]
        combine_height[1,...] -= np.mean(combine_height[1,num_obs == 2,:], axis=0)[np.newaxis,:]
        new_height = np.nansum(combine_height, axis=0) / num_obs[:,np.newaxis]
        ensemble = np.append(ensemble, new_height[:, :, np.newaxis], axis=2)
        ensemble = np.delete(ensemble, merge_stats, axis=2)
    merged_stat = ensemble.squeeze()
    return(merged_stat)

def determine_station_order():
    global station_ensembles, station_list, settings
    lat = np.zeros(len(station_list))
    lon = np.zeros(len(station_list))
    for idx,stat in enumerate(station_list):
        lat[idx] = stat['lat']
        lon[idx] = stat['lon']
    time_acc = np.isfinite(station_ensembles[:,0,:])
    merge_order = []
    run_flag    = True
    while run_flag:
        latmid,lonmid,merged_stats,virstat_time_acc = compute_midpoint_min_overlap(lat, lon, time_acc)
        lat = np.delete(lat,merged_stats)
        lon = np.delete(lon,merged_stats)
        merge_order.append(merged_stats)
        time_acc = np.delete(time_acc,merged_stats,axis=1)
        lat = np.append(lat,latmid)
        lon = np.append(lon,lonmid)
        time_acc = np.append(time_acc,virstat_time_acc[:,np.newaxis],axis=1)
        if len(lat)==1: run_flag=False
    station_order = np.array(merge_order)
    return(station_order)

def compute_midpoint_min_overlap(lat,lon,time_acc):
    lonmat,latmat = np.meshgrid(lon,lat)
    distmat = np.arcsin(np.sqrt(np.sin(np.deg2rad(0.5 * (latmat - latmat.T))) ** 2 + np.cos(np.deg2rad(latmat.T)) * np.cos(np.deg2rad(latmat)) * np.sin(np.deg2rad(0.5 * (lonmat - lonmat.T))) ** 2))
    distmat = distmat + np.tri(N=len(distmat),k=0)*1e10
    ovl_notfound = True
    while ovl_notfound:
        merged_stats = np.array(np.unravel_index(np.argmin(distmat), distmat.shape))
        n_ovl = np.sum(time_acc[:,merged_stats[0]] & time_acc[:,merged_stats[1]])
        if n_ovl >= 20*12: ovl_notfound=False
        else: distmat[merged_stats[0], merged_stats[1]] = 1e10
     # Midpoint computation
    Bx = np.cos(np.deg2rad(lat[merged_stats[1]])) * np.cos(np.deg2rad(lon[merged_stats[1]] - lon[merged_stats[0]]))
    By = np.cos(np.deg2rad(lat[merged_stats[1]])) * np.sin(np.deg2rad(lon[merged_stats[1]] - lon[merged_stats[0]]))
    latmid = np.rad2deg(np.arctan2(np.sin(np.deg2rad(lat[merged_stats[0]])) + np.sin(np.deg2rad(lat[merged_stats[1]])), np.sqrt((np.cos(np.deg2rad(lat[merged_stats[0]])) + Bx) ** 2 + By ** 2)))
    lonmid = np.rad2deg(np.deg2rad(lon[merged_stats[0]]) + np.arctan2(By, np.cos(np.deg2rad(lat[merged_stats[0]])) + Bx))
    # Accepted points of new virtual station
    virstat_time_acc = time_acc[:,merged_stats[0]] | time_acc[:,merged_stats[1]]
    return(latmid,lonmid,merged_stats,virstat_time_acc)

def corr_gia_grd_ens(ens):
    print('   Ensemble ' + str(ens + 1))
    global station_ensembles_gia, station_ensembles_gia_grd, station_ensembles_full, settings, stat_coords
    global mask_SA,area
    # Read GIA
    fh = Dataset(settings['fn_gia_ens_rsl'],'r')
    fh.set_auto_mask(False)
    rsl_gia_ens = fh['rsl'][ens,...]
    rsl_gia_SA  = (area*rsl_gia_ens*mask_SA).sum()/(area*mask_SA).sum()
    rsl_gia_lcl = rsl_gia_ens[stat_coords[:,0],stat_coords[:,1]]
    fh.close()

    # Read GRD
    fn = settings['dir_grd_ens'] + 'grd_'+str(ens)+'.nc'
    fh = Dataset(fn,'r')
    fh.set_auto_mask(False)
    rsl_grd_ens = fh['rsl'][:]
    rsl_grd_SA  = (area*rsl_grd_ens*mask_SA).sum(axis=(1,2))/(area*mask_SA).sum()
    rsl_grd_lcl = rsl_grd_ens[:,stat_coords[:,0],stat_coords[:,1]]
    fh.close()

    # Correct ensembles for GIA
    station_ensembles_gia[:,ens,:] += (rsl_gia_SA - rsl_gia_lcl)[np.newaxis,:]*(settings['time']-settings['time'].mean())[:,np.newaxis]
    station_ensembles_gia_grd[:,ens,:] += (rsl_gia_SA - rsl_gia_lcl)[np.newaxis,:]*(settings['time']-settings['time'].mean())[:,np.newaxis]
    station_ensembles_full[:,ens,:] += (rsl_gia_SA - rsl_gia_lcl)[np.newaxis,:]*(settings['time']-settings['time'].mean())[:,np.newaxis]

    # Correct for GRD
    grd_anom_int = interp1d(np.arange(1900,2019)+0.5,rsl_grd_SA[:,np.newaxis] - rsl_grd_lcl,kind='linear',fill_value='extrapolate',axis=0)(settings['time'])
    station_ensembles_gia_grd[:,ens,:] += grd_anom_int
    station_ensembles_full[:,ens,:]    += grd_anom_int
    return

def compute_station_ensembles():
    global settings, station_list, station_ensembles
    # For each station, produce an ensemble of (uncorrected) sea-level data
    station_ensembles = mp_empty_float([len(settings['time']),len(settings['ens_range']),len(station_list)])
    print(' Reading sea-level data...')
    for idx, stat in enumerate(station_list):
        print('   '+stat['name']+'...')
        if stat['type'] == 'PSMSL': station_ensembles[:,:,idx] = compute_station_psmsl(stat)
        elif stat['name'] == 'Dakar': station_ensembles[:,:,idx] = compute_station_dakar(stat)
        elif stat['name'] == 'Kerguelen': station_ensembles[:,:,idx] = compute_station_kerguelen(stat)
        elif stat['name'] == 'Falklands': station_ensembles[:,:,idx] = compute_saltmarsh_falklands(stat)
    return()

def compute_station_psmsl(stat):
    global settings
    statlist_PSMSL = np.loadtxt(settings['fn_statlist_PSMSL'],usecols=(0,1,2),delimiter=';',dtype=object)
    statlist_PSMSL[:,0] = statlist_PSMSL[:,0].astype(int)
    statlist_PSMSL[:,1:] = statlist_PSMSL[:,1:].astype(float)
    statlist_PSMSL[:,2] = np.mod(statlist_PSMSL[:,2],360)
    n_stats = len(stat['rlr_num'])
    if 'metric_num' in stat: n_stats+=len(stat['metric_num'])
    tg_indiv = np.zeros([n_stats,len(settings['time'])])*np.nan
    for idx, statnum in enumerate(stat['rlr_num']):
        lat = statlist_PSMSL[statlist_PSMSL[:,0] == statnum,1][0]
        lon = statlist_PSMSL[statlist_PSMSL[:,0] == statnum,2][0]
        # Download file
        url = 'https://www.psmsl.org/data/obtaining/rlr.monthly.data/'+str(statnum)+'.rlrdata'
        tg_indiv[idx,:] = proc_tg_file(lat, lon, url, statnum, settings)
    if 'metric_num' in stat:
        for idx, statnum in enumerate(stat['metric_num']):
            lat = statlist_PSMSL[statlist_PSMSL[:, 0] == statnum, 1][0]
            lon = statlist_PSMSL[statlist_PSMSL[:, 0] == statnum, 2][0]
            url = 'https://www.psmsl.org/data/obtaining/met.monthly.data/'+str(statnum)+'.metdata'
            tg_indiv[idx+len(stat['rlr_num']),:] = proc_tg_file(lat, lon, url, statnum, settings)
    # Merge individual stations using virtual station method
    while tg_indiv.shape[0] > 1:
        nstats = tg_indiv.shape[0]
        ovl_matrix = np.ones([nstats, nstats], dtype=np.int)*-1
        for i in range(nstats):
            for j in range(nstats):
                if i!=j: ovl_matrix[i, j] = np.sum(np.isfinite(tg_indiv[i,:]*tg_indiv[j,:]))
        stats_to_merge = np.where(ovl_matrix == ovl_matrix.max())[0]
        combine_height = np.array([tg_indiv[stats_to_merge[0],:], tg_indiv[stats_to_merge[1],:]])
        num_obs = np.sum(np.isfinite(combine_height),axis=0) * 1.0
        num_obs[num_obs == 0] = np.nan
        combine_height[0]-= np.mean(combine_height[0, num_obs == 2])
        combine_height[1]-= np.mean(combine_height[1, num_obs == 2])
        newstat = np.nansum(combine_height, 0) / num_obs
        tg_indiv = np.delete(tg_indiv, stats_to_merge, 0)
        tg_indiv = np.append(tg_indiv, newstat[np.newaxis,:], axis=0)
    stat_indiv = np.squeeze(tg_indiv)
    # Generate random noise with GGM model
    stat_ensemble = hector.generate_ggm_noise(settings['time'], stat_indiv, len(settings['ens_range']))
    return(stat_ensemble)

def compute_station_dakar(stat):
    global settings
    dakar_raw = np.loadtxt(settings['fn_dakar'])
    time_dakar = np.round(dakar_raw[:,0] + dakar_raw[:,1]/12 - (1/24),4)
    tg_time = np.zeros(len(time_dakar))
    tg_height = dakar_raw[:,2]
    tg_height -=tg_height.mean()
    for i in range(len(time_dakar)): tg_time[i] = settings['time'][np.argmin(np.abs(settings['time'] - time_dakar[i]))]
    tg_tseries = rm_meteo_nodal(14.683333,np.mod(-17.416667,360),tg_time,tg_height,settings)
    stat_ensemble = hector.generate_ggm_noise(settings['time'], tg_tseries, len(settings['ens_range']))
    return(stat_ensemble)

def compute_station_kerguelen(stat):
    # Testut et al., 2006
    kerguelen_trend  = 1.1
    kerguelen_sterr = 0.7
    acc_idx = (settings['time']>1962) & (settings['time']<2006)
    tg_tseries = np.zeros([len(settings['time']),len(settings['ens_range'])]) * np.nan
    tg_tseries[acc_idx,:] = np.random.normal(loc=kerguelen_trend,scale=kerguelen_sterr,size=len(settings['ens_range']))[np.newaxis,:] * (settings['time'][acc_idx] - settings['time'][acc_idx].mean())[:,np.newaxis]
    return(tg_tseries)

def compute_saltmarsh_falklands(stat):
    global settings
    # -------------------------------------------------------------------------------------
    # Read Falklands Salt Marsh data:
    # Two possible runs
    # 1. Just Swan inlet salt marsh data, but with all the points
    # 2. Swan inlet, except points without excellent modern analogues, but with Port Louis.
    # -------------------------------------------------------------------------------------

    # Read Swan Inlet data
    tg_tseries = np.zeros([len(settings['time']),len(settings['ens_range'])]) * np.nan
    falklands_raw = np.loadtxt(settings['fn_falklands'],skiprows=1,usecols=(1,2,3,4,5))
    swan_time       = falklands_raw[:,0]
    swan_height     = 1000*falklands_raw[:,3]
    swan_time_ste   = (falklands_raw[:,1] - falklands_raw[:,2])/4
    swan_height_ste = 1000*(falklands_raw[:,4])/2

    # Read PSMSL Stanley 1,2
    psmsl_stanley_1 = np.loadtxt(settings["dir_data"]+'TideGauges/rlr_monthly/data/1082.rlrdata',delimiter=';')
    acc_idx = (((psmsl_stanley_1[:, 3] == 10) | (psmsl_stanley_1[:, 3] == 0)) & (psmsl_stanley_1[:, 1] > -10000))
    psmsl_stanley_1 = psmsl_stanley_1[acc_idx,:2]

    psmsl_stanley_2 =np.loadtxt(settings["dir_data"]+'TideGauges/rlr_monthly/data/1796.rlrdata',delimiter=';')
    acc_idx = (((psmsl_stanley_2[:, 3] == 10) | (psmsl_stanley_2[:, 3] == 0)) & (psmsl_stanley_2[:, 1] > -10000))
    psmsl_stanley_2 = psmsl_stanley_2[acc_idx,:2]

    # Read PLW data: Stanley and define Port Louis
    tg_louis = np.zeros([4, 3])
    tg_louis[0, :] = [1842.5,-1752.855,41]
    tg_louis[1, :] = [1982.0,-1696.2,28]
    tg_louis[2, :] = [1984.5,-1745.5  ,35]
    tg_louis[3, :] = [2009.17,-1622.9 ,73]
    tg_govjetty_plw = np.loadtxt(settings['dir_sl_plw']+'slc.to.rossbm.govjetty.ac.sl')
    tg_govjetty_plw[:,1] *= 10
    tg_fipass_plw = np.loadtxt(settings['dir_sl_plw']+'slc.to.rossbm.fipass.ac.sl')
    tg_fipass_plw[:,1] *= 10

    # Correct PSMSL data to align with Ross benchmark
    rlr2ross = psmsl_stanley_1[:,1].mean() - tg_govjetty_plw[:,1].mean()
    psmsl_stanley_1[:,1] -= rlr2ross
    psmsl_stanley_2[:,1] -= rlr2ross

    # Correct Swan inlet to tie the 2006 index point to Ross BM
    swan_height += (psmsl_stanley_2[np.floor(psmsl_stanley_2[:,0]) == 2006,1].mean() - swan_height[0])

    plot=False
    if plot:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(8,4))
        plt.plot(tg_govjetty_plw[:,0],tg_govjetty_plw[:,1],'k')
        plt.plot(tg_fipass_plw[:,0],tg_fipass_plw[:,1],'k')
        plt.plot(psmsl_stanley_2[:,0],psmsl_stanley_2[:,1],'C0')

        plt.plot(tg_louis[:,0],tg_louis[:,1],'o',color='C1')
        plt.plot(swan_time,swan_height,'s',color='C2')

        plt.ylabel('Height relative to Ross BM (cm)',fontsize=9)
        plt.xticks(fontsize=9)
        plt.yticks(fontsize=9)
        plt.grid()
        plt.tight_layout()
        plt.savefig('FL_tg.png')

    if settings['falklands_test']:
        # Remove problematic idx points and insert Port Louis
        acc_idx = np.ones(len(swan_time), dtype=np.bool)
        acc_idx[1:6] = False
        acc_idx[-2] = False
        swan_time = swan_time[acc_idx]
        swan_height = swan_height[acc_idx]
        swan_time_ste = swan_time_ste[acc_idx]
        swan_height_ste = swan_height_ste[acc_idx]
        # Pt Louis
        swan_time = np.append(swan_time,tg_louis[1:,0])
        swan_time_ste = np.append(swan_time_ste,np.zeros(len(tg_louis[1:,0])))
        swan_height = np.append(swan_height,tg_louis[1:,1])
        swan_height_ste = np.append(swan_height_ste,tg_louis[1:,2])

    height_rnd = swan_height[:,np.newaxis] + swan_height_ste[:,np.newaxis] * np.random.normal(loc=1,scale=1,size=([len(swan_time),len(settings['ens_range'])]))
    time_rnd   = swan_time[:,np.newaxis] + swan_time_ste[:,np.newaxis] *np.random.normal(loc=1,scale=1,size=([len(swan_time),len(settings['ens_range'])]))
    for ens in settings['ens_range']:
        sort_idx = np.argsort(time_rnd[:,ens])
        height_rnd[:,ens] = height_rnd[sort_idx,ens]
        time_rnd[:,ens] = time_rnd[sort_idx,ens]

    for ens in settings['ens_range']:
        tg_tseries[:,ens] = interp1d(time_rnd[:,ens],height_rnd[:,ens],kind='linear',bounds_error=False,fill_value=(height_rnd[0,ens],height_rnd[-1,ens]))(settings['time'])
    return(tg_tseries)

def proc_tg_file(lat,lon,url,statnum,settings):
    # Download file
    fn_tg = settings['dir_scratch'] + str(statnum) + 'txt'
    null = urllib.request.urlretrieve(url, fn_tg)
    # Read file
    tgdata = np.loadtxt(fn_tg, delimiter=';')
    # Flagged values out
    acc_idx = (((tgdata[:, 3] == 10) | (tgdata[:, 3] == 0)) & (tgdata[:, 1] > -10000)) & (tgdata[:, 0] >= settings['time'][0] - 0.01) & (tgdata[:, 0] <= settings['time'][-1] + 0.01)
    tgdata = tgdata[acc_idx, :]

    # Interpolate data on common time steps
    tg_time = np.zeros(len(tgdata))
    for i in range(len(tgdata)): tg_time[i] = settings['time'][np.argmin(np.abs(settings['time'] - tgdata[i, 0]))]
    tg_height = tgdata[:, 1]
    tg_height -= tg_height.mean()
    os.remove(fn_tg)
    tg_tseries = rm_meteo_nodal(lat,lon,tg_time,tg_height,settings)
    return(tg_tseries)

def rm_meteo_nodal(lat,lon,tg_time,tg_height,settings):
    # Remove meteo and nodal cycle
    fh = Dataset(settings['fn_nodal_amp'], 'r')
    fh.set_auto_mask(False)
    lat_idx_nod = np.argmin(np.abs(lat - fh['y'][:]))
    lon_idx_nod = np.argmin(np.abs(lon - fh['x'][:]))
    nodal_cycle = -1.0 * fh['rsl_eq'][lat_idx_nod, lon_idx_nod] * np.cos(2 * np.pi / 18.612958 * (tg_time - 1922.7))
    fh.close()
    tg_height -= nodal_cycle

    fh = Dataset(settings['fn_ERA_monthly'], 'r')
    fh.set_auto_mask(False)
    lat_idx_era = np.argmin(np.abs(lat - fh['lat'][:]))
    lon_idx_era = np.argmin(np.abs(lon - fh['lon'][:]))
    mslp = fh['mslp'][:, lat_idx_era, lon_idx_era]
    uws = fh['uws'][:, lat_idx_era, lon_idx_era]
    vws = fh['vws'][:, lat_idx_era, lon_idx_era]
    fh.close()

    time_acc = np.in1d(settings['time'], tg_time)
    amat = np.ones([len(tg_time), 9])
    amat[:, 1] = tg_time - tg_time.mean()
    amat[:, 2] = np.sin(2 * np.pi * tg_time)
    amat[:, 3] = np.cos(2 * np.pi * tg_time)
    amat[:, 4] = np.sin(4 * np.pi * tg_time)
    amat[:, 5] = np.cos(4 * np.pi * tg_time)
    amat[:, 6] = mslp[time_acc] - mslp[time_acc].mean()
    amat[:, 7] = uws[time_acc] - uws[time_acc].mean()
    amat[:, 8] = vws[time_acc] - vws[time_acc].mean()
    sol = np.linalg.lstsq(amat, tg_height, rcond=None)[0]
    sol[:2] = 0
    tg_height -= (amat @ sol)
    tg_tseries = np.zeros(len(settings['time']))*np.nan
    tg_tseries[time_acc] = tg_height
    return(tg_tseries)

def save_ensembles():
    global station_ensembles, station_ensembles_gia, station_ensembles_gia_grd, station_ensembles_resvlm, station_ensembles_full, settings, station_list
    global SA_ensembles, SA_ensembles_gia, SA_ensembles_gia_grd, SA_ensembles_resvlm, SA_ensembles_full
    np.savez_compressed(settings['fn_sl_ensembles_stat'],station_ensembles=station_ensembles,station_ensembles_gia=station_ensembles_gia,station_ensembles_gia_grd=station_ensembles_gia_grd,station_ensembles_resvlm=station_ensembles_resvlm,station_ensembles_full=station_ensembles_full,station_list=station_list)
    np.savez_compressed(settings['fn_sl_ensembles_SA'],SA_ensembles=SA_ensembles,SA_ensembles_gia=SA_ensembles_gia,SA_ensembles_gia_grd=SA_ensembles_gia_grd,SA_ensembles_resvlm=SA_ensembles_resvlm,SA_ensembles_full=SA_ensembles_full)
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



