# ---------------------------------------
# Compute statistics from ensemble of
# station and basin-mean sea-level curves
# ---------------------------------------
import numpy as np
import os
import mod_gentools as gentools
from netCDF4 import Dataset
import mod_hector as hector
from scipy import signal
def main():
    global settings
    set_settings()
    load_ensembles()
    read_gia_likelihood()
    station_stats, SA_stats = compute_stats()
    np.save(settings['fn_sl_stats_stat'],station_stats)
    np.save(settings['fn_sl_stats_SA'],SA_stats)
    return

def set_settings():
    global settings
    settings = {}
    settings['dir_data']    = os.getenv('HOME')+'/Data/'
    settings['dir_scratch'] = os.getenv('HOME')+'/Scratch/'
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['dir_grd_ens'] = os.getenv('HOME')+'/Data/Budget_20c/grd/'
    settings['dir_sl_ensembles'] = settings['dir_project']+'/Data/sl_ensembles/'

    settings['fn_mask']      = settings['dir_data'] + 'GRACE/JPL_mascon/mask.npy'
    settings['fn_gia_ens_rsl'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rsl_ens_05.nc'

    settings['ens_range'] = np.arange(100)
    settings['falklands_test'] = False
    settings['orig_test'] = False

    settings['fn_station_list'] = settings['dir_project'] + 'Data/station_list.npy'
    settings['fn_residual_vlm'] = settings['dir_project'] + 'Data/residual_vlm.npy'
    settings['fn_gia_ens_rad'] = settings['dir_data'] + 'GIA/Caron/Ensemble/rad_ens_05.nc'

    settings['fn_statlist_PSMSL'] = settings['dir_project'] + 'Data/sealevel/statlist_PSMSL.txt'
    settings['fn_ERA_monthly'] = settings['dir_data'] + 'SA/ERA_monthly.nc'
    settings['fn_nodal_amp'] = settings['dir_data'] + 'Nodal/Nodal.nc'
    settings['fn_dakar']       = settings['dir_data'] + 'TideGauges/Dakar/Dakar_monthly_completedwithPSMSL_Thomas.dat'
    settings['fn_falklands']   = settings['dir_project'] + 'Data/sealevel/Falklands_2020_10.txt'

    if settings['falklands_test']:
        settings['fn_sl_ensembles_stat'] = settings['dir_sl_ensembles'] + 'sl_ensembles_stat_falklands_test.npz'
        settings['fn_sl_ensembles_SA'] = settings['dir_sl_ensembles'] + 'sl_ensembles_SA_falklands_test.npz'
        settings['fn_sl_stats_stat'] = settings['dir_project'] + 'Data/sl_stats_stat_falklands_test.npy'
        settings['fn_sl_stats_SA'] = settings['dir_project'] + 'Data/sl_stats_SA_falklands_test.npy'
    elif settings['orig_test']:
        settings['fn_sl_ensembles_stat'] = settings['dir_sl_ensembles'] + 'sl_ensembles_stat_orig_test.npz'
        settings['fn_sl_ensembles_SA'] = settings['dir_sl_ensembles'] + 'sl_ensembles_SA_orig_test.npz'
        settings['fn_sl_stats_stat'] = settings['dir_project'] + 'Data/sl_stats_stat_orig_test.npy'
        settings['fn_sl_stats_SA'] = settings['dir_project'] + 'Data/sl_stats_SA_orig_test.npy'
    else:
        settings['fn_sl_ensembles_stat'] = settings['dir_sl_ensembles'] + 'sl_ensembles_stat.npz'
        settings['fn_sl_ensembles_SA'] = settings['dir_sl_ensembles'] + 'sl_ensembles_SA.npz'
        settings['fn_sl_stats_stat'] = settings['dir_project'] + 'Data/sl_stats_stat.npy'
        settings['fn_sl_stats_SA'] = settings['dir_project'] + 'Data/sl_stats_SA.npy'


    settings['time'] = np.arange(1900+1/24,2019+1/24,1/12)
    return

def load_ensembles():
    global station_ensembles, station_ensembles_gia, station_ensembles_gia_grd, station_ensembles_resvlm, station_ensembles_full, settings, station_list
    global SA_ensembles, SA_ensembles_gia, SA_ensembles_gia_grd, SA_ensembles_resvlm, SA_ensembles_full

    stat_ens = np.load(settings['fn_sl_ensembles_stat'],allow_pickle=True)
    station_ensembles         = stat_ens['station_ensembles']
    station_ensembles_gia     = stat_ens['station_ensembles_gia']
    station_ensembles_gia_grd = stat_ens['station_ensembles_gia_grd']
    station_ensembles_resvlm  = stat_ens['station_ensembles_resvlm']
    station_ensembles_full    = stat_ens['station_ensembles_full']
    station_list    = stat_ens['station_list']

    del stat_ens
    SA_ens = np.load(settings['fn_sl_ensembles_SA'],allow_pickle=True)
    SA_ensembles         = SA_ens['SA_ensembles']
    SA_ensembles_gia     = SA_ens['SA_ensembles_gia']
    SA_ensembles_gia_grd = SA_ens['SA_ensembles_gia_grd']
    SA_ensembles_resvlm  = SA_ens['SA_ensembles_resvlm']
    SA_ensembles_full    = SA_ens['SA_ensembles_full']
    return

def read_gia_likelihood():
    global likelihood
    likelihood = Dataset(settings['fn_gia_ens_rad'],'r').variables['probability'][settings['ens_range']]._get_data()
    likelihood /= likelihood.sum()
    return

def compute_stats():
    global station_ensembles, station_ensembles_gia, station_ensembles_gia_grd, station_ensembles_resvlm, station_ensembles_full, settings, station_list
    global SA_ensembles, SA_ensembles_gia, SA_ensembles_gia_grd, SA_ensembles_resvlm, SA_ensembles_full

    # For each station
    station_stats = np.zeros(len(station_list),dtype=dict)
    for idx, stat in enumerate(station_list):
        station_stats[idx] = {}
        station_stats[idx]['name'] = stat['name']
        station_stats[idx]['no_corr'] = compute_indiv_stat(station_ensembles[:,:,idx], do_global=False)
        station_stats[idx]['gia'] = compute_indiv_stat(station_ensembles_gia[:,:,idx], do_global=False)
        station_stats[idx]['gia_grd'] = compute_indiv_stat(station_ensembles_gia_grd[:,:,idx], do_global=False)
        station_stats[idx]['resvlm'] = compute_indiv_stat(station_ensembles_resvlm[:,:,idx], do_global=False)
        station_stats[idx]['full'] = compute_indiv_stat(station_ensembles_full[:,:,idx], do_global=False)

    # SA mean
    SA_stats = {}
    SA_stats['no_corr'] = compute_indiv_stat(SA_ensembles, do_global=True)
    SA_stats['gia'] = compute_indiv_stat(SA_ensembles_gia, do_global=True)
    SA_stats['gia_grd'] = compute_indiv_stat(SA_ensembles_gia_grd, do_global=True)
    SA_stats['resvlm'] = compute_indiv_stat(SA_ensembles_resvlm, do_global=True)
    SA_stats['full'] = compute_indiv_stat(SA_ensembles_full, do_global=True)
    return(station_stats, SA_stats)


def compute_indiv_stat(ensemble, do_global):
    global settings, likelihood
    stats = {}
    stats['tseries'] = np.nansum(ensemble * likelihood, axis=1)
    if do_global:
        # Compute trend over 1901-2010
        t_acc = (settings['time'] > 1901) & (settings['time'] < 2011) & (np.isfinite(ensemble[:,0]))
    else:
        t_acc = np.isfinite(ensemble[:,0])
        stats['tseries'][~t_acc] = np.nan
    # Trend from Hector
    time = settings['time'][t_acc]
    ggm_trend = hector.est_trend(time,stats['tseries'][t_acc],NoiseModel='GGM',AR=1)['trend']

    # Trend from ensemble
    amat = np.ones([len(time), 2])
    amat[:, 1] = time - time.mean()
    amat_T = amat.T
    amat_sq = np.linalg.inv(np.dot(amat_T, amat))

    trend_array = np.zeros(len(settings['ens_range']))
    # Loop over all ensemble members
    for i in range(len(settings['ens_range'])):
        sol = np.dot(amat_sq, np.dot(amat_T, ensemble[t_acc,i]))
        trend_array[i] = sol[1]

    # Statistics from ensemble
    sort_idx = np.argsort(trend_array)
    sort_cdf = np.cumsum((likelihood / likelihood.sum())[sort_idx])
    stats['trend'] = np.zeros(3)
    stats['trend'][0] = np.abs(trend_array[sort_idx][np.argmin(np.abs(sort_cdf - 0.05))])
    stats['trend'][1] = ggm_trend[0]
    stats['trend'][2] = np.abs(trend_array[sort_idx][np.argmin(np.abs(sort_cdf - 0.95))])
    # stats['trend'][0] = ggm_trend[0] - stats['trend'][0]
    # stats['trend'][2] = ggm_trend[0] + stats['trend'][2]

    # Trend PDF
    n_steps = 100
    stats['trend_pdf'] = np.zeros([n_steps,2])
    pdf_boundaries = np.linspace(-2,5,n_steps+1)
    stats['trend_pdf'][:,0] = (pdf_boundaries[:n_steps]+pdf_boundaries[1:])/2
    for i in range(n_steps):
        acc_idx = (trend_array>=pdf_boundaries[i]) & (trend_array<pdf_boundaries[i+1])
        stats['trend_pdf'][i,1] = np.sum(likelihood[acc_idx])

    # Smoothing test
    b,a = signal.butter(2, 2/10)
    stats['trend_pdf'][:,1] = signal.filtfilt(b,a,np.hstack([np.zeros(100),stats['trend_pdf'][:,1],np.zeros(100)]),padtype=None,method='gust')[100:-100]
    stats['trend_pdf'][:,1] = np.fmax(stats['trend_pdf'][:,1],np.zeros(len(stats['trend_pdf'][:,1])))
    return(stats)

# if __name__ == '__main__':
#     main()
