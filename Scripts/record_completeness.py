# Calculate record completeness for Table 1
import os
import numpy as np

def main():
    settings = {}
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['fn_station_list'] = settings['dir_project'] + 'Data/station_list.npy'
    settings['fn_sl_stats_stat'] = settings['dir_project'] + 'Data/sl_stats_stat.npy'
    settings['time'] = np.arange(1900+1/24,2019+1/24,1/12)

    sl_stats_stat = np.load(settings['fn_sl_stats_stat'],allow_pickle=True)
    station_list = np.load(settings['fn_station_list'],allow_pickle=True)

    for idx, stat in enumerate(station_list):
        record_completeness(stat['name'], sl_stats_stat[idx]['no_corr']['tseries'], settings)
    return

def record_completeness(statname,tseries,settings):
    has_vals = np.where(np.isfinite(tseries))[0]
    t_tot  = len(settings['time'][has_vals[0]:has_vals[-1]+1])
    frac   = len(has_vals)/t_tot
    ystart = np.floor(settings['time'][has_vals[0]])
    ystop  = np.floor(settings['time'][has_vals[-1]])
    print(statname)
    print(ystart)
    print(ystop)
    print(frac)
    print('------')
    return