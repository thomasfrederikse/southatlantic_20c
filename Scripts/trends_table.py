# Write trends to LaTeX table
import os
import numpy as np

def main():
    settings = {}
    settings['dir_data'] = os.getenv('HOME') + '/Data/'
    settings['dir_project'] = os.getenv('HOME') + '/Projects/2020_SA/'
    settings['fn_station_list'] = settings['dir_project'] + 'Data/station_list.npy'
    settings['fn_sl_stats_stat'] = settings['dir_project'] + 'Data/sl_stats_stat.npy'
    settings['fn_sl_stats_SA'] = settings['dir_project'] + 'Data/sl_stats_SA.npy'

    sl_stats_stat = np.load(settings['fn_sl_stats_stat'], allow_pickle=True)
    sl_stats_SA = np.load(settings['fn_sl_stats_SA'], allow_pickle=True).all()

    station_list = np.load(settings['fn_station_list'], allow_pickle=True)
    print_table = []
    # Rij 1 GMSL
    for stat in sl_stats_stat:
        print_table.append(stat['name']+'&'+print_sl_trend(stat['no_corr']['trend'])+'&'+print_sl_trend(stat['gia']['trend'])+'&'+print_sl_trend(stat['gia_grd']['trend'])+'&'+print_sl_trend(stat['resvlm']['trend'])+'&'+print_sl_trend(stat['full']['trend']) +'\\\\')
    print_table.append("\hline")
    print_table.append('South Atlantic &' + print_sl_trend(sl_stats_SA['no_corr']['trend']) + '&' + print_sl_trend(sl_stats_SA['gia']['trend']) + '&' + print_sl_trend(sl_stats_SA['gia_grd']['trend']) + '&' + print_sl_trend(sl_stats_SA['resvlm']['trend']) + '&' + print_sl_trend(sl_stats_SA['full']['trend']) + '\\\\')
    for i in print_table:
        print(i)

    return



def print_sl_trend(trend):
    trend_mn = "{:.2f}".format(trend[1])
    trend_lo = "{:.2f}".format(trend[0])
    trend_hi = "{:.2f}".format(trend[2])
    print_str = trend_mn + ' &['+trend_lo+' '+trend_hi+']'
    return(print_str)
