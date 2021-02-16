# ------------------
# Write station list
# ------------------
import numpy as np
import os

def main():
    global settings
    settings = {}
    settings['dir_project'] = os.getenv('HOME')+'/Projects/2020_SA/'
    settings['fn_station_list'] = settings['dir_project'] + 'Data/station_list.npy'
    station_list = np.zeros(8,dtype=dict)
    for n in range(len(station_list)): station_list[n] = {}
    n=0
    station_list[n]['name'] = 'Buenos Aires'
    station_list[n]['type'] = 'PSMSL'
    station_list[n]['rlr_num'] = [157, 832]
    station_list[n]['gps_list'] = ['IGM1','BUE1','BUE2','MA02','LPGS']
    station_list[n]['lat'] = -34.6
    station_list[n]['lon'] = np.mod(-58.366667,360)

    n=1
    station_list[n]['name'] = 'Montevideo'
    station_list[n]['type'] = 'PSMSL'
    station_list[n]['rlr_num'] = [431, 434, 764]
    station_list[n]['gps_list'] = ['UYMO','MTV1','MTV2','UYLP']
    station_list[n]['lat'] = -34.9
    station_list[n]['lon'] = np.mod(-56.25,360)

    n=2
    station_list[n]['name'] = 'Mar del Plata'
    station_list[n]['type'] = 'PSMSL'
    station_list[n]['rlr_num']  = [223, 819, 857]
    station_list[n]['metric_num'] = [177]
    station_list[n]['gps_list'] = ['BCAR','MPL2','MPLA']
    station_list[n]['lat'] = -38.033333
    station_list[n]['lon'] = np.mod(-57.516667,360)

    n=3
    station_list[n]['name'] = 'Puerto Madryn'
    station_list[n]['type'] = 'PSMSL'
    station_list[n]['rlr_num'] = [501, 867]
    station_list[n]['gps_list'] = ['RWSN']
    station_list[n]['lat'] = -42.766667
    station_list[n]['lon'] = np.mod(-65.033333,360)

    n=4
    station_list[n]['name'] = 'Dakar'
    station_list[n]['type'] = 'Custom'
    station_list[n]['gps_list'] = ['DAKR','DAKA','FG02']
    station_list[n]['lat'] = 14.683333
    station_list[n]['lon'] = np.mod(-17.416667,360)

    n=5
    station_list[n]['name'] = 'South Africa'
    station_list[n]['type'] = 'PSMSL'
    station_list[n]['rlr_num'] = [284, 820, 826, 836, 910, 911, 914, 950, 1192, 1195]
    station_list[n]['gps_list'] = ['ANTH', 'BENI', 'BETH', 'BFTN', 'BISO', 'BRIT', 'BRNK', 'BWES','CALV', 'CPNT', 'CTWN', 'DEAR', 'DRBA', 'DRBN', 'ELDN', 'EMLO','ERAS', 'FG08', 'GDAL', 'GEO1', 'GEOA', 'GREY', 'GRHM', 'GRNT','HARB', 'HEID', 'HNUS', 'HRAC', 'HRAO', 'IXOP', 'KLEY','KMAN', 'KRUG', 'KSTD', 'LGBN', 'LSMH', 'MALM', 'MBRG', 'MFKG','MRIV', 'NSPT', 'NYLS', 'OKNY', 'PBWA', 'PELB', 'PMBG', 'POTG','PRE1', 'PRET', 'PSKA', 'PTBG', 'QTWN', 'RBAY', 'SBOK','SCO1', 'SPRT', 'STBS', 'STNG', 'SUT1', 'SUTH', 'SUTM', 'SUTV','TDOU', 'ULDI', 'UMTA', 'UPTA', 'VERG', 'WIND', 'WORC']
    station_list[n]['lat'] = -33.905278
    station_list[n]['lon'] = 18.434722

    n=6
    station_list[n]['name'] = 'Kerguelen'
    station_list[n]['type'] = 'Custom'
    station_list[n]['gps_list']  = ['KERG','KRGG','KETG']
    station_list[n]['DORIS_list'] = ['KETB']
    station_list[n]['lat'] = -49.351501
    station_list[n]['lon'] = 70.255501

    n=7
    station_list[n]['name'] = 'Falklands'
    station_list[n]['type'] = 'Custom'
    station_list[n]['gps_list'] = ['FALK','LKTH']
    station_list[n]['lat'] = -51.692358
    station_list[n]['lon'] = np.mod(-57.820514,360)


    np.save(settings['fn_station_list'],station_list)
    return

if __name__ == '__main__':
    main()
