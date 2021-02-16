# Code and Data supplement for "Constraining 20th-century sea-level rise in the South Atlantic Ocean" 

This repository contains the scripts and data generated for the manuscript "Constraining 20th-century sea-level rise in the South Atlantic Ocean" from by Thomas Frederikse [^1], Surendra Adhikari [^1], Tim J. Daley [^2], Soenke Dangendorf [^3], Roland Gehrels [^4], Felix Landerer [^1], Marta Marcos [^5,6], Thomas L. Newton [^2], Graham Rush [^4], Aimee B.A. Slangen [^7], Guy Woeppelmann [^8]

1 Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California, USA
2 School of Geography, Earth and Environmental Sciences, Plymouth University, Plymouth,  UK
3 Old Dominion University, Norfolk, Virginia, USA & University of Siegen, Siegen, Germany
4 Department of Environment and Geography, University of York, Heslington, York, UK
5 IMEDEA (UIB-CSIC), Esporles, Spain
6 Department of Physics, University of the Balearic Islands, Palma, Spain
7 NIOZ Royal Netherlands Institute for Sea Research, department of Estuarine and Delta Systems, and Utrecht University, Yerseke, The Netherlands
8 LIENSs, Universite de La Rochelle - CNRS, La Rochelle, France

© 2021 All rights reserved

Part of this work was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004)

## Directory `Data`

### GRD_rsl_uniform.grd:
netCDF file with uniform realtive sea-level fingerprints for Greenland, West- and East-Antarctica and glaciers. All fingerprints are normalized to a unit barystatic contribution. Multiple free software packages is available to view and modify NetCDF files. For example [python](https://unidata.github.io/netcdf4-python/), [Julia](https://github.com/Alexander-Barth/NCDatasets.jl), [GMT](https://www.generic-mapping-tools.org/), [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html), and many others. 

### steric_trends_1957_2018.nc:
netCDF file with linear trends over 1957-2018 in steric sea level from various data sets. All units are mm. Please properly cite the underlying data source, see the metadata of the netCDF file and/or the manuscript.  Multiple free software packages is available to view and modify NetCDF files. For example [python](https://unidata.github.io/netcdf4-python/), [Julia](https://github.com/Alexander-Barth/NCDatasets.jl), [GMT](https://www.generic-mapping-tools.org/), [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html), and many others. 

### Swan_inlet_sl_data.txt
Text file with the sea-level index points from the Swan Inlet salt marsh. Please cite Newton, T. L., Gehrels, W. R., Fyfe, R. M., & Daley, T. J. (2020). Reconstructing Sea-level change in the Falkland Islands (Islas Malvinas) using salt-marsh foraminifera, diatoms and testate amoebae. Marine Micropaleontology, 101923. [doi:10.1016/j.marmicro.2020.101923](https://doi.org/10.1016/j.marmicro.2020.101923) when using this dataset.

## Directory `Scripts`
This directory contains the Python scripts used to compute all the results of the paper. Please note that these scripts depend on some 3rd-party codes, such as [MIDAS](http://geodesy.unr.edu/) and [Hector](http://segal.ubi.pt/hector/), which have to be installed manually, as well as external data sets. Please see the paper and references herein for details on these data sets. 
