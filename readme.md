# Code and Data supplement for "Constraining 20th-century sea-level rise in the South Atlantic Ocean" 

This repository contains the scripts and data generated for the manuscript  
"Constraining 20th-century sea-level rise in the South Atlantic Ocean"   
Thomas Frederikse<sup>1</sup>, Surendra Adhikari<sup>1</sup>, Tim J. Daley<sup>2</sup>, Soenke Dangendorf<sup>3</sup>, Roland Gehrels<sup>4</sup>, Felix Landerer<sup>1</sup>, Marta Marcos<sup>5,6</sup>, Thomas L. Newton<sup>2</sup>, Graham Rush<sup>4</sup>, Aimee B.A. Slangen<sup>7</sup>, Guy Woeppelmann <sup>8</sup>

<sup>1</sup> Jet Propulsion Laboratory, California Institute of Technology, Pasadena, California, USA  
<sup>2</sup> School of Geography, Earth and Environmental Sciences, Plymouth University, Plymouth,  UK  
<sup>3</sup> Old Dominion University, Norfolk, Virginia, USA & University of Siegen, Siegen, Germany  
<sup>4</sup> Department of Environment and Geography, University of York, Heslington, York, UK  
<sup>5</sup> IMEDEA (UIB-CSIC), Esporles, Spain  
<sup>6</sup> Department of Physics, University of the Balearic Islands, Palma, Spain  
<sup>7</sup> NIOZ Royal Netherlands Institute for Sea Research, department of Estuarine and Delta Systems, and Utrecht University, Yerseke, The Netherlands  
<sup>8</sup> LIENSs, Universite de La Rochelle - CNRS, La Rochelle, France  

© 2021 All rights reserved

Part of this work was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (80NM0018D0004)

## Directory `Data`

### GRD_rsl_uniform.grd:
netCDF file with uniform realtive sea-level fingerprints for Greenland, West- and East-Antarctica and glaciers. All fingerprints are normalized to a unit barystatic contribution. Multiple free software packages are available to view and modify NetCDF files. For example [python](https://unidata.github.io/netcdf4-python/), [Julia](https://github.com/Alexander-Barth/NCDatasets.jl), [GMT](https://www.generic-mapping-tools.org/), [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html), and many others. 

### steric_trends_1957_2018.nc:
netCDF file with linear trends over 1957-2018 in steric sea level from various data sets. All units are mm. When using the data in this file, please properly cite the underlying data source. For each of the solutions, the appropriate reference is included in the metadata of the netCDF file and can also be found in the manuscript. Multiple free software packages are available to view and modify NetCDF files. For example [python](https://unidata.github.io/netcdf4-python/), [Julia](https://github.com/Alexander-Barth/NCDatasets.jl), [GMT](https://www.generic-mapping-tools.org/), [ncview](http://meteora.ucsd.edu/~pierce/ncview_home_page.html), and many others. 

### Swan_inlet_sl_data.txt
Text file with the sea-level index points from the Swan Inlet salt marsh. Please cite Newton, T. L., Gehrels, W. R., Fyfe, R. M., & Daley, T. J. (2020). Reconstructing Sea-level change in the Falkland Islands (Islas Malvinas) using salt-marsh foraminifera, diatoms and testate amoebae. Marine Micropaleontology, 101923. [doi:10.1016/j.marmicro.2020.101923](https://doi.org/10.1016/j.marmicro.2020.101923) when using this dataset.

## Directory `Scripts`
This directory contains the Python scripts used to compute all the results of the paper. Please note that these scripts depend on some 3rd-party codes, such as [MIDAS](http://geodesy.unr.edu/) and [Hector](http://segal.ubi.pt/hector/), which have to be installed manually, as well as external data sets. Please see the paper and references herein for details on these data sets. 
