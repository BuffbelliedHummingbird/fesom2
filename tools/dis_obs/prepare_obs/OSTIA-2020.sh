#!/bin/bash

# Download from website:
# https://data.marine.copernicus.eu/product/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/files?subdataset=METOFFICE-GLO-SST-L4-REP-OBS-SST_202003&path=SST_GLO_SST_L4_REP_OBSERVATIONS_010_011%2FMETOFFICE-GLO-SST-L4-REP-OBS-SST_202003%2F

# Move from Downloads to Albedo:
rsync -rav ...

# combine all in directory to 1-year-file:
cdo mergetime *.nc OSTIA_SST_20200101_20201231_regular_daily.nc

# regrid:
cdo remapbil,/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/core2_griddes_nodes.nc OSTIA_SST_20020101_20201231_regular_daily.nc OSTIA_SST_20020101_20201231_daily.nc

# rename dimensions:
ncrename -d ncells,nodes_2d  OSTIA_SST_20020101_20201231_daily.nc
ncrename -d time,T  OSTIA_SST_20020101_20201231_daily.nc

# rename variables:
ncrename -v analysed_sst,sst  OSTIA_SST_20020101_20201231_daily.nc
ncrename -v analysis_error,std  OSTIA_SST_20020101_20201231_daily.nc

# convert data type:
ncap2 -O -s 'sst=float(sst)' OSTIA_SST_20200101_20201231_daily.nc OSTIA_SST_20200101_20201231_daily_2.nc
ncap2 -O -s 'std=float(std)' OSTIA_SST_20200101_20201231_daily_2.nc OSTIA_SST_20200101_20201231_daily_3.nc
ncap2 -O -s 'time=float(time)' OSTIA_SST_20200101_20201231_daily_3.nc OSTIA_SST_20200101_20201231_daily_4.nc

# from Kelvin to Celsius:
ncap2 -s 'sst=sst-273.15f' OSTIA_SST_20200101_20201231_daily_4.nc OSTIA_SST_20200101_20201231_daily_5.nc

# change missing value:
ncatted -a _FillValue,sst,m,d,-999 OSTIA_SST_20200101_20201231_daily_5.nc
ncatted -a _FillValue,std,m,d,-999 OSTIA_SST_20200101_20201231_daily_5.nc
