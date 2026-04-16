#!/bin/bash

# Download via ftp:
# Link: ftp://oc-cci-data:ELaiWai8ae@ftp.rsg.pml.ac.uk
# ftp ftp.rsg.pml.ac.uk
# name? --> oc-cci-data
# pw? --> ELaiWai8ae

# navigate, and set ???? to year
cd occci-v6.0/geographic/netcdf/daily/chlor_a/????
# to download
prompt noprompt # to avoid asking for each file
mget FILENAME

for((year=2010;year<=2020;year++))
do
wget -r -nH --cut-dirs=6 -nc ftp://oc-cci-data:ELaiWai8ae@ftp.rsg.pml.ac.uk//occci-v6.0/geographic/netcdf/daily/chlor_a/${year}
done
# the -nc is for only updating new files (instead of overwriting existing)

# interpolation to FESOM grid
gridfile=/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/core2_griddes_nodes.nc

# interpolation of monthly files
for((yy=2010;yy<=2020;yy++))
do
  yystr=`printf %04d $yy`
  for((mm=1;mm<=12;mm++))
  do
    mmstr=`printf %02d $mm`
    cdo remapbil,${gridfile} original/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_4km_GEO_PML_OCx-${yystr}${mmstr}-fv6.0.nc CORE2/ESACCI-OC-L3S-CHLOR_A-MERGED-1M_MONTHLY_CORE2-${yystr}${mmstr}-fv6.0.nc
  done
done

# interpolation of daily files
# using gnudate ("date" command on linux, "gdate" on Mac)

for((year=2010;year<=2020;year++)); do
   mkdir /albedo/work/projects/p_recompdaf/frbunsen/data/BGC/chlsurf_OC-CCI/CORE2_daily/${year}
   let ndays=($(date +%s -d ${year}1231)-$(date +%s -d ${year}0101))/86400+1
   for((i=1;i<=${ndays};i++))
      do datestr=`date +%Y%m%d -d "$((${year}-1))1231 +$i days"`;
      file_in=/albedo/work/projects/p_recompdaf/frbunsen/data/BGC/chlsurf_OC-CCI/original_daily/${year}/ESACCI-OC-L3S-CHLOR_A-MERGED-1D_DAILY_4km_GEO_PML_OCx-${datestr}-fv6.0.nc;
      file_out=/albedo/work/projects/p_recompdaf/frbunsen/data/BGC/chlsurf_OC-CCI/CORE2_daily/${year}/CCI-OC-${datestr}.nc;
      echo $datestr;
      cdo remapbil,/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/core2_griddes_nodes.nc ${file_in} ${file_out};
   done
done

# data is saved at: /albedo/work/projects/p_recompdaf/frbunsen/data/BGC/chlsurf_OC-CCI/CORE2_daily
# --> go there!

# merge daily files to yearly files with daily data
for((year=2010;year<=2020;year++))
do
# combine all in directory to 1-year-file:
cdo mergetime ${year}/*.nc CCI-OC-${year}-fullyear-nod2.nc

# rename dimensions:
ncrename -d ncells,nodes_2D  CCI-OC-${year}-fullyear-nod2.nc

# change missing value:
cdo setmissval,-999 CCI-OC-${year}-fullyear-nod2.nc CCI-OC-${year}-fullyear.nc
rm CCI-OC-${year}-fullyear-nod2.nc

done

# script to distribute observations onto FESOM-PEs is in: fesom/tools/dis_obs/chl_CCI/
# --> go there!

# submit script to distribute observations
for((year=2010;year<=2020;year++))
do
./distribute_obs ${year} > protocol_${year}.out
IsInLine=$( tail -1 protocol_${year}.out | grep -c 'END')
if [ ${IsInLine} -eq 1 ]; then
echo ${year} ' END'
fi
done
