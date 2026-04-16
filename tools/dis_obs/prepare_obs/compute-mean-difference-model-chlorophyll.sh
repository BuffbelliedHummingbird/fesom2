### compute with output of FREE (after having postprocessed chlorophyll)

griddes="/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/core2_griddes_nod2.nc"

for((year=2010;year<=2020;year++))
do

### Extract the bias field
ncatted -a ancillary_variables,obsdiff,d,, extract-chl-surf-log-diff.${year}.nc
ncks -v obsdiff extract-chl-surf-log-diff.${year}.nc chlobsdiff.${year}.nc

### Interpolate to 0.5x grid:
time cdo -remapycon,global_1 -setgrid,$griddes chlobsdiff.${year}.nc chlobsdiff.${year}.regular.nc
# timer: real 17m54.222s

### Smooth:

# # Test with data for one day:
# cdo -selday,1,1 chlobsdiff.${year}.regular.nc chlobsdiff.${year}.regular.1.nc
# time cdo -smooth,radius=1000km chlobsdiff.${year}.regular.1.nc chlobsdiff.${year}.smooth.1.nc

# time cdo -smooth,radius=1000km chlobsdiff.${year}.regular.nc chlobsdiff.${year}.smooth.nc

# --> skip here and smooth later!

done

# calculating daily climatology 2010-2020

# check if fix for time is required:
year1=2010
cdo -P 128 -cat 'chlobsdiff.????.regular.nc' chlobsdiff.cat.nc
cdo settaxis,${year1}-01-01,00:00:00,1day chlobsdiff.cat.nc chlobsdiff.catfix.nc
cdo -P 128 -ydaymean chlobsdiff.catfix.nc chlobsdiff.clim.nccdo settaxis,${year}-01-01,00:00:00,1day chlobsdiff.cat.nc chlobsdiff.catfix.nc

# if no fix required, simply:
# time cdo -ydaymean -cat 'chlobsdiff.????.smooth.nc' chlobsdiff.clim.nc

# extend climatology circularly - check which of following commands works ...

# cdo cat -selday,-15/-1 chlobsdiff.clim.nc chlobsdiff.clim.nc -selday,1/15 chlobsdiff.clim.nc chlobsdiff.clim-extend.nc
# alternatively extend climatology:
# cdo cat -shifttime,-15days -selday,-15/-1 chlobsdiff.clim.nc chlobsdiff.clim.nc -shifttime,15days -selday,1/15 chlobsdiff.clim.nc chlobsdiff.clim-extend.nc
# cdo cat chlobsdiff.clim.nc -shifttime,15days -selday,4018,4032 chlobsdiff.clim.nc chlobsdiff.clim-extend.nc

# seems to work:
cdo -shifttime,15days chlobsdiff.clim.nc chlobsdiff.clim.plus15.nc
cdo -shifttime,-15days chlobsdiff.clim.nc chlobsdiff.clim.min15.nc
cdo -seltimestep,1/15 chlobsdiff.clim.min15.nc chlobsdiff.clim-m15end.nc
cdo -seltimestep,352/366 chlobsdiff.clim.plus15.nc chlobsdiff.clim-p15start.nc
cdo mergetime chlobsdiff.clim-m15end.nc chlobsdiff.clim.nc chlobsdiff.clim-p15start.nc chlobsdiff.clim-extend.nc

# smooth climatology:
cdo runmean,31 chlobsdiff.clim-extend.nc chlobsdiff.clim-extend.runmean.nc

# back to irregular grid:
path='/albedo/work/projects/p_recompdaf/frbunsen/modelruns/fesom2/update_BGC_tests/FREE/work/01'

cdo -P 128 runmean,31 chlobsdiff.clim-extend.nc chlobsdiff.clim-extend.runmean.nc
cdo -P 128 -smooth,radius=1000km chlobsdiff.clim-extend.runmean.nc chlobsdiff.clim-extend.runmean.smooth.nc

cdo -P 128 remapbil,${griddes} chlobsdiff.clim-extend.runmean.smooth.nc chlobsdiff.runmean.smooth.core2.nc
cdo -P 128 remapbil,${griddes} chlobsdiff.clim-extend.runmean.nc chlobsdiff.runmean.core2.nc
