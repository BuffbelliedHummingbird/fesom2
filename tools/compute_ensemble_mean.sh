files='ls 02/*.fesom.2016.nc'
for file in $files
do
   # get names of output variables
   # remove first 3 letters (02/) and last 14 letters (.fesom.2016.nc)
   var=$(echo $file | sed 's/^...//;s/..............$//')
   echo $var

   # calculate ensemble statistics:
   cdo ensmean ??/"$var".fesom.2016.nc ensemble_means/"$var".mean.2016.nc
done
