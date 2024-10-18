#!/bin/bash

source ../../env/albedo/shell
source ./librunscript.sh # librunscript defines some helper functions
module load cdo

echo "Let's get started."

# Directories:
exp_name='dbfail'
start_dir=${PWD}
run_dir=/albedo/work/projects/p_recompdaf/frbunsen/modelruns/fesom2/${exp_name}/work/

echo ${run_dir}

# Settings:
NENS=2           # number of ensembles --> namelist.fesom.pdaf
fes_numproc=72   # ntasks divided by NENS

# Simulation start and end date
run_start_date="2010-01-01"
run_end_date="2020-08-01"

# Repeat from begin of year onwards:
repeat_date="2020-01-01"

# Restart frequency
# For runs without restart, leave this variable empty
rst_freq="8 month"
rst_freq_number=`echo ${rst_freq} | cut -c 1` # 1st character (number)
rst_freq_unit=`echo ${rst_freq} | cut -c 3`   # 3rd character (d,m,y)
# Number of restart legs to be run in one single slurm job
run_num_legs=1

#-----------------------------------------------------------------------

# Set $force_run_from_scratch to 'true' if you want to remove files
# from previous test runs from $run_dir.
# NOTE: If set to 'true', the run directory $run_dir is cleaned!
force_run_from_scratch=False

# This file is used to store information about restarts
ece_info_file="ece.info"

# NOTE: For reproduction runs, reset ece.info file manually to run_repro_date!

# ----------------------------------------------------------------------
# *** Read platform dependent configuration
# ----------------------------------------------------------------------

# File for standard output.
# NOTE: This will be modified for restart jobs!

# Resubmit this job for automatic restarts? [true/false]
# Also, add options for the resubmit command here.
resubmit_job=true
resubmit_opt=""

# ----------------------------------------------------------------------
# *** Go to run directory
#     Everything is done from here.
# ----------------------------------------------------------------------
cd ${run_dir}

# ----------------------------------------------------------------------
# *** Determine the time span of this run and whether it's a restart leg
# ----------------------------------------------------------------------

# Regularise the format of the start and end date of the simulation
run_start_date=$(date -uR -d "${run_start_date}")
run_end_date=$(date -uR -d "${run_end_date}")
repeat_date=$(date -uR -d "${repeat_date}")


# Update ece.info file:
echo "#"                                                     | tee -a ${run_dir}/${ece_info_file}
echo "# Decided to repeat last leg at `date '+%F %T'` from"  | tee -a ${run_dir}/${ece_info_file}
echo "leg_end_date=\"${repeat_date}\""                       | tee -a ${run_dir}/${ece_info_file}

# Check for restart information file and set the current leg start date

if ! [ -r ${ece_info_file} ]; then
    leg_is_restart=false
    leg_start_date=${run_start_date}
    leg_number=1
else
    leg_is_restart=true
    . ./${ece_info_file} # set leg start and end date from ece.info file
    leg_start_date=${leg_end_date}
    leg_number=$((leg_number+1))
fi

# Compute the end date of the current leg
if [ -n "${rst_freq}" ]; then # restart freq is not zero
    leg_end_date=$(date -uR -d "${leg_start_date} + ${rst_freq}")
else
    leg_end_date=${run_end_date}
fi

if [ $(date -d "${leg_end_date}" +%s) -gt $(date -d "${run_end_date}" +%s) ]; then
    leg_end_date=${run_end_date}
fi

# Some time variables needed later
leg_length_sec=$(( $(date -d "${leg_end_date}" +%s) - $(date -d "${leg_start_date}" +%s) ))
leg_start_sec=$(( $(date -d "${leg_start_date}" +%s) - $(date -d "${run_start_date}" +%s) ))
leg_end_sec=$(( $(date -d "${leg_end_date}" +%s) - $(date -d "${run_start_date}" +%s) ))
leg_start_date_yyyymmdd=$(date -u -d "${leg_start_date}" +%Y%m%d)
leg_start_date_yyyy=$(date -u -d "${leg_start_date}" +%Y)
leg_end_date_yyyy=$(date -u -d "${leg_end_date}" +%Y)

leg_start_dayofyear=$(date -d "${leg_start_date}" +%j)
step_null=$((32*(${leg_start_dayofyear}-1))) # 32 is number of time steps per day

# Compute days since beginning of first leg:
days_since_DAstart=$(( ($(date -d "${leg_start_date}" +%s) - $(date -d "${run_start_date}" +%s))/86400+1))

# Check whether there's actually time left to simulate - exit otherwise
if [ ${leg_length_sec} -le 0 ]; then
    echo "Leg start date equal to or after end of simulation."
    echo "Nothing left to do. Exiting."
    exit 0
fi

# Update namelists from jobscript settings
sed -i `grep -n dim_ens              ${start_dir}/namelist.fesom.pdaf | cut -d ':' -f 1`"c dim_ens=$NENS"                            ${start_dir}/namelist.fesom.pdaf

sed -i `grep -n restart_length_unit  ${start_dir}/namelist.config     | cut -d ':' -f 1`"c restart_length_unit='${rst_freq_unit}'"   ${start_dir}/namelist.config

sed -i `grep -n restart_length=      ${start_dir}/namelist.config     | cut -d ':' -f 1`"c restart_length=${rst_freq_number}"        ${start_dir}/namelist.config

sed -i `grep -n run_length=          ${start_dir}/namelist.config     | cut -d ':' -f 1`"c run_length=${rst_freq_number}"            ${start_dir}/namelist.config

sed -i `grep -n run_length_unit      ${start_dir}/namelist.config     | cut -d ':' -f 1`"c run_length_unit='${rst_freq_unit}'"       ${start_dir}/namelist.config

sed -i `grep -n this_is_pdaf_restart ${start_dir}/namelist.fesom.pdaf | cut -d ':' -f 1`"c this_is_pdaf_restart=.true."              ${start_dir}/namelist.fesom.pdaf

sed -i `grep -n step_null            ${start_dir}/namelist.fesom.pdaf | cut -d ':' -f 1`"c step_null=${step_null}"                   ${start_dir}/namelist.fesom.pdaf

sed -i `grep -n days_since_DAstart   ${start_dir}/namelist.fesom.pdaf | cut -d ':' -f 1`"c days_since_DAstart=${days_since_DAstart}" ${start_dir}/namelist.fesom.pdaf

# Path to JRA forcing
yearnew=$leg_start_date_yyyy
if (( yearnew >= 2020 && yearnew <= 2023 )); then
  echo "Setting to JRA 1.5 forcing"
  jrapath='/albedo/work/projects/p_pool_recom/forcing/JRA55-do-v1.5.0.1/'
elif (( yearnew >= 1958 && yearnew <= 2019 )); then
  echo "Setting to JRA 1.4 forcing"
  jrapath='/albedo/work/projects/p_pool_fesom1/forcing/JRA55-do-v1.4.0/'
else
  echo "No atmospheric forcing for this year!"
fi

sed -i `grep -n  nm_xwind_file  ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_xwind_file ='${jrapath}/uas.'"  ${start_dir}/namelist.forcing
sed -i `grep -n  nm_ywind_file  ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_ywind_file ='${jrapath}/vas.'"  ${start_dir}/namelist.forcing
sed -i `grep -n  nm_humi_file   ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_humi_file  ='${jrapath}/huss.'" ${start_dir}/namelist.forcing
sed -i `grep -n  nm_qsr_file    ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_qsr_file   ='${jrapath}/rsds.'" ${start_dir}/namelist.forcing
sed -i `grep -n  nm_qlw_file    ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_qlw_file   ='${jrapath}/rlds.'" ${start_dir}/namelist.forcing
sed -i `grep -n  nm_tair_file   ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_tair_file  ='${jrapath}/tas.'"  ${start_dir}/namelist.forcing
sed -i `grep -n  nm_prec_file   ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_prec_file  ='${jrapath}/prra.'" ${start_dir}/namelist.forcing
sed -i `grep -n  nm_snow_file   ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_snow_file  ='${jrapath}/prsn.'" ${start_dir}/namelist.forcing
sed -i `grep -n  nm_mslp_file   ${start_dir}/namelist.forcing | cut -d ':' -f 1`"c  nm_mslp_file  ='${jrapath}/psl.'"  ${start_dir}/namelist.forcing

# ------------------------------------------------------------------
# *** Prepare the run directory for a run from scratch
# ------------------------------------------------------------------


# Prepare directories
echo "Preparing start directory"

# Prepare start directory
cp ${start_dir}/../../bin/fesom.x ${start_dir}/
rm ${run_dir}/../start/*
cp -r ${start_dir}/* ${run_dir}/../start/

# Prepare work directories
for((i=1;i<=$NENS;i++))
do
    ENSstr=`printf %02d $i`
    ENS4str=`printf %04d $i`
    echo 'Preparing work directory '${ENS4str}
    
    cd ${run_dir}/${ENSstr}
    
    # Reset the clock file
    # rm fesom.clock
    clockfile="fesom.clock"
    yearnew=$leg_start_date_yyyy
    yearold=$(((${leg_start_date_yyyy}-1)))
    dayold=$(date -d "12/31/$yearold" +%j)
    daynew=1
    echo "83700   ${dayold}   ${yearold}" >> $clockfile
    echo "0       ${daynew}   ${yearnew}" >> $clockfile
    
    # Delete files from failed experiment
    # rm *${yearnew}*
    
    # if [ -d ${run_dir}/${ENSstr}/atmos ]; then
    #    rm ${run_dir}/${ENSstr}/atmos/*${yearnew}*
    # fi
    
    # Link atmosphere output for compatibility with old FESOM-PDAF code
    # ln -s atmos/*${yearold}* .

    # *** Link executable of model component
    rm fesom.x
    ln -s ${start_dir}/fesom.x .
    
    # *** Copy namelist files
    cp ${start_dir}/namelist.* .
    
    # set path in namelist for ensemble members
    sed -i `grep -n ResultPath    ${run_dir}/${ENSstr}/namelist.config     | cut -d ':' -f 1`"c ResultPath='${run_dir}/${ENSstr}/'"     ${run_dir}/${ENSstr}/namelist.config
    sed -i `grep -n DAoutput_path ${run_dir}/${ENSstr}/namelist.fesom.pdaf | cut -d ':' -f 1`"c DAoutput_path='${run_dir}/${ENSstr}/'"  ${run_dir}/${ENSstr}/namelist.fesom.pdaf

done # ((i=1;i<=$NENS;i++))
cd ${start_dir}
