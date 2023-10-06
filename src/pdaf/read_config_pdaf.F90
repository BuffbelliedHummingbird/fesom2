! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE read_config_pdaf()

! !DESCRIPTION:
! This routine read the namelist file with
! parameters controlling data assimilation with
! PDAF.
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3

! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_model, n_modeltasks, task_id
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, &
       offset, screen, filtertype, subtype, dim_ens, &
       delt_obs_ocn, &
       dim_bias, DA_couple_type,write_monthly_mean, &
       incremental, type_forget, peak_obs_error, &
       forget, locweight, local_range, srange, &
       type_trans, type_sqrt, step_null, &
       eff_dim_obs, loc_radius, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, varscale, &
       twin_experiment, dim_obs_max, use_global_obs, DAoutput_path, &
       ASIM_START_USE_CLIM_STATE, this_is_pdaf_restart, &
       path_atm_cov, days_since_DAstart, &
       ! Temp-Salt-Profiles:
       path_obs_rawprof, file_rawprof_prefix, file_rawprof_suffix, &
       proffiles_o, start_year_o, end_year_o
  USE mod_nc_out_variables, &
       ONLY: write_ens
  USE obs_sst_pdafomi, &
       ONLY: assim_o_sst, rms_obs_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       sst_exclude_ice, sst_exclude_diff, bias_obs_sst, sst_fixed_rmse
  USE obs_sss_smos_pdafomi, &
       ONLY: assim_o_sss, rms_obs_sss, path_obs_sss, file_sss_prefix, file_sss_suffix, &
       sss_exclude_ice, sss_exclude_diff, bias_obs_sss, sss_fixed_rmse
  USE obs_sss_cci_pdafomi, &
       ONLY: assim_o_sss_cci, rms_obs_sss_cci, path_obs_sss_cci, file_sss_cci_prefix, file_sss_cci_suffix, &
       sss_cci_exclude_ice, sss_cci_exclude_diff, bias_obs_sss_cci, sss_cci_fixed_rmse
  USE obs_ssh_cmems_pdafomi, &
       ONLY: assim_o_ssh, rms_obs_ssh, path_obs_ssh, file_ssh_prefix, file_ssh_suffix, &
       ssh_exclude_ice, ssh_exclude_diff, bias_obs_ssh, ssh_fixed_rmse, lradius_ssh
  USE obs_TSprof_EN4_pdafomi, &
       ONLY: ASSIM_O_en4_t, ASSIM_O_en4_s, & 
       path_obs_prof, file_prof_prefix, file_prof_suffix, &
       rms_obs_S, rms_obs_T, &
       file_syntobs_prof, prof_exclude_diff, bias_obs_prof
  USE mod_atmos_ens_stochasticity, &
       ONLY: disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       atmos_stochasticity_ON, write_atmos_st, &
       varscale_wind, varscale_tair, &
       varscale_humi, varscale_qlw
  USE g_clock, &
       ONLY: yearold, yearnew


  IMPLICIT NONE
!EOP

! Local variables
  CHARACTER(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file
  CHARACTER(len=32)  :: handle             ! Handle for command line parser
  LOGICAL :: printconfig = .TRUE.          ! Print information on all configuration parameters
  CHARACTER(len=4) :: year_string          ! Current year as string
  
       
  NAMELIST /pdaf/ filtertype, subtype, dim_ens, screen, &
       incremental, type_forget, forget, dim_bias, &
       local_range, locweight, srange, DA_couple_type, &
       n_modeltasks, peak_obs_error, use_global_obs, &
       path_init, file_init, step_null, printconfig, &
       file_inistate, read_inistate, write_ens, varscale, &
       type_trans, type_sqrt, dim_lag, &
       loctype, loc_ratio, delt_obs_ocn, &     
       dim_obs_max, &
       twin_experiment, &
       write_monthly_mean, &
       DAoutput_path, &
       ASIM_START_USE_CLIM_STATE, this_is_pdaf_restart, &
       days_since_DAstart, &
       ! Salt SMOS:
       ASSIM_o_sss, path_obs_sss, file_sss_prefix, file_sss_suffix, &
       rms_obs_sss, sss_fixed_rmse, &
       sss_exclude_ice, sss_exclude_diff, &
       ! Salt CCI:
       ASSIM_o_sss_cci, path_obs_sss_cci, file_sss_cci_prefix, file_sss_cci_suffix, &
       rms_obs_sss_cci, sss_cci_fixed_rmse, &
       sss_cci_exclude_ice, sss_cci_exclude_diff, &
       ! SSH:
       ASSIM_o_ssh, path_obs_ssh, file_ssh_prefix, file_ssh_suffix, &
       rms_obs_ssh, ssh_fixed_rmse, bias_obs_ssh, lradius_ssh, &
       ! SST:
       ASSIM_o_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       rms_obs_sst, sst_fixed_rmse, bias_obs_sst, &
       sst_exclude_ice, sst_exclude_diff, &
       ! Profiles:
       ASSIM_o_en4_t, ASSIM_o_en4_S, &
       path_obs_prof, file_prof_prefix, file_prof_suffix, &
       rms_obs_S, rms_obs_T, &
       path_obs_rawprof, file_rawprof_prefix, proffiles_o, &
       start_year_o, end_year_o, &
       ! Atmosphere cov:
       path_atm_cov

  NAMELIST /atmos_stoch/ disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       varscale_wind, varscale_tair, varscale_humi, varscale_qlw, &
       write_atmos_st
       
!~   NAMELIST /pdaf_output/
       
! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  WRITE(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  OPEN (20,file=nmlfile)
  READ (20,NML=pdaf)
  CLOSE (20)

  OPEN(30,file=nmlfile)
  READ(30,NML=atmos_stoch)
  CLOSE(30)
  
!~   OPEN(40,file=nmlfile)
!~   READ(40,NML=pdaf_output)
!~   CLOSE(40)

! *** Add trailing slash to paths ***
  CALL add_slash(path_obs_sst)
  CALL add_slash(path_obs_sss)
  CALL add_slash(path_obs_ssh)
  
  CALL add_slash(path_obs_prof)
  CALL add_slash(path_init)
  
! *** Is atmospheric stochasticity used at all?
IF (disturb_humi   .OR. &
    disturb_mslp   .OR. &
    disturb_xwind  .OR. &
    disturb_ywind  .OR. &
    disturb_qlw    .OR. &
    disturb_qsr    .OR. &
    disturb_tair   .OR. &
    disturb_prec   .OR. &
    disturb_snow   ) THEN
    
    atmos_stochasticity_ON = .TRUE.
ELSE
    atmos_stochasticity_ON = .FALSE.
ENDIF

IF (write_ens) WRITE (*,*) 'FESOM-PDAF ', '*** WARNING *** ', 'Ensemble output for day-average (m) NOT AVAILABLE. Hint: Activate model ouput in fvom_main.F90 instead.'

! Observation file prefixes:
WRITE(year_string,'(i4.4)') yearnew

file_sst_prefix = 'OSTIA_SST_'//TRIM(year_string)//'0101_'//TRIM(year_string)//'1231_daily_dist72_'
file_sss_prefix = 'SMOS_SSS_'//TRIM(year_string)//'_dist72_'
file_sss_cci_prefix = 'CCI_SSS_'//TRIM(year_string)//'_dist72_'

! *** Print configuration variables ***
  showconf: IF (printconfig .AND. mype_model==0 .AND. task_id==1) THEN

     WRITE (*,'(/a,1x,a)') 'FESOM-PDAF','-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)') 'FESOM-PDAF','PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','filtertype  ', filtertype
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','subtype     ', subtype
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','n_modeltasks', n_modeltasks
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_ens     ', dim_ens
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','delt_obs_ocn', delt_obs_ocn
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','step_null   ', step_null
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF', 'days_since_DAstart', days_since_DAstart
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','screen      ', screen
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','incremental ', incremental
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','type_forget ', type_forget
     WRITE (*,'(a,5x,a,f10.4)') 'FESOM-PDAF','forget      ', forget
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','varscale   ', varscale
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_bias    ', dim_bias
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','type_trans  ', type_trans
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','local_range ', local_range
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','locweight   ', locweight
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','loctype     ', loctype
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','srange      ', srange
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','lradius_ssh ', lradius_ssh
     WRITE (*,'(a,5x,a,es10.2)') 'FESOM-PDAF','loc_ratio   ', loc_ratio
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','proffiles_o  ', proffiles_o
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','start_year_o ', start_year_o
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','end_year_o   ', end_year_o
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','peak_obs_error', peak_obs_error
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','bias_obs_sst  ', bias_obs_sst
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','bias_obs_prof ', bias_obs_prof
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sst_exclude_ice', sst_exclude_ice
     WRITE (*,'(a,5x,a,f11.3)') 'FESOM-PDAF','sst_exclude_diff', sst_exclude_diff
     WRITE (*,'(a,5x,a,f11.3)') 'FESOM-PDAF','prof_exclude_diff', prof_exclude_diff
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','use_global_obs', use_global_obs
     WRITE (*,'(a,5x,a,i10)')     'FESOM-PDAF','use_global_obs', use_global_obs
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_lag     ', dim_lag
     WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','DA_couple_type  ', DA_couple_type
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_sst     ', TRIM(path_obs_sst)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sst_prefix  ', TRIM(file_sst_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sst_suffix  ', TRIM(file_sst_suffix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_prof    ', TRIM(path_obs_prof)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_prof_prefix ', TRIM(file_prof_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_prof_suffix ', TRIM(file_prof_suffix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_rawprof    ', TRIM(path_obs_rawprof)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_rawprof_prefix ', TRIM(file_rawprof_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_rawprof_suffix ', TRIM(file_rawprof_suffix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','DAoutput_path ', TRIM(DAoutput_path)
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_sss   ', assim_o_sss
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_sss_cci', assim_o_sss_cci
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_ssh   ', assim_o_ssh
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_sic   ', assim_o_sic
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_sit   ', assim_o_sit
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_siu   ', assim_o_siu
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_siv   ', assim_o_siv
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_sst   ', assim_o_sst
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_en4_t ', assim_o_en4_t
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','assim_o_en4_s ', assim_o_en4_s
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_sst ', rms_obs_sst
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_sss ', rms_obs_sss
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_sss_cci ', rms_obs_sss_cci
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_ssh ', rms_obs_ssh
!~      WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_sic ', rms_obs_sic
!~      WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_sit ', rms_obs_sit
!~      WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_siu ', rms_obs_siu
!~      WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_siv ', rms_obs_siv
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_T   ', rms_obs_T
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','rms_obs_S   ', rms_obs_S
     WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF','bias_obs_ssh', bias_obs_ssh
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sic_fixed_rmse', sic_fixed_rmse
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sit_fixed_rmse', sit_fixed_rmse
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','siu_fixed_rmse', siu_fixed_rmse
!~      WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','siv_fixed_rmse', siv_fixed_rmse
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','ssh_fixed_rmse', ssh_fixed_rmse
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sst_fixed_rmse', sst_fixed_rmse
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sss_fixed_rmse', sss_fixed_rmse
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','sss_cci_fixed_rmse', sss_cci_fixed_rmse
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_sss     ', TRIM(path_obs_sss)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sss_prefix  ', TRIM(file_sss_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sss_suffix  ', TRIM(file_sss_suffix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_sss_cci ', TRIM(path_obs_sss_cci)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sss_cci_prefix  ', TRIM(file_sss_cci_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_sss_cci_suffix  ', TRIM(file_sss_cci_suffix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_obs_ssh     ', TRIM(path_obs_ssh)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_ssh_prefix  ', TRIM(file_ssh_prefix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_ssh_suffix  ', TRIM(file_ssh_suffix)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_init   ', TRIM(path_init)
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_init   ', TRIM(file_init)
     IF (filtertype==100 .or. twin_experiment) THEN
        WRITE (*,'(a,5x,a,i10)')   'FESOM-PDAF','dim_obs_max ', dim_obs_max
     END IF
     IF (read_inistate) THEN
        WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','file_inistate ', TRIM(file_inistate)
     ENDIF
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','write_ens            ', write_ens
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','write_monthly_mean', write_monthly_mean
     
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','twin_experiment', twin_experiment
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','ASIM_START_USE_CLIM_STATE', ASIM_START_USE_CLIM_STATE 
     WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF','this_is_pdaf_restart', this_is_pdaf_restart
     WRITE (*,'(a,5x,a,a)')     'FESOM-PDAF','path_atm_cov  ', TRIM(path_atm_cov)

WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'atmos_stochasticity_ON', atmos_stochasticity_ON
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'write_atmos_st', write_atmos_st
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_xwind', disturb_xwind
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_ywind', disturb_ywind
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_humi', disturb_humi
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_qlw', disturb_qlw
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_qsr', disturb_qsr
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_tair', disturb_tair
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_prec', disturb_prec
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_snow', disturb_snow
WRITE (*,'(a,5x,a,l)')     'FESOM-PDAF', 'disturb_mslp', disturb_mslp
WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF', 'varscale_wind ', varscale_wind
WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF', 'varscale_tair ', varscale_tair
WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF', 'varscale_humi ', varscale_humi
WRITE (*,'(a,5x,a,es10.2)')'FESOM-PDAF', 'varscale_qlw  ', varscale_qlw


     WRITE (*,'(a,1x,a)') 'FESOM-PDAF','-- End of PDAF configuration overview --'

  END IF showconf

END SUBROUTINE read_config_pdaf
! ==============================================================================
!BOP
!
! !ROUTINE: add_slash --- Add trailing slash to path string
!
! !INTERFACE:
SUBROUTINE add_slash(path)

! !DESCRIPTION:
! This routine ensures that a string defining a path
! has a trailing slash.
!
! !USES:
  IMPLICIT NONE

! !ARGUMENTS:
  CHARACTER(len=100) :: path  ! String holding the path
!EOP

! *** Local variables ***
  INTEGER :: strlength

! *** Add trailing slash ***
  strlength = LEN_TRIM(path)

  IF (path(strlength:strlength) /= '/') THEN
     path = TRIM(path) // '/'
  END IF

END SUBROUTINE add_slash
