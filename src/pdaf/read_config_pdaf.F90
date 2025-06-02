! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE read_config_pdaf()

! !DESCRIPTION:
! This routine read the namelist file with
! parameters controlling data assimilation with
! PDAF.
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3

! !USES:


! general assimilation settings
  USE mod_parallel_pdaf, &
       ONLY: mype_model, n_modeltasks, task_id
  USE mod_assim_pdaf, &
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, &
       offset, screen, filtertype, subtype, &
       delt_obs_ocn, &
       dim_bias, DA_couple_type, &
       incremental, type_forget, &
       forget, locweight, local_range, srange, &
       type_trans, type_sqrt, step_null, &
       eff_dim_obs, loc_radius, loctype, loc_ratio, &
       path_init, file_init, file_inistate, read_inistate, &
       twin_experiment, dim_obs_max, use_global_obs, DAoutput_path, &
       ASIM_START_USE_CLIM_STATE, this_is_pdaf_restart, &
       path_atm_cov, days_since_DAstart, assimilateBGC, assimilatePHY, &
       start_from_ENS_spinup, resetforget, &
       ! initial ensemble perturbation
       varscale, perturb_ssh, perturb_u, &
       perturb_v, perturb_temp, perturb_salt, &
       perturb_DIC, perturb_Alk, perturb_DIN, perturb_O2, &
       ! Temp-Salt-Profiles:
       path_obs_rawprof, file_rawprof_prefix, file_rawprof_suffix, &
       proffiles_o, start_year_o, end_year_o
       ! output
  USE mod_nc_out_variables, &
       ONLY: setoutput
       ! postprocessing
  USE mod_postprocess, &
       ONLY: isPP, yearPP, pathsim
       ! clock
  USE g_clock, &
       ONLY: yearold, yearnew
  
  
  ! physics observations
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

  ! biogeochem observations    
  USE obs_chl_cci_pdafomi, &
       ONLY: assim_o_chl_cci, rms_obs_chl_cci, path_obs_chl_cci, file_chl_cci_prefix, file_chl_cci_suffix, &
       chl_cci_exclude_ice, chl_cci_exclude_diff, bias_obs_chl_cci, chl_cci_fixed_rmse, chl_logarithmic, &
       path_bias, file_bias_prefix
  USE obs_DIC_glodap_pdafomi, &
       ONLY: assim_o_DIC_glodap, path_obs_DIC_glodap, &
       rms_obs_DIC_glodap, DIC_glodap_exclude_diff
  USE obs_Alk_glodap_pdafomi, &
       ONLY: assim_o_Alk_glodap, path_obs_Alk_glodap, &
       rms_obs_Alk_glodap, Alk_glodap_exclude_diff
  USE obs_pco2_SOCAT_pdafomi, &
       ONLY: assim_o_pCO2_SOCAT, path_obs_pCO2_SOCAT, &
       rms_obs_pCO2_SOCAT, pCO2_SOCAT_exclude_diff
  USE obs_o2_comf_pdafomi, &
       ONLY: assim_o_o2_comf, path_obs_o2_comf, &
       rms_obs_o2_comf, o2_comf_exclude_diff
  USE obs_n_comf_pdafomi, &
       ONLY: assim_o_n_comf, path_obs_n_comf, &
       rms_obs_n_comf, n_comf_exclude_diff
  USE obs_o2_argo_pdafomi, &
       ONLY: assim_o_o2_argo, path_obs_o2_argo, &
       rms_obs_o2_argo, o2_argo_exclude_diff
  USE obs_N_argo_pdafomi, &
       ONLY: assim_o_N_argo, path_obs_N_argo, &
       rms_obs_N_argo, N_argo_exclude_diff
  USE obs_o2_merged_pdafomi, &
       ONLY: assim_o_o2_merged, path_obs_o2_merged, &
       rms_obs_o2_merged, o2_merged_exclude_diff, &
       o2_merged_excl_absolute, o2_merged_excl_relative
  USE obs_N_merged_pdafomi, &
       ONLY: assim_o_N_merged, path_obs_N_merged, &
       rms_obs_N_merged, N_merged_exclude_diff, &
       n_merged_excl_absolute, n_merged_excl_relative
  
  ! ensemble initialization    
  USE mod_atmos_ens_stochasticity, &
       ONLY: disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       atmos_stochasticity_ON, write_atmos_st, &
       varscale_wind, varscale_tair, &
       varscale_humi, varscale_qlw
  USE mod_perturbation_pdaf, &
       ONLY: perturb_scale, &
       perturb_params_bio, perturb_params_phy
  


  IMPLICIT NONE


! Local variables
  CHARACTER(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file
  CHARACTER(len=32)  :: handle             ! Handle for command line parser
  LOGICAL :: printconfig = .TRUE.          ! Print information on all configuration parameters
  CHARACTER(len=4) :: year_string          ! Current year as string
  
       
  NAMELIST /pdaf/ filtertype, subtype, screen, &
       incremental, type_forget, forget, resetforget, dim_bias, &
       local_range, locweight, srange, DA_couple_type, &
       n_modeltasks, use_global_obs, &
       path_init, file_init, step_null, printconfig, &
       file_inistate, read_inistate, &
       type_trans, type_sqrt, dim_lag, &
       loctype, loc_ratio, delt_obs_ocn, &     
       dim_obs_max, &
       twin_experiment, &
       DAoutput_path, &
       ASIM_START_USE_CLIM_STATE, this_is_pdaf_restart, &
       days_since_DAstart, start_from_ENS_spinup, &
       assimilatePHY, &
       ! initial ensemble perturbation:
       varscale, perturb_ssh, perturb_u, &
       perturb_v, perturb_temp, perturb_salt, &
       perturb_DIC, perturb_Alk, perturb_DIN, perturb_O2, &
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
       ! parameter perturbation:
       perturb_scale, perturb_params_bio, perturb_params_phy, &
       ! BioGeoChemistry:
       assimilateBGC, &
       ! Chl-a CCI:
       assim_o_chl_cci, rms_obs_chl_cci, path_obs_chl_cci, &
       file_chl_cci_prefix, file_chl_cci_suffix, &
       chl_cci_exclude_ice, chl_cci_exclude_diff, &
       bias_obs_chl_cci, chl_cci_fixed_rmse, chl_logarithmic, &
       path_bias, file_bias_prefix, &
       ! DIC GLODAP:
       assim_o_DIC_glodap, path_obs_DIC_glodap, &
       rms_obs_DIC_glodap, DIC_glodap_exclude_diff, &
       ! Alk GLODAP:
       assim_o_Alk_glodap, path_obs_Alk_glodap, &
       rms_obs_Alk_glodap, Alk_glodap_exclude_diff, &
       ! pCO2 SOCAT:
       assim_o_pCO2_SOCAT, path_obs_pCO2_SOCAT, &
       rms_obs_pCO2_SOCAT, pCO2_SOCAT_exclude_diff, &
       ! O2 COMFORT:
       assim_o_o2_comf, path_obs_o2_comf, &
       rms_obs_o2_comf, o2_comf_exclude_diff, &
       ! DIN COMFORT:
       assim_o_n_comf, path_obs_n_comf, &
       rms_obs_n_comf, n_comf_exclude_diff, &
       ! O2 ARGO:
       assim_o_o2_argo, path_obs_o2_argo, &
       rms_obs_o2_argo, o2_argo_exclude_diff, &
       ! DIN ARGO:
       assim_o_N_argo, path_obs_N_argo, &
       rms_obs_N_argo, N_argo_exclude_diff, &
       ! O2 MERGED:
       assim_o_o2_merged, path_obs_o2_merged, &
       rms_obs_o2_merged, o2_merged_exclude_diff, &
       o2_merged_excl_absolute, o2_merged_excl_relative, &
       ! DIN MERGED:
       assim_o_n_merged, path_obs_n_merged, &
       rms_obs_n_merged, n_merged_exclude_diff, &
       n_merged_excl_absolute, n_merged_excl_relative

  NAMELIST /atmos_stoch/ &
       path_atm_cov, &
       disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       varscale_wind, varscale_tair, varscale_humi, varscale_qlw, &
       write_atmos_st
       
  NAMELIST /pp/ &
       isPP, yearPP, pathsim
       
  NAMELIST /pdafoutput/ &
       setoutput
              
! ****************************************************
! ***   Initialize PDAF parameters from namelist   ***
! ****************************************************

! *** Read namelist file ***
#ifdef DEBUG
  WRITE(*,*) 'Read PDAF namelist file: ',nmlfile
#endif

  OPEN (20,file=nmlfile)
  READ (20,NML =pdaf)
  CLOSE(20)

  OPEN (30,file=nmlfile)
  READ (30,NML =atmos_stoch)
  CLOSE(30)
  
  OPEN (40,file=nmlfile)
  READ (40,NML =pp)
  CLOSE(40)
  
  OPEN (50,file=nmlfile)
  READ (50,NML =pdafoutput)
  CLOSE(50)

! *** Add trailing slash to paths ***
  CALL add_slash(path_obs_sst)
  CALL add_slash(path_obs_sss)
  CALL add_slash(path_obs_ssh)
  
  CALL add_slash(path_obs_prof)
  CALL add_slash(path_init)
  
  CALL add_slash(path_obs_chl_cci)
  CALL add_slash(path_obs_DIC_glodap)
  CALL add_slash(path_obs_Alk_glodap)
  CALL add_slash(path_obs_pCO2_SOCAT)
  CALL add_slash(path_obs_O2_comf)
  CALL add_slash(path_obs_N_comf)
  
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

! Observation file prefixes:
WRITE(year_string,'(i4.4)') yearnew

file_sst_prefix = 'OSTIA_SST_'//TRIM(year_string)//'0101_'//TRIM(year_string)//'1231_daily_dist72_'
file_sss_prefix = 'SMOS_SSS_'//TRIM(year_string)//'_dist72_'
file_sss_cci_prefix = 'CCI_SSS_'//TRIM(year_string)//'_dist72_'
file_chl_cci_prefix = 'CCI_OC_'//TRIM(year_string)//'_dist72_'

! *** Print configuration variables ***
  showconf: IF (printconfig .AND. mype_model==0 .AND. task_id==1) THEN

     ! Overview of PDAF configuration
     WRITE (*,'(/a,1x,a)')           'FESOM-PDAF',   '-- Overview of PDAF configuration --'
     WRITE (*,'(a,3x,a)')            'FESOM-PDAF',   'PDAF [namelist: pdaf]:'
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'filtertype  ',         filtertype
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'subtype     ',         subtype
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'n_modeltasks',         n_modeltasks
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_ens     ',         dim_ens
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'delt_obs_ocn',         delt_obs_ocn
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'step_null   ',         step_null
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'days_since_DAstart',   days_since_DAstart
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'screen      ',         screen
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'incremental ',         incremental
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'type_forget ',         type_forget
     WRITE (*,'(a,5x,a20,1x,f10.4)') 'FESOM-PDAF',   'forget      ',         forget
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'resetforget ',         resetforget
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_bias    ',         dim_bias
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'type_trans  ',         type_trans
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'local_range ',         local_range
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'locweight   ',         locweight
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'loctype     ',         loctype
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'srange      ',         srange
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'loc_ratio   ',         loc_ratio
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'use_global_obs',       use_global_obs
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_lag     ',         dim_lag
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'DA_couple_type  ',     DA_couple_type
     
     ! Pre-processing of TS-profile data
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'proffiles_o  ',        proffiles_o
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'start_year_o ',        start_year_o
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'end_year_o   ',        end_year_o
     
     ! Physics observation-type settings
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'bias_obs_sst  ',       bias_obs_sst
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'bias_obs_prof ',       bias_obs_prof
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sst_exclude_ice',      sst_exclude_ice
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'sst_exclude_diff',     sst_exclude_diff
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'prof_exclude_diff',    prof_exclude_diff
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_sst     ',    TRIM(path_obs_sst)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sst_prefix  ',    TRIM(file_sst_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sst_suffix  ',    TRIM(file_sst_suffix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_prof    ',    TRIM(path_obs_prof)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_prof_prefix ',    TRIM(file_prof_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_prof_suffix ',    TRIM(file_prof_suffix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_rawprof    ', TRIM(path_obs_rawprof)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_rawprof_prefix ', TRIM(file_rawprof_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_rawprof_suffix ', TRIM(file_rawprof_suffix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'DAoutput_path ',       TRIM(DAoutput_path)
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_sss   ',       assim_o_sss
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_sss_cci',      assim_o_sss_cci
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_ssh   ',       assim_o_ssh
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_sst   ',       assim_o_sst
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_en4_t ',       assim_o_en4_t
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_en4_s ',       assim_o_en4_s
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_sst ',         rms_obs_sst
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_sss ',         rms_obs_sss
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_sss_cci ',     rms_obs_sss_cci
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_ssh ',         rms_obs_ssh
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'lradius_ssh ',         lradius_ssh
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_T   ',         rms_obs_T
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'rms_obs_S   ',         rms_obs_S
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'bias_obs_ssh',         bias_obs_ssh
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'ssh_fixed_rmse',       ssh_fixed_rmse
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sst_fixed_rmse',       sst_fixed_rmse
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sss_fixed_rmse',       sss_fixed_rmse
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'sss_cci_fixed_rmse',   sss_cci_fixed_rmse
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_sss     ',    TRIM(path_obs_sss)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_prefix  ',    TRIM(file_sss_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_suffix  ',    TRIM(file_sss_suffix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_sss_cci ',    TRIM(path_obs_sss_cci)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_cci_prefix',  TRIM(file_sss_cci_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_sss_cci_suffix',  TRIM(file_sss_cci_suffix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_ssh     ',    TRIM(path_obs_ssh)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_ssh_prefix  ',    TRIM(file_ssh_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_ssh_suffix  ',    TRIM(file_ssh_suffix)
     
     ! BGC observation-type settings
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_chl_cci',      assim_o_chl_cci
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_chl_cci_prefix',  TRIM(file_chl_cci_prefix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_chl_cci_suffix',  TRIM(file_chl_cci_suffix)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_chl_cci ',    TRIM(path_obs_chl_cci)
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'chl_cci_exclude_ice',  chl_cci_exclude_ice
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'chl_cci_fixed_rmse',   chl_cci_fixed_rmse
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_chl_cci',      rms_obs_chl_cci
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'chl_cci_exclude_diff', chl_cci_exclude_diff
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'bias_obs_chl_cci',     bias_obs_chl_cci
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'chl_logarithmic',      chl_logarithmic
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_bias',            TRIM(path_bias)
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_bias_prefix',     TRIM(file_bias_prefix)
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_DIC_glodap',      assim_o_DIC_glodap
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_DIC_glodap',     path_obs_DIC_glodap
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_DIC_glodap',      rms_obs_DIC_glodap
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'DIC_glodap_exclude_diff', DIC_glodap_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_Alk_glodap',      assim_o_Alk_glodap
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_Alk_glodap',     path_obs_Alk_glodap
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_Alk_glodap',      rms_obs_Alk_glodap
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'Alk_glodap_exclude_diff', Alk_glodap_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_pCO2_SOCAT',      assim_o_pCO2_SOCAT
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_pCO2_SOCAT',     path_obs_pCO2_SOCAT
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_pCO2_SOCAT',      rms_obs_pCO2_SOCAT
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'pCO2_SOCAT_exclude_diff', pCO2_SOCAT_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_n_comf',       assim_o_n_comf
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_n_comf',      path_obs_n_comf
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_n_comf',       rms_obs_n_comf
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_comf_exclude_diff',  N_comf_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_o2_comf',      assim_o_o2_comf
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_o2_comf',     path_obs_o2_comf
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_o2_comf',      rms_obs_o2_comf
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_comf_exclude_diff', o2_comf_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_o2_argo',      assim_o_o2_argo
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_o2_argo',     path_obs_o2_argo
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_o2_argo',      rms_obs_o2_argo
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_argo_exclude_diff', o2_argo_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_N_argo',       assim_o_N_argo
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_N_argo',      path_obs_N_argo
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_N_argo',       rms_obs_N_argo
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_argo_exclude_diff',  N_argo_exclude_diff
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_o2_merged',       assim_o_o2_merged
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_o2_merged',      path_obs_o2_merged
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_o2_merged',       rms_obs_o2_merged
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_merged_exclude_diff',  o2_merged_exclude_diff
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_merged_excl_absolute', o2_merged_excl_absolute
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'o2_merged_excl_relative', o2_merged_excl_relative
     
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assim_o_N_merged',       assim_o_n_merged
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_obs_N_merged',      path_obs_n_merged
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'rms_obs_N_merged',       rms_obs_n_merged
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_merged_exclude_diff',  n_merged_exclude_diff
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_merged_excl_absolute', n_merged_excl_absolute
     WRITE (*,'(a,5x,a20,1x,f11.3)') 'FESOM-PDAF',   'N_merged_excl_relative', n_merged_excl_relative
          
     ! initial ensemble covariance
     WRITE (*,'(a,5x,a20,1x,1x,a)')     'FESOM-PDAF',   'path_init   ',         TRIM(path_init)
     WRITE (*,'(a,5x,a20,1x,1x,a)')     'FESOM-PDAF',   'file_init   ',         TRIM(file_init)
     WRITE (*,'(a,5x,a20,1x,1x,l)')     'FESOM-PDAF',   'this_is_pdaf_restart', this_is_pdaf_restart
     WRITE (*,'(a,5x,a20,1x,1x,l)')     'FESOM-PDAF',   'start_from_ENS_spinup',start_from_ENS_spinup
     WRITE (*,'(a,5x,a20,1x,1x,l)')     'FESOM-PDAF',   'ASIM_START_USE_CLIM_STATE', ASIM_START_USE_CLIM_STATE
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale    ', varscale
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_ssh',  perturb_ssh
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_u',    perturb_u
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_v',    perturb_v
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_temp', perturb_temp
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_salt', perturb_salt
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_DIC',  perturb_DIC
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_Alk',  perturb_Alk
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_DIN',  perturb_DIN
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_O2',   perturb_O2
     
     ! Twin experiment
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'twin_experiment', twin_experiment
     IF (filtertype==100 .or. twin_experiment) THEN
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'dim_obs_max ',    dim_obs_max
     END IF
     IF (read_inistate) THEN
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'file_inistate ',  TRIM(file_inistate)
     ENDIF

     ! Atmospheric perturbation
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'path_atm_cov  ',       TRIM(path_atm_cov)
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'atmos_stochasticity_ON',atmos_stochasticity_ON
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'write_atmos_st',       write_atmos_st
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_xwind',        disturb_xwind
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_ywind',        disturb_ywind
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_humi',         disturb_humi
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_qlw',          disturb_qlw
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_qsr',          disturb_qsr
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_tair',         disturb_tair
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_prec',         disturb_prec
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_snow',         disturb_snow
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'disturb_mslp',         disturb_mslp
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_wind ',       varscale_wind
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_tair ',       varscale_tair
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_humi ',       varscale_humi
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'varscale_qlw  ',       varscale_qlw
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_params_bio',   perturb_params_bio
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'perturb_params_phy',   perturb_params_phy
     WRITE (*,'(a,5x,a20,1x,es10.2)')'FESOM-PDAF',   'perturb_scale',        perturb_scale
     
     ! Updated state variables
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assimilatePHY',        assimilatePHY
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'assimilateBGC',        assimilateBGC
     
     ! Postprocessing
     WRITE (*,'(a,5x,a20,1x,l)')     'FESOM-PDAF',   'isPP',                 isPP
     WRITE (*,'(a,5x,a20,1x,i10)')   'FESOM-PDAF',   'yearPP',               yearPP
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'pathsim(1)',           TRIM(pathsim(1))
     WRITE (*,'(a,5x,a20,1x,a)')     'FESOM-PDAF',   'pathsim(2)',           TRIM(pathsim(2))

     WRITE (*,'(a,1x,a)') 'FESOM-PDAF','-- End of PDAF configuration overview --'
     

  END IF showconf

END SUBROUTINE read_config_pdaf



! ==============================================================================
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

! *** Local variables ***
  INTEGER :: strlength

! *** Add trailing slash ***
  strlength = LEN_TRIM(path)

  IF (path(strlength:strlength) /= '/') THEN
     path = TRIM(path) // '/'
  END IF

END SUBROUTINE add_slash
