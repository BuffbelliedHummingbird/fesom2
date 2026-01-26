! *****************************
! -----------------------------
! *** prestep_pdaf          ***
! -----------------------------
SUBROUTINE prestep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the prepoststep
! routine corresponding to the selected filter algorithm.
!
! This variant is adapted for multiple consecutive calls of the PDAF 
! assimilation routine. It is executed at the beginning of multiple
! consecutive assimilation steps.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022    - Frauke       - Adapted for FESOM2.0
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, COMM_filter, writepe, mype_world
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: step_null, filtertype, dim_lag, eff_dim_obs, loctype, &
       offset, proffiles_o, state_fcst, state_fcst_SSH_p, &
       monthly_state_f, monthly_state_a, monthly_state_m, &
       monthly_state_sf, monthly_state_sa, monthly_state_sm, &
       endday_of_month_in_year, startday_of_month_in_year, &
       depth_excl, depth_excl_no, this_is_pdaf_restart, mesh_fesom, nlmax, &
       dim_fields, dim_fields_glob, offset, offset_glob, nfields, id, &
       timemean, timemean_s, delt_obs_ocn, &
       days_since_DAstart, forget, &
       stdev_SSH_f_p, &
       factor_mass, factor_conc, DAoutput_path, &
       area_surf_glob, inv_area_surf_glob, &
       volo_full_glob, inv_volo_full_glob, &
       cellvol, &
       compute_monthly_aa, compute_monthly_ff, &
       compute_monthly_sa, compute_monthly_sf, &
       compute_monthly_mm, compute_monthly_sm, &
       resetforget, &
       count_lim_salt0_g , count_lim_salt0_p, &
       count_lim_absvel_g, count_lim_absvel_p, &
       count_lim_ssh_g   , count_lim_ssh_p, &
       count_lim_tempM2_g, count_lim_tempM2_p
  USE mod_atmos_ens_stochasticity, &
      ONLY: stable_rmse
  USE g_PARSUP, &
       ONLY: MPI_DOUBLE_PRECISION, MPI_SUM, MPIerr, MPI_STATUS_SIZE, &
       MPI_INTEGER, MPI_MAX, MPI_MIN, mydim_nod2d, MPI_COMM_FESOM, &
       myList_edge2D, myDim_edge2D, myList_nod2D
  USE o_ARRAYS, ONLY: hnode_new
  USE g_comm_auto, ONLY: gather_nod
  USE recom_config, &
       ONLY: tiny_chl, tiny, chl2N_max, chl2N_max_d, NCmax, &      
       NCmax_d, SiCmax, Redfield, SecondsPerDay
  USE mod_nc_out_routines, &
       ONLY: netCDF_out
  USE mod_nc_out_variables, &
       ONLY: sfields, nfields_3D, ids_3D, &
       w_dayensm, w_daymemb, w_monensm, w_monmemb, w_mm, w_sm, &
       mm, aa, ff, ii, sa, sf, si, sm, oo, dd, ee
  ! mean state forecast for observation exclusion criteria
  USE obs_TSprof_EN4_pdafomi, &
       ONLY: assim_o_en4_t, assim_o_en4_s, prof_exclude_diff, mean_temp_p
  USE obs_sst_pdafomi, &
       ONLY: assim_o_sst, sst_exclude_ice, sst_exclude_diff, &
             mean_ice_p, mean_sst_p
  USE obs_sss_smos_pdafomi, &
        ONLY: assim_o_sss, sss_exclude_ice, sss_exclude_diff, &
              mean_sss_p
  USE obs_sss_cci_pdafomi, &
        ONLY: assim_o_sss_cci, sss_cci_exclude_ice, sss_cci_exclude_diff, &
              mean_sss_cci_p
  USE obs_chl_cci_pdafomi, &
        ONLY: assim_o_chl_cci, chl_cci_exclude_ice, chl_cci_exclude_diff, &
              mean_chl_cci_p
  USE obs_o2_merged_pdafomi, &
        ONLY: assim_o_o2_merged, o2_merged_excl_absolute, o2_merged_excl_relative, &
              mean_O2_p
  USE obs_n_merged_pdafomi, &
        ONLY: assim_o_n_merged, n_merged_excl_absolute, n_merged_excl_relative, &
              mean_n_p
              
  USE g_clock, &
        ONLY: dayold, yearold, check_fleapyr, daynew, yearnew, &
              num_day_in_month, fleapyear, month, cyearnew, &
              day_in_month, timenew
  USE mod_assim_pdaf, &
        ONLY: debug_id_nod2
  USE mod_carbon_fluxes_diags
  USE netcdf
  USE g_events

  IMPLICIT NONE
  SAVE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step           ! Current time step, starting from 0 at beginning of year
                                        ! (When the routine is called before
                                        ! the analysis, -step is provided.)
  INTEGER, INTENT(in) :: dim_p          ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens        ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p      ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p      ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
                                        ! The array 'state_p' is not generally not initialized in the case of SEIK.
                                        ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag           ! PDAF status flag

! CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)

! *** Local variables ***
  INTEGER :: i,j,k,member,ed,s,n,f,ids,nz   ! Counters
  REAL :: invdim_ens                  ! Inverse ensemble size

  REAL :: diffm, aux                ! temporary arrays
  CHARACTER(len=1) :: typestr       ! Character indicating call type (intial, forecast, analysis)
  REAL :: min_eff_dim_obs, max_eff_dim_obs       ! Stats on effective observation dimensions
  REAL :: min_eff_dim_obs_g, max_eff_dim_obs_g   ! Stats on effective observation dimensions
  REAL :: sum_eff_dim_obs, avg_eff_dim_obs_g     ! Stats on effective observation dimensions
  LOGICAL :: now_to_write_monthly
  
  REAL :: tiny_N                 ! Min PhyN
  REAL :: tiny_N_d               ! Min DiaN
  REAL :: tiny_C                 ! Min PhyC
  REAL :: tiny_C_d               ! Min DiaC
  REAL :: tiny_Si                ! Min DiaSi
  REAL :: tiny_R                 ! Min ZoC
  
  REAL, ALLOCATABLE :: stdev_p(:)  ! ensemble standard deviation at grid proints
  REAL :: stdevglob_temp           ! global full ocean average of ensemble standard deviation for temperature field
                                                                    ! regional average of local ensemble standard deviation
  REAL :: stdev_surf_p(nfields,nlmax), stdev_surf_g(nfields,nlmax)  ! on surfaces
  REAL :: stdev_volo_p(nfields),       stdev_volo_g(nfields)        ! on full ocean volume
  
  INTEGER, parameter :: int0 = 0
  
  ! for carbon mass conservation diagnostics
  REAL :: weights
  REAL, allocatable  :: data3_g(:,:)                 ! Temporary array for global 3D-fields
  character(len=200) :: filename                     ! Full name of output file
  integer            :: fileid                       ! nc-file ID for output file
  character(len=100) :: varname

  ! variables for debugging:
  LOGICAL :: debug                  = .false.
  LOGICAL :: write_debug            = .false.
  LOGICAL :: simplify_debug         = .false.
  LOGICAL :: simplify_output        = .false.
  INTEGER :: fileID_debug
  CHARACTER(len=3) :: day_string
  INTEGER :: myDebug_id(1)
  LOGICAL :: debugging_monthlymean  = .false.
  
  simplify_debug  = .true. ! remove functionality for debugging purposes
  simplify_output = .true. ! remove output functionality for debugging purposes
  
  ! set debug output
  debug = .false.
  IF (.not. debug) THEN
     write_debug = .false.
  ELSE
     IF (mype_world>0) THEN
        write_debug = .false.
     ELSE
        write_debug = .true.
     ENDIF
  ENDIF

! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter==0) THEN
     IF ((step-step_null)==0) THEN
       IF (.not.(this_is_pdaf_restart)) THEN
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2,1x,a)') 'FESOM-PDAF', 'prestep_pdaf: WARNING: Nothing is done during initialization at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0),'(prepoststep_pdaf instead?)'
        WRITE (typestr,'(a1)') 'i'
       ELSE
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF','prestep_pdaf: WARNING: Nothing is done during initialization at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0),'(prepoststep_pdaf instead?)'
        WRITE (typestr,'(a1)') 'i'
       END IF
     ELSE IF ((step-step_null)>0) THEN
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF', 'prestep_pdaf: Nothing is done after physics assimilation at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0)
        WRITE (typestr,'(a1)') 'a'
     ELSE IF ((step-step_null)<0) THEN
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF', 'prestep_pdaf: Analyze forecast state ensemble before physics assimilation at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0)
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF ! IF (mype_filter==0)
  
  IF ((step-step_null)<0) THEN ! begin of pre-step
  
  ! variables allocated and saved during forecast; and deallocated after analysis
  IF (.not. ALLOCATED(stdev_SSH_f_p))    ALLOCATE(stdev_SSH_f_p(dim_fields(id%SSH)))
  IF (.not. ALLOCATED(state_fcst_SSH_p)) ALLOCATE(state_fcst_SSH_p(dim_fields(id%SSH),dim_ens))
  
! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! monthly event
  call monthly_event_assimstep(now_to_write_monthly)

! ****************************
! *** Corrections          ***
! ****************************

    IF ((step-step_null)<0) THEN
    ! Corrections
    
    IF (.not. simplify_debug) THEN ! simplify_debug-01
    ! *** store forecast state fields temporarily to compare with analysis afterwards ***    
      DO member = 1, dim_ens
        DO i = 1, dim_fields(id% SSH)
           state_fcst_SSH_p(i,member) = ens_p(i+offset(id% SSH),member)
        ENDDO
      ENDDO
    ENDIF ! simplify_debug-01
    
    END IF ! Corrections
   
! *******************************
! *** Compute ensemble mean   ***
! *******************************
  IF (.not. simplify_debug) THEN ! simplify_debug-04
  IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: compute ensemble mean'
  
  ! Local: 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i,member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)
  ENDIF ! simplify_debug-04


! *********************************************************************
! *** Store ensemble mean values for observation exclusion criteria ***
! *********************************************************************
! save values at forecast phase

  IF (.not. simplify_debug) THEN ! simplify_debug-05
  IF ((step-step_null)<0) THEN

     ! -- Sea-ice concentration --
     ! save mean_ice_p
     IF (    (sst_exclude_ice     .and. assim_o_sst     )&
        .OR. (sss_exclude_ice     .and. assim_o_sss     )&
        .OR. (sss_cci_exclude_ice .and. assim_o_sss_cci )&
        .OR. (chl_cci_exclude_ice .and. assim_o_chl_cci )) THEN 
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast SEA-ICE for observation exclusion'
        IF (ALLOCATED(mean_ice_p)) DEALLOCATE(mean_ice_p)
        ALLOCATE (mean_ice_p(dim_fields(id% a_ice)))
        mean_ice_p = state_p(offset(id% a_ice)+ 1 : &
                             offset(id% a_ice)+ dim_fields(id% a_ice))
     END IF

     ! -- SST --
     ! save mean_sst_p
     IF ((sst_exclude_diff > 0.0) .and. assim_o_sst) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast SST for observation exclusion'
        IF (ALLOCATED(mean_sst_p)) DEALLOCATE(mean_sst_p)
        ALLOCATE (mean_sst_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sst_p(i) = state_p(offset(id% temp) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! -- SSS (CASE SMOS) --
     ! save mean_sss_p
     IF ((sss_exclude_diff > 0.0) .and. assim_o_sss) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast SSS for observation exclusion'
        IF (ALLOCATED(mean_sss_p)) DEALLOCATE(mean_sss_p)
        ALLOCATE (mean_sss_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_p(i) = state_p(offset(id% salt) + (i-1) * (nlmax) + 1)
        END DO
     END IF

     ! -- SSS (CASE CCI) --
     ! save mean_sss_cci_p
     IF ((sss_cci_exclude_diff > 0.0) .and. assim_o_sss_cci) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast SSS for observation exclusion'
        IF (ALLOCATED(mean_sss_cci_p)) DEALLOCATE(mean_sss_cci_p)
        ALLOCATE (mean_sss_cci_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_cci_p(i) = state_p(offset(id% salt) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! -- Chlorophyll --
     ! save mean_chl_cci_p
     IF ((chl_cci_exclude_diff > 0.0) .and. assim_o_chl_cci) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast CHL for observation exclusion'
        IF (ALLOCATED(mean_chl_cci_p)) DEALLOCATE(mean_chl_cci_p)
        ALLOCATE (mean_chl_cci_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_chl_cci_p(i) = state_p(offset(id% PhyChl) + (i-1) * (nlmax) + 1) &
                            + state_p(offset(id% DiaChl) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! -- 3D temperature field --
     ! save mean_temp_p
     IF ((assim_o_en4_t .OR. assim_o_en4_s) &
         .AND. &
         (prof_exclude_diff > 0.0)) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean temperature (3D) for observation exclusion'
        ! Store mean temperature for profile assimilation
        IF (ALLOCATED(mean_temp_p)) DEALLOCATE(mean_temp_p)
        ALLOCATE (mean_temp_p(dim_fields(id%temp)))
        mean_temp_p = state_p(offset(id%temp)+1 : offset(id%temp)+dim_fields(id%temp))
     END IF
     
     ! -- oxygen --
     ! save mean_o2_p
     IF (((o2_merged_excl_absolute > 0.0) .or. (o2_merged_excl_relative > 0.0)) &
        .and. assim_o_o2_merged) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast oxygen for observation exclusion'
        IF (ALLOCATED(mean_O2_p)) DEALLOCATE(mean_O2_p)
        ALLOCATE (mean_O2_p(dim_fields(id%O2)))
        mean_O2_p = state_p(offset(id%O2)+1 : offset(id%O2)+dim_fields(id%O2))
     END IF
     
     ! -- nitrate --
     ! save mean_n_p
     IF (((n_merged_excl_absolute > 0.0) .or. (n_merged_excl_relative > 0.0)) &
        .and. assim_o_n_merged) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- prestep_pdaf: save ensemble mean forecast DIN for observation exclusion'
        IF (ALLOCATED(mean_n_p)) DEALLOCATE(mean_n_p)
        ALLOCATE (mean_n_p(dim_fields(id%DIN)))
        mean_n_p = state_p(offset(id%DIN)+1 : offset(id%DIN)+dim_fields(id%DIN))
     END IF

  END IF ! forecast phase
  END IF ! simplify_debug-05
  
! ************************************************
! *** Carbon sources minus sinks diagnostics   ***
! ************************************************
  
  IF (.not. simplify_debug) THEN ! simplify_debug-06
  ! factor to convert concentration to mass
  factor_mass = mesh_fesom%areasvol(:nlmax,:myDim_nod2D) * hnode_new(:nlmax,:myDim_nod2D) / SecondsPerDay
  factor_conc = 1.0 / SecondsPerDay

  IF ((step-step_null)<0) THEN
  ! forecast phase
  ! get fmass and fconc before analysis step
  
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- prestep_pdaf: compute carbon diagnostics at forecast'
  
    DO i = 1, myDim_nod2D
      DO k = 1, nlmax
      s = (i-1) * (nlmax) + k ! index in state vector
      ! DIC
      cffields(id_s_asml_dic)% fconc (k, i) = state_p(s + offset(id% DIC))
      ! Alk
      cffields(id_s_asml_alk)% fconc (k, i) = state_p(s + offset(id% Alk))
      ! Living carbon biomass
      cffields(id_s_asml_livingmatter)% fconc (k, i) = &
                               (state_p(s + offset(id% PhyC)) &
                              + state_p(s + offset(id% DiaC)) &
                              + state_p(s + offset(id% Zo1C)) &
                              + state_p(s + offset(id% Zo2C)) &
                              + state_p(s + offset(id% PhyCalc)))
      ! Dead organic carbon
      cffields(id_s_asml_deadmatter)% fconc (k, i) = &
                               (state_p(s + offset(id% DOC))     &
                              + state_p(s + offset(id% DetC))    &
                              + state_p(s + offset(id% DetCalc)) &
                              + state_p(s + offset(id% Det2C))   &
                              + state_p(s + offset(id% Det2Calc)))
      ENDDO ! k=1,nlmax
    ENDDO ! i=1,my_Dim_nod2D
    
    ! convert concentration to mass
    DO s=1, size(cffieldsasml)
       i = cffieldsasml(s)
       cffields(i)% fmass = cffields(i)% fconc * factor_mass
       cffields(i)% fconc = cffields(i)% fconc * factor_conc
    ENDDO

  ENDIF ! (forecast phase)
  ENDIF ! simplify_debug-06
  

! *****************************************************************
! *** Compute ensemble spread (STD) for different fields        ***
! *****************************************************************
  
  IF (.not. simplify_debug) THEN ! simplify_debug-07  
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- prestep_pdaf: compute ensemble standard deviation'
  
  ! Compute standard deviation of ensemble at grid points
  ! ---------------------------------------------------------------------------------------
  ! --- stdev_p (dim_p)    | standard deviation of ensemble at grid points for state vector
  ! ---------------------------------------------------------------------------------------
  ALLOCATE(stdev_p(dim_p))
  stdev_p(:) = 0.0
  DO member=1, dim_ens
     DO i=1, dim_p
        stdev_p(i) = stdev_p(i) & 
                + ((ens_p(i,member) - state_p(i)) * (ens_p(i,member) - state_p(i)))
     ENDDO ! j=1, dim_p
  ENDDO ! member=1, dim_ens
  stdev_p = SQRT(invdim_ens * stdev_p)
  
  ! if forecast: STD of SSH is saved and used for corrections at next analysis step
  IF ((step-step_null) < 0) then
     stdev_SSH_f_p = stdev_p( offset(id%SSH)+1 : offset(id%SSH)+dim_fields(id%SSH) )
  endif
  
  ! -----------------------------------------------------------------------------------------------------
  ! --- stdev_surf_g (nfields)    | layerwise surface mean of grid-point ensemble STD for each field area-weighted
  ! -----------------------------------------------------------------------------------------------------
  ! Compute pe-local surface mean of ensemble STD for each field
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- prestep_pdaf: compute ensemble standard deviation surface mean'
  stdev_surf_p = 0.0
  stdev_surf_g = 0.0
  
  DO f=1,nfields
     DO n=1,myDim_nod2D
           IF (sfields(f)%ndims == 2) THEN
           ! 3D fields
             DO nz=1,nlmax
                stdev_surf_p(f,nz) =  stdev_surf_p(f,nz) &
                                   +  mesh_fesom%areasvol(nz,n) * stdev_p( offset(f) + (n-1)*(nlmax) + nz )
             ENDDO ! nz,nlmax
           ELSE
           ! surface fields
             stdev_surf_p(f,1) =  stdev_surf_p(f,1) &
                               +  mesh_fesom%areasvol(1, n) * stdev_p( offset(f) + n)
           ENDIF
     ENDDO ! n, myDim_nod2D
  ENDDO ! f, nfields
  
  ! Reduce to global mean
  CALL MPI_Allreduce (stdev_surf_p, stdev_surf_g, nfields*nlmax, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_FESOM, MPIerr)
  DO nz=1,nlmax
     stdev_surf_g(:,nz) = stdev_surf_g(:,nz) * inv_area_surf_glob(nz)
  ENDDO
  
  ! -----------------------------------------------------------------------------------------------------
  ! --- stdev_volo_g (nfields)    | global mean of grid-point ensemble STD for each field volume-weighted
  ! -----------------------------------------------------------------------------------------------------
  ! Compute pe-local mean of ensemble STD for each field
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- prestep_pdaf: compute ensemble standard deviation global ocean mean'
  stdev_volo_p = 0.0
  stdev_volo_g = 0.0
  
  DO f=1,nfields
     DO n=1,myDim_nod2D
           IF (sfields(f)%ndims == 2) THEN
           ! 3D fields
             DO nz=1,nlmax
                stdev_volo_p(f) =  stdev_volo_p(f) &
                                +  cellvol(nz,n) * stdev_p( offset(f) + (n-1)*(nlmax) + nz )
             ENDDO ! nz,nlmax
           ENDIF
     ENDDO ! n,myDim_nod2D
  ENDDO ! f,nfields
  
  ! Reduce to global mean
  CALL MPI_Allreduce (stdev_volo_p, stdev_volo_g, nfields, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_FESOM, MPIerr)
  DO f=1,nfields
     IF (sfields(f)%ndims == 2) THEN
     ! 3D fields
       stdev_volo_g(f) = stdev_volo_g(f) * inv_volo_full_glob
     ELSE
     ! surface fields
       stdev_volo_g(f) = stdev_surf_g(f,1)
     ENDIF
  ENDDO
  
  ! Global 3D mean of temperature field used to tune ensemble inflation
  stdevglob_temp = stdev_volo_g(id%temp)

  ! Display RMS errors
  IF (mype_filter==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'FESOM-PDAF', 'prestep_pdaf: Ensemble standard deviation:'
     WRITE (*,'(a,7x,    a14,   a14,   a14,   a14,  a14, /a, 10x,70a)') &
          'FESOM-PDAF', 'CO2f','pCO2','temp','DIC','Alk', &
          'FESOM-PDAF', ('-',i=1,70)
     WRITE (*,'(a,10x,  5es14.4, 3x,a13,a1,/a, 10x,70a)')  &
          'FESOM-PDAF', stdev_surf_g(id% CO2f  ,1), &
                        stdev_surf_g(id% pCO2s ,1), &
                        stdev_surf_g(id% temp  ,1), &
                        stdev_surf_g(id% DIC   ,1), &
                        stdev_surf_g(id% Alk   ,1), &
                       'surface STDEV', typestr, 'FESOM-PDAF', ('-',i=1,70)
     WRITE (*,'(a,10x,  5es14.4, 3x,a13,a1,/a, 10x,70a)')  &
          'FESOM-PDAF', stdev_surf_g(id% CO2f  ,11), &
                        stdev_surf_g(id% pCO2s ,11), &
                        stdev_surf_g(id% temp  ,11), &
                        stdev_surf_g(id% DIC   ,11), &
                        stdev_surf_g(id% Alk   ,11), &
                       '90-100m STDEV', typestr, 'FESOM-PDAF', ('-',i=1,70)
     WRITE (*,'(a,10x,  5es14.4, 3x,a13,a1,/a, 10x,70a)')  &
          'FESOM-PDAF', stdev_volo_g(id% CO2f ), &
                        stdev_volo_g(id% pCO2s), &
                        stdev_volo_g(id% temp ), &
                        stdev_volo_g(id% DIC  ), &
                        stdev_volo_g(id% Alk  ), &
                       'vol oce STDEV', typestr, 'FESOM-PDAF', ('-',i=1,70)
  END IF
  END IF ! simplify_debug-07
  
! *******************************
! *** Reset forgetting factor ***
! *******************************
! Forgetting factor increases after the start of the assimilation (days_since_DAstart).
! In case of model restarts:
!  -  days_since_DAstart is set by slurm-job-script
!  -  current forgetting factor and target temperature ensemble standard deviation are read from atmos-perturbation file
  
  IF (.not. simplify_debug) THEN ! simplify_debug-08
  IF ((step-step_null)<0) THEN
  ! forecast phase   
     
     IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' days_since_DAstart ', days_since_DAstart
     
     IF (resetforget) THEN
     ! reset forgetting factor:
     ! set value for 1st half-month of assimilation
     IF     (days_since_DAstart==  1) THEN
       forget=0.95
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
      ! set value for 2nd half of 1st month of assimilation
     ELSEIF (days_since_DAstart== 16) THEN
       forget=0.96
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
      ! set value for 2nd month of assimilation 
     ELSEIF (days_since_DAstart== 32) THEN
       forget=0.97
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! set value for 3rd month of assimilation 
     ELSEIF (days_since_DAstart== 60) THEN
       forget=0.98
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! set value for 4th-17th months of assimilation
     ELSEIF (days_since_DAstart== 90) THEN
       forget=0.99
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! during month 17, save temperature ensemble standard deviation ("stable_rmse")
     ! month 17
     ELSEIF ((days_since_DAstart >= 485) .and. (days_since_DAstart <= 516)) THEN
       stable_rmse = stable_rmse + stdevglob_temp/31
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Saving ensemble spread to adapt forget at day ', days_since_DAstart
     
     ! after month 17, tune ensemble inflation based on the temperature ensemble standard deviation
     ! reset forgetting factor if ensemble standard deviation becomes larger/smaller than target value
     
     ! after month 17, higher ensemble standard deviation
     ELSEIF ((days_since_DAstart >= 516) .and. (stdevglob_temp > 1.1*stable_rmse)) THEN
       forget=1.00
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF '   , 'Resetting Forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' Target RMSE: ', stable_rmse, &
                                      ' New forget: ' , forget
                                      
     ! after month 17, lower ensemble standard deviation
     ELSEIF ((days_since_DAstart >= 516) .and. (stdevglob_temp < 0.9*stable_rmse)) THEN
       forget=0.99
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF '   , 'Resetting Forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' Target RMSE: ', stable_rmse, &
                                      ' New forget: ' , forget
     
     ! after month 17, ensemble standard deviation around target value
     ELSE
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF '   , 'Keeping Forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' Target RMSE: ', stable_rmse, &
                                      ' Forget:      ', forget
     
     ENDIF ! days_since_DAstart == ... (whether to reset forgetting factor according to scheme)
     ENDIF ! resetforget               (whether to use reset scheme for forgetting factor)
     
     ! day count during daily forecast phase:
     days_since_DAstart=days_since_DAstart+1
     
  ENDIF ! (forecast phase)
  ENDIF ! simplify_debug-08
   
! *****************************
! *** Compute monthly means ***
! *****************************
  
  IF (.not. simplify_debug) THEN ! simplify_debug-11
  debugging_monthlymean = .false.
  
     ! include state into monthly mean
     IF ((step-step_null) < 0) THEN
     ! *** forecasted fields ***
       IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,'(a, 8x, a,i7,a)') 'FESOM-PDAF', '--- prestep_pdaf: step ', step, 'adding forecast state to monthly mean.'
       IF (compute_monthly_ff) monthly_state_f  = monthly_state_f  + state_p
       IF (compute_monthly_sf) monthly_state_sf = monthly_state_sf + stdev_p
     END IF
     
     IF (now_to_write_monthly) THEN
     ! computing monthly mean at last day of month
     weights =  1.0/REAL(num_day_in_month(fleapyear,month))
     IF ((step-step_null) < 0) THEN
     ! *** forecasted state fields ***
       IF (compute_monthly_ff) monthly_state_f  = monthly_state_f  * weights
       IF (compute_monthly_sf) monthly_state_sf = monthly_state_sf * weights
     END IF
     ENDIF ! now_to_write_monthly
  ENDIF ! simplify_debug-11

! **************************
! *** Write output files ***
! **************************
! note: after monthly output is written, reset monthly fields to zero

  IF (.not. simplify_output) THEN ! simplify_debug (daily output)
  IF (.not. now_to_write_monthly) THEN
  IF ((step-step_null) < 0) THEN
        ! *** write forecast state fields ***
        ! during forecast phase, additionally, datetime and forgetting factor are written to output file
        ! ensemble mean
        IF (w_dayensm) CALL netCDF_out('ff',state_p,int0, now_to_write_monthly, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g, forget=forget)
        IF (w_dayensm) CALL netCDF_out('sf',stdev_p,int0, now_to_write_monthly)
        ! ensemble members
        IF (w_daymemb) THEN
          DO member = 1, dim_ens
            CALL netCDF_out('ff',ens_p(:,member), member, now_to_write_monthly)
          ENDDO
        ENDIF
  ENDIF
  ENDIF
  ENDIF ! simplify_debug (daily output)
  
  IF (.not. simplify_debug) THEN ! simplify_debug-13
  ! monthly output
  IF (now_to_write_monthly) THEN
  ! end of month: pass monthly output in addition to daily output
  IF ((step-step_null) < 0) THEN
        ! *** write forecast state fields ***
        ! ensemble mean, adding monthly mean of forecast states
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('ff',state_p, int0, now_to_write_monthly, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g, forget=forget, m_state_p=monthly_state_f )
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('sf',stdev_p, int0, now_to_write_monthly,                                                                      m_state_p=monthly_state_sf)
        ! ensemble members, adding snapshot
        IF (w_daymemb .or. w_monmemb) THEN
          DO member = 1, dim_ens
            CALL netCDF_out('ff',ens_p(:,member), member, now_to_write_monthly, m_state_p=ens_p(:,member))
          ENDDO
        ENDIF
  END IF
  END IF
  END IF ! simplify_debug-13

  IF (.not. simplify_debug) THEN ! simplify_debug-14
  ! at last day of month, reset monthly_state to zero (has been written)
  IF (now_to_write_monthly) THEN
     IF ((step-step_null) < 0) THEN
     ! *** forecasted state fields ***
       IF (compute_monthly_ff) monthly_state_f  = 0.0D0
       IF (compute_monthly_sf) monthly_state_sf = 0.0D0
     END IF
  ENDIF ! now_to_write_monthly
  ENDIF ! simplify_debug-14
  
! ********************
! *** finishing up ***
! ********************

  IF (allocated(stdev_p)) deallocate(stdev_p)
  
  ENDIF  ! end of pre-step

END SUBROUTINE prestep_pdaf
