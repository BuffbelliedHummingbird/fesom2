! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
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
       resetforget
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
  INTEGER, INTENT(in) :: step        ! Current time step, starting from 0 at beginning of year
                                     ! (When the routine is called before
                                     ! the analysis, -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

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
  INTEGER, allocatable :: count_lim_salt0_g(:) , count_lim_salt0_p(:)       ! Count how many excessively large updates are limited to treshold
  INTEGER, allocatable :: count_lim_absvel_g(:), count_lim_absvel_p(:)      ! Count how many excessively large updates are limited to treshold
  INTEGER, allocatable :: count_lim_ssh_g(:)   , count_lim_ssh_p(:)         ! Count how many excessively large updates are limited to treshold
  INTEGER, allocatable :: count_lim_tempM2_g(:), count_lim_tempM2_p(:)      ! Count how many excessively large updates are limited to treshold
  
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
  LOGICAL :: debug
  LOGICAL :: write_debug
  INTEGER :: fileID_debug
  CHARACTER(len=3) :: day_string
  INTEGER :: myDebug_id(1)
  LOGICAL :: debugging_monthlymean = .false.

! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter==0) THEN
     IF ((step-step_null)==0) THEN
       IF (.not.(this_is_pdaf_restart)) THEN
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF', 'Analyze initial state ensemble at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0)
        WRITE (typestr,'(a1)') 'i'
       ELSE
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF','This is a PDAF restart. No initial fields at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0)
        WRITE (typestr,'(a1)') 'i'
       END IF
     ELSE IF ((step-step_null)>0) THEN
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF', 'Analyze assimilated state ensemble at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0)
        WRITE (typestr,'(a1)') 'a'
     ELSE IF ((step-step_null)<0) THEN
        WRITE (*,'(a, 8x,a,1x,i7,2x,i4,a1,i2.2,a1,i2.2,2x,i2.2,a1,i2.2)') 'FESOM-PDAF', 'Analyze forecast state ensemble at step', &
        step, yearnew,'-',month,'-',day_in_month,FLOOR(timenew/3600.0),':',INT(MOD(timenew,3600.0)/60.0)
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF ! IF (mype_filter==0)
  
  ! variables allocated and saved during forecast; and deallocated after analysis
  IF (.not. ALLOCATED(stdev_SSH_f_p))    ALLOCATE(stdev_SSH_f_p(dim_fields(id%SSH)))
  IF (.not. ALLOCATED(state_fcst_SSH_p)) ALLOCATE(state_fcst_SSH_p(dim_fields(id%SSH),dim_ens))
  
  IF ((step-step_null)==0) THEN
  ! allocate monthly states at initial time; never de-allocated; monthly reset to zero
  DO s=1, nfields
     !                                                  -- output? -----------        -- ensemble mean? ----------            -- monthly? ----------------------
     compute_monthly_ff = compute_monthly_ff .or. ( sfields(s)%output(ff,oo) .and. (.not. sfields(s)%output(ff,ee)) .and. (.not. sfields(s)%output(ff,dd)))  ! have monthly forecast (ff) states to write?
     compute_monthly_aa = compute_monthly_aa .or. ( sfields(s)%output(aa,oo) .and. (.not. sfields(s)%output(aa,ee)) .and. (.not. sfields(s)%output(aa,dd)))  ! have monthly analysis (aa) states to write?
     compute_monthly_mm = compute_monthly_mm .or. ( sfields(s)%output(mm,oo) .and. (.not. sfields(s)%output(mm,ee)) .and. (.not. sfields(s)%output(mm,dd)))  ! have monthly daymean  (mm) states to write?
     compute_monthly_sf = compute_monthly_sf .or. ( sfields(s)%output(sf,oo) .and. (.not. sfields(s)%output(sf,ee)) .and. (.not. sfields(s)%output(sf,dd)))  ! have monthly forecast (sf) STD to write?
     compute_monthly_sa = compute_monthly_sa .or. ( sfields(s)%output(sa,oo) .and. (.not. sfields(s)%output(sa,ee)) .and. (.not. sfields(s)%output(sa,dd)))  ! have monthly analysis (sa) STD to write?
     compute_monthly_sm = compute_monthly_sm .or. ( sfields(s)%output(sm,oo) .and. (.not. sfields(s)%output(sm,ee)) .and. (.not. sfields(s)%output(sm,dd)))  ! have monthly daymean  (sm) STD to write?
  ENDDO

  IF ((.not. ALLOCATED(monthly_state_a)))     ALLOCATE(monthly_state_a(dim_p))
  IF ((.not. ALLOCATED(monthly_state_f)))     ALLOCATE(monthly_state_f(dim_p))
  IF ((.not. ALLOCATED(monthly_state_m)))     ALLOCATE(monthly_state_m(dim_p))
  
  IF ((.not. ALLOCATED(monthly_state_sa)))    ALLOCATE(monthly_state_sa(dim_p))
  IF ((.not. ALLOCATED(monthly_state_sf)))    ALLOCATE(monthly_state_sf(dim_p))
  IF ((.not. ALLOCATED(monthly_state_sm)))    ALLOCATE(monthly_state_sm(dim_p))
 
  ! allocate correction-counters at initial time; never de-allocated; reset to zero during each analysis
  IF (.not. ALLOCATED(count_lim_salt0_g))   ALLOCATE(count_lim_salt0_g  (dim_ens))
  IF (.not. ALLOCATED(count_lim_salt0_p))   ALLOCATE(count_lim_salt0_p  (dim_ens))
  IF (.not. ALLOCATED(count_lim_absvel_g))  ALLOCATE(count_lim_absvel_g (dim_ens))
  IF (.not. ALLOCATED(count_lim_absvel_p))  ALLOCATE(count_lim_absvel_p (dim_ens))
  IF (.not. ALLOCATED(count_lim_ssh_g))     ALLOCATE(count_lim_ssh_g    (dim_ens))
  IF (.not. ALLOCATED(count_lim_ssh_p))     ALLOCATE(count_lim_ssh_p    (dim_ens))
  IF (.not. ALLOCATED(count_lim_tempM2_g))  ALLOCATE(count_lim_tempM2_g (dim_ens))
  IF (.not. ALLOCATED(count_lim_tempM2_p))  ALLOCATE(count_lim_tempM2_p (dim_ens))

  ! initialize numbers
  invdim_ens = 1.0 / REAL(dim_ens)
  
  ! init monthly state
  IF ((step-step_null)==0) THEN
     monthly_state_a=  0.0D0
     monthly_state_m=  0.0D0
     monthly_state_f=  0.0D0
     monthly_state_sa= 0.0D0
     monthly_state_sm= 0.0D0
     monthly_state_sf= 0.0D0
  ENDIF
  ENDIF ! ((step-step_null)==0)

! ****************************
! *** Perform pre/poststep ***
! ****************************

  ! monthly event
  call monthly_event_assimstep(now_to_write_monthly)

! ****************************
! *** Corrections          ***
! ****************************
   
    IF ((step-step_null)<0) THEN
    ! *** store forecast state fields temporarily to compare with analysis afterwards ***    
      DO member = 1, dim_ens
        DO i = 1, dim_fields(id% SSH)
           state_fcst_SSH_p(i,member) = ens_p(i+offset(id% SSH),member)
        ENDDO
      ENDDO
    
    ELSE IF ((step-step_null)>0) THEN
   ! *** correcting assimilated state fields ***
   
   ! *** salinity must be > 0 ***
   count_lim_salt0_p = 0
   
   DO member = 1, dim_ens
      DO i = 1, dim_fields(id% salt)
           IF (ens_p(i+ offset(id% salt),member) < 0.0D0) THEN
               ens_p(i+ offset(id% salt),member) = 0.0D0
               
               count_lim_salt0_p(member) = count_lim_salt0_p(member)+1
           END IF
      END DO
   END DO
   
   CALL MPI_Allreduce(count_lim_salt0_p, count_lim_salt0_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   IF (mype_filter == 0) &
            WRITE(*, *) 'FESOM-PDAF', &
            '--- Updated salinity limited to zero: ', (count_lim_salt0_g(member), member = 1, dim_ens)
   
   ! *** SSH state update must be <= 2*sigma
   count_lim_ssh_p = 0
   
   DO member = 1, dim_ens
   DO i = 1, dim_fields(id% SSH)
       
       diffm = ens_p(offset(id% SSH)+i,member) - state_fcst_SSH_p(i,member)
   
       IF (ABS(diffm) > 2.0*stdev_SSH_f_p(i)) THEN
           ens_p(offset(id% SSH)+i,member) = state_fcst_SSH_p(i,member) + SIGN(2.0*stdev_SSH_f_p(i),diffm)
           count_lim_ssh_p(member) = count_lim_ssh_p(member)+1
       END IF
       
   END DO
   END DO
   
   CALL MPI_Allreduce(count_lim_ssh_p, count_lim_ssh_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   IF (mype_filter == 0) &
            WRITE(*,*) 'FESOM-PDAF', &
            '--- SSH updates limited to 2x standard deviation: ', (count_lim_ssh_g(member), member = 1, dim_ens)


   ! *** temperature must be > -2 degC ***
   count_lim_tempM2_p = 0
   
   DO member = 1, dim_ens
      DO i = 1, dim_fields(id% temp)
           IF (ens_p(i+ offset(id% temp),member) < -2.0) THEN
               ens_p(i+ offset(id% temp),member) = -2.0
               
               count_lim_tempM2_p(member) = count_lim_tempM2_p(member)+1
           END IF
      END DO
   END DO
   
   CALL MPI_Allreduce(count_lim_tempM2_p, count_lim_tempM2_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   IF (mype_filter == 0) &
            WRITE(*,*) 'FESOM-PDAF', &
            '--- Updated temperature limited to -2 degC: ', (count_lim_tempM2_g(member), member = 1, dim_ens)
            
   ! *** BGC fields must be larger than "tiny" ***
   IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- reset BGC to tiny'
   
   tiny_N   = tiny_chl/chl2N_max      ! Chl2N_max   = 0.00001/ 3.15d0 [mg CHL/mmol N] Maximum CHL-a:N ratio = 0.3 gCHL gN^{-1}
   tiny_N_d = tiny_chl/chl2N_max_d    ! Chl2N_max_d = 0.00001/ 4.2d0
   tiny_C   = tiny_N  /NCmax          ! NCmax       = 0.2d0           [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
   tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d     = 0.2d0 
   tiny_Si  = tiny_C_d/SiCmax         ! SiCmax      = 0.8d0
   tiny_R   = tiny * Redfield
   
   DO member = 1, dim_ens
    DO i = 1, myDim_nod2D
     DO k = 1, mesh_fesom%nlevels_nod2D(i)-1 ! loop through all wet nodes
     
           s = (i-1) * (nlmax) + k ! index in state vector
           
           DO f=1,nfields
             ! biogeochemical model tracers
             IF ((sfields(f)%bgc) .and. (sfields(f)% trnumfesom > 0)) THEN
                ens_p(offset(f)+s,member) = max(tiny,ens_p(offset(f)+s,member))
             ENDIF
           ENDDO ! f=1,nfields
           
           ! small phytoplankton
           ens_p(s+ offset(id% PhyN   ),member) = max(tiny_N  ,ens_p(s+ offset(id% PhyN   ),member))
           ens_p(s+ offset(id% PhyC   ),member) = max(tiny_C  ,ens_p(s+ offset(id% PhyC   ),member))
           ens_p(s+ offset(id% PhyChl ),member) = max(tiny_chl,ens_p(s+ offset(id% PhyChl ),member))
           ens_p(s+ offset(id% PhyCalc),member) = max(tiny    ,ens_p(s+ offset(id% PhyCalc),member))
           
           ! diatoms
           ens_p(s+ offset(id% DiaN   ),member) = max(tiny_N_d,ens_p(s+ offset(id% DiaN   ),member))
           ens_p(s+ offset(id% DiaC   ),member) = max(tiny_C_d,ens_p(s+ offset(id% DiaC   ),member))
           ens_p(s+ offset(id% DiaChl ),member) = max(tiny_chl,ens_p(s+ offset(id% DiaChl ),member))
           ens_p(s+ offset(id% DiaSi  ),member) = max(tiny_Si ,ens_p(s+ offset(id% DiaSi  ),member))
           
           ! zooplankton 1
           ens_p(s+ offset(id% Zo1N   ),member) = max(tiny    ,ens_p(s+ offset(id% Zo1N   ),member))
           ens_p(s+ offset(id% Zo1C   ),member) = max(tiny_R  ,ens_p(s+ offset(id% Zo1C   ),member))
           
           ! zooplankton 2
           ens_p(s+ offset(id% Zo2N   ),member) = max(tiny    ,ens_p(s+ offset(id% Zo2N   ),member))
           ens_p(s+ offset(id% Zo2C   ),member) = max(tiny_R  ,ens_p(s+ offset(id% Zo2C   ),member))

      ENDDO ! k=1,nlmax
    ENDDO ! i=1,my_Dim_nod2D
   ENDDO ! member=1,dim_ens

   END IF ! Corrections

! *******************************
! *** Compute ensemble mean   ***
! *******************************
  IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- compute ensemble mean'
  
  ! Local: 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i,member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


! *********************************************************************
! *** Store ensemble mean values for observation exclusion criteria ***
! *********************************************************************
! save values at forecast phase

  IF ((step-step_null)<0) THEN

     ! -- Sea-ice concentration --
     ! save mean_ice_p
     IF (    (sst_exclude_ice     .and. assim_o_sst     )&
        .OR. (sss_exclude_ice     .and. assim_o_sss     )&
        .OR. (sss_cci_exclude_ice .and. assim_o_sss_cci )&
        .OR. (chl_cci_exclude_ice .and. assim_o_chl_cci )) THEN 
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast SEA-ICE for observation exclusion'
        IF (ALLOCATED(mean_ice_p)) DEALLOCATE(mean_ice_p)
        ALLOCATE (mean_ice_p(dim_fields(id% a_ice)))
        mean_ice_p = state_p(offset(id% a_ice)+ 1 : &
                             offset(id% a_ice)+ dim_fields(id% a_ice))
     END IF

     ! -- SST --
     ! save mean_sst_p
     IF ((sst_exclude_diff > 0.0) .and. assim_o_sst) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast SST for observation exclusion'
        IF (ALLOCATED(mean_sst_p)) DEALLOCATE(mean_sst_p)
        ALLOCATE (mean_sst_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sst_p(i) = state_p(offset(id% temp) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! -- SSS (CASE SMOS) --
     ! save mean_sss_p
     IF ((sss_exclude_diff > 0.0) .and. assim_o_sss) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast SSS for observation exclusion'
        IF (ALLOCATED(mean_sss_p)) DEALLOCATE(mean_sss_p)
        ALLOCATE (mean_sss_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_p(i) = state_p(offset(id% salt) + (i-1) * (nlmax) + 1)
        END DO
     END IF

     ! -- SSS (CASE CCI) --
     ! save mean_sss_cci_p
     IF ((sss_cci_exclude_diff > 0.0) .and. assim_o_sss_cci) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast SSS for observation exclusion'
        IF (ALLOCATED(mean_sss_cci_p)) DEALLOCATE(mean_sss_cci_p)
        ALLOCATE (mean_sss_cci_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_cci_p(i) = state_p(offset(id% salt) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! -- Chlorophyll --
     ! save mean_chl_cci_p
     IF ((chl_cci_exclude_diff > 0.0) .and. assim_o_chl_cci) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast CHL for observation exclusion'
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
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean temperature (3D) for observation exclusion'
        ! Store mean temperature for profile assimilation
        IF (ALLOCATED(mean_temp_p)) DEALLOCATE(mean_temp_p)
        ALLOCATE (mean_temp_p(dim_fields(id%temp)))
        mean_temp_p = state_p(offset(id%temp)+1 : offset(id%temp)+dim_fields(id%temp))
     END IF
     
     ! -- oxygen --
     ! save mean_o2_p
     IF (((o2_merged_excl_absolute > 0.0) .or. (o2_merged_excl_relative > 0.0)) &
        .and. assim_o_o2_merged) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast oxygen for observation exclusion'
        IF (ALLOCATED(mean_O2_p)) DEALLOCATE(mean_O2_p)
        ALLOCATE (mean_O2_p(dim_fields(id%O2)))
        mean_O2_p = state_p(offset(id%O2)+1 : offset(id%O2)+dim_fields(id%O2))
     END IF
     
     ! -- nitrate --
     ! save mean_n_p
     IF (((n_merged_excl_absolute > 0.0) .or. (n_merged_excl_relative > 0.0)) &
        .and. assim_o_n_merged) THEN
        IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- save ensemble mean forecast DIN for observation exclusion'
        IF (ALLOCATED(mean_n_p)) DEALLOCATE(mean_n_p)
        ALLOCATE (mean_n_p(dim_fields(id%DIN)))
        mean_n_p = state_p(offset(id%DIN)+1 : offset(id%DIN)+dim_fields(id%DIN))
     END IF

  END IF ! forecast phase
  
! ************************************************
! *** Carbon sources minus sinks diagnostics   ***
! ************************************************
  
  ! factor to convert concentration to mass
  factor_mass = mesh_fesom%areasvol(:nlmax,:myDim_nod2D) * hnode_new(:nlmax,:myDim_nod2D) / SecondsPerDay
  factor_conc = 1.0 / SecondsPerDay

  IF ((step-step_null)<0) THEN
  ! forecast phase
  ! get fmass and fconc before analysis step
  
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- compute carbon diagnostics at forecast'
  
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
  
  IF ((step-step_null)>0) THEN
  ! analysis phase
  ! get amass and aconc after analysis step
  
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- compute carbon diagnostics at analysis'
  
     DO i = 1, myDim_nod2D
      DO k = 1, nlmax
      s = (i-1) * (nlmax) + k ! index in state vector
      ! DIC
      cffields(id_s_asml_dic)% aconc (k, i) = state_p(s + offset(id% DIC))
      ! Alk
      cffields(id_s_asml_alk)% aconc (k, i) = state_p(s + offset(id% Alk))
      ! Living carbon biomass
      cffields(id_s_asml_livingmatter)% aconc (k, i) = &
                               (state_p(s + offset(id% PhyC)) &
                              + state_p(s + offset(id% DiaC)) &
                              + state_p(s + offset(id% Zo1C)) &
                              + state_p(s + offset(id% Zo2C)) &
                              + state_p(s + offset(id% PhyCalc)))
      ! Dead organic carbon
      cffields(id_s_asml_deadmatter)% aconc (k, i) = &
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
       cffields(i)% amass = cffields(i)% aconc * factor_mass
       cffields(i)% aconc = cffields(i)% aconc * factor_conc
    ENDDO
    
    ! save the difference of forecast and analysis
    CALL carbonfluxes_diags_output_timemean_asml()
  
  ENDIF ! (analysis phase)
  

! *****************************************************************
! *** Compute ensemble spread (STD) for different fields        ***
! *****************************************************************
  
  ! Set debug output
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
  
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation'
  
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
      WRITE(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation surface mean'
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
      WRITE(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation global ocean mean'
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
          'FESOM-PDAF', 'Ensemble standard deviation:'
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
  
! *******************************
! *** Reset forgetting factor ***
! *******************************
! Forgetting factor increases after the start of the assimilation (days_since_DAstart).
! In case of model restarts:
!  -  days_since_DAstart is set by slurm-job-script
!  -  current forgetting factor and target temperature ensemble standard deviation are read from atmos-perturbation file
  
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


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

  IF (loctype==1 .AND. ((step-step_null) > 0)) THEN

     max_eff_dim_obs = 0.0
     min_eff_dim_obs = 1.0e16
     sum_eff_dim_obs = 0.0

     DO i = 1, myDim_nod2D
        IF (eff_dim_obs(i) > max_eff_dim_obs) max_eff_dim_obs = eff_dim_obs(i)
        IF (eff_dim_obs(i) < min_eff_dim_obs) min_eff_dim_obs = eff_dim_obs(i)
        sum_eff_dim_obs = sum_eff_dim_obs + eff_dim_obs(i)
     END DO
     IF (npes_filter>1) THEN
        CALL MPI_Reduce(sum_eff_dim_obs, avg_eff_dim_obs_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
             0, COMM_filter, MPIerr)
        CALL MPI_Reduce(max_eff_dim_obs, max_eff_dim_obs_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             0, COMM_filter, MPIerr)
        CALL MPI_Reduce(min_eff_dim_obs, min_eff_dim_obs_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
             0, COMM_filter, MPIerr)
     ELSE
        ! This is a work around for working with nullmpi.F90
        avg_eff_dim_obs_g = sum_eff_dim_obs
        min_eff_dim_obs_g = min_eff_dim_obs
        max_eff_dim_obs_g = max_eff_dim_obs
     END IF

     IF (mype_filter==0) THEN
        avg_eff_dim_obs_g = avg_eff_dim_obs_g / REAL(mesh_fesom%nod2d)

        WRITE (*, '(a, 8x, a)') &
             'FESOM-PDAF', '--- Effective observation dimensions for local analysis:'
        WRITE (*, '(a, 12x, a, f12.2)') &
             'FESOM-PDAF', 'min. effective observation dimension:       ', min_eff_dim_obs_g
        WRITE (*, '(a, 12x, a, f12.2)') &
             'FESOM-PDAF', 'max. effective observation dimension:       ', max_eff_dim_obs_g
        WRITE (*, '(a, 12x, a, f12.2)') &
             'FESOM-PDAF', 'avg. effective observation dimension:       ', avg_eff_dim_obs_g
     END IF
  END IF
  
! ***************************
! *** Compute daily means ***
! ***************************
! daily means ("m"-state) are averaged over one analysis step and the consecutive model forecast steps of that day
! during analysis step, add to m-fields
  IF (w_mm) THEN
  IF ((step-step_null) > 0) THEN
     timemean = timemean + state_p / delt_obs_ocn
  ENDIF ! step > 0
  ENDIF ! w_mm
  IF (w_sm) THEN
  IF ((step-step_null) > 0) THEN
     timemean_s = timemean_s + stdev_p / delt_obs_ocn
  ENDIF ! step > 0
  ENDIF ! w_sm
  
! *****************************
! *** Compute monthly means ***
! *****************************

  debugging_monthlymean = .false.
  
     ! include state into monthly mean
     IF ((step-step_null) > 0) THEN
     ! *** analyzed fields ***
       IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'adding to monthly mean.'
       IF (compute_monthly_aa) monthly_state_a  = monthly_state_a  + state_p
       IF (compute_monthly_sa) monthly_state_sa = monthly_state_sa + stdev_p
     ELSE IF ((step-step_null) < 0) THEN
     ! *** forecasted fields ***
       IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'adding to monthly mean.'
       IF (compute_monthly_ff) monthly_state_f  = monthly_state_f  + state_p
       IF (compute_monthly_sf) monthly_state_sf = monthly_state_sf + stdev_p
     END IF
     
     IF (now_to_write_monthly) THEN
     ! computing monthly mean at last day of month
     weights =  1.0/REAL(num_day_in_month(fleapyear,month))
     IF ((step-step_null) > 0) THEN
     ! *** analyzed state fields ***
       IF (compute_monthly_aa) monthly_state_a  = monthly_state_a  * weights
       IF (compute_monthly_sa) monthly_state_sa = monthly_state_sa * weights
     ELSE IF ((step-step_null) < 0) THEN
     ! *** forecasted state fields ***
       IF (compute_monthly_ff) monthly_state_f  = monthly_state_f  * weights
       IF (compute_monthly_sf) monthly_state_sf = monthly_state_sf * weights
     END IF
     ENDIF ! now_to_write_monthly

! **************************
! *** Write output files ***
! **************************
! note: after monthly output is written, reset monthly fields to zero

  ! *** write initial state fields ***
  IF ((step - step_null)==0 .and. ( .not. this_is_pdaf_restart)) THEN
      ! ensemble mean
      IF (w_dayensm) CALL netCDF_out('ii',state_p, int0, now_to_write_monthly, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g)
      IF (w_dayensm) CALL netCDF_out('si',stdev_p, int0, now_to_write_monthly)
      ! ensemble members
      IF (w_daymemb) THEN
        DO member = 1, dim_ens
           CALL netCDF_out('ii',ens_p(:,member), member, now_to_write_monthly)
        ENDDO
      ENDIF
  ENDIF

  ! daily output
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
  ELSE IF ((step-step_null) > 0) THEN
        ! *** write analysis ***
        ! ensemble mean
        IF (w_dayensm) CALL netCDF_out('aa',state_p   , int0, now_to_write_monthly, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g)
        IF (w_dayensm) CALL netCDF_out('sa',stdev_p   , int0, now_to_write_monthly)
        ! ensemble members
        IF (w_daymemb) THEN
          DO member = 1, dim_ens
            CALL netCDF_out('aa',ens_p(:,member), member, now_to_write_monthly)
          ENDDO
        ENDIF
  END IF
  ENDIF
  
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
  ELSE IF ((step-step_null) > 0) THEN
        ! *** write analysis and "m"-state fields ***
        ! ensemble mean, adding monthly mean of analysis
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('aa',state_p   , int0, now_to_write_monthly, stdev_surf_g=stdev_surf_g, stdev_volo_g=stdev_volo_g, m_state_p=monthly_state_a )
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('sa',stdev_p   , int0, now_to_write_monthly,                                                       m_state_p=monthly_state_sa)
        ! ensemble members, adding snapshot
        IF (w_daymemb .or. w_monmemb) THEN
          DO member = 1, dim_ens
            CALL netCDF_out('aa',ens_p(:,member), member, now_to_write_monthly, m_state_p=ens_p(:,member))
          ENDDO
        ENDIF
  END IF
  END IF

! at last day of month, reset monthly_state to zero (has been written)
IF (now_to_write_monthly) THEN
   IF ((step-step_null) > 0) THEN
   ! *** assimilated state fields ***
     IF (compute_monthly_aa) monthly_state_a  = 0.0D0
     IF (compute_monthly_sa) monthly_state_sa = 0.0D0
   ELSE IF ((step-step_null) < 0) THEN
   ! *** forecasted state fields ***
     IF (compute_monthly_ff) monthly_state_f  = 0.0D0
     IF (compute_monthly_sf) monthly_state_sf = 0.0D0
   END IF
ENDIF ! now_to_write_monthly

! ********************
! *** finishing up ***
! ********************

  ! variables deallocated after analysis step
  IF ((step-step_null) >= 0) THEN
     IF (allocated(stdev_SSH_f_p))    deallocate(stdev_SSH_f_p)
     IF (allocated(state_fcst_SSH_p)) deallocate(state_fcst_SSH_p)
     IF (allocated(mean_O2_p ))       deallocate(mean_O2_p )
     IF (allocated(mean_n_p))         deallocate(mean_n_p)
     IF (allocated(mean_chl_cci_p))   deallocate(mean_chl_cci_p)
     IF (allocated(mean_temp_p))      deallocate(mean_temp_p)
     IF (allocated(mean_sss_cci_p))   deallocate(mean_sss_cci_p)
     IF (allocated(mean_sss_p))       deallocate(mean_sss_p)     
     IF (allocated(mean_sst_p))       deallocate(mean_sst_p)     
     IF (allocated(mean_ice_p))       deallocate(mean_ice_p)     
  ENDIF
  DEALLOCATE(stdev_p)

END SUBROUTINE prepoststep_pdaf
