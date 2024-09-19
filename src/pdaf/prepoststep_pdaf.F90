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
       offset, proffiles_o, state_fcst, &
       monthly_state_f, monthly_state_a, write_monthly_mean,mon_snapshot_mem, &
       num_day_in_month, endday_of_month_in_year, startday_of_month_in_year, &
       depth_excl, depth_excl_no, this_is_pdaf_restart, mesh_fesom, nlmax, &
       dim_fields, dim_fields_glob, offset, offset_glob, nfields, id, &
       timemean, monthly_state_m, delt_obs_ocn, &
       monthly_state_ens_f, monthly_state_ens_a, days_since_DAstart, forget, &
       stdev_SSH_f_p
  USE mod_nc_out_variables, &
      ONLY: sfields
  USE mod_atmos_ens_stochasticity, &
      ONLY: stable_rmse
  USE g_PARSUP, &
       ONLY: MPI_DOUBLE_PRECISION, MPI_SUM, MPIerr, MPI_STATUS_SIZE, &
       MPI_INTEGER, MPI_MAX, MPI_MIN, mydim_nod2d, MPI_COMM_FESOM, &
       myList_edge2D, myDim_edge2D, myList_nod2D
  USE o_ARRAYS, ONLY: hnode_new
  USE recom_config, ONLY: tiny_chl, tiny, chl2N_max, chl2N_max_d, NCmax, &      
       NCmax_d, SiCmax, Redfield
  USE mod_nc_out_routines, &
       ONLY: netCDF_out, netCDF_STD_out
  USE mod_nc_out_variables, &
       ONLY: write_pos_da, write_ens, nfields_3D, ids_3D
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
  USE g_clock, &
        ONLY: dayold, yearold, check_fleapyr,daynew,yearnew
  USE mod_assim_pdaf, &
        ONLY: debug_id_nod2

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
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
  INTEGER :: i, j,k, member, ed, s  ! Counters
  INTEGER :: nlayermax              ! Max number of FESOM layers
  REAL :: invdim_ens                ! Inverse ensemble size
  REAL :: rmse_surf_p(nfields)             ! PE-local ensemble spread surface
  REAL :: rmse_surf_g(nfields)             ! Global ensemble spread surface

  REAL :: diffm                  ! temporary 
  CHARACTER(len=1) :: typestr       ! Character indicating call type
  REAL :: min_eff_dim_obs, max_eff_dim_obs       ! Stats on effective observation dimensions
  REAL :: min_eff_dim_obs_g, max_eff_dim_obs_g   ! Stats on effective observation dimensions
  REAL :: sum_eff_dim_obs, avg_eff_dim_obs_g     ! Stats on effective observation dimensions
  INTEGER :: dim_nod2D_ice_p        ! PE-local ice node dimension
  LOGICAL :: now_to_write_monthly, flag_print
  INTEGER :: month_iter, whichmonth,fleap,kk,row,nod
  REAL :: tmp_a, tmp_f
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
  
  REAL, ALLOCATABLE :: stdev_p(:)         ! ensemble standard deviation at grid proints
  REAL, ALLOCATABLE :: stdevprof_p(:)     ! ensemble standard deviation, horizontally-averaged profile
  REAL, ALLOCATABLE :: stdevprof_g(:)
  REAL, ALLOCATABLE :: meanprof_p(:)      ! mean state, horizontally-averaged profile
  REAL, ALLOCATABLE :: meanprof_g(:)
  REAL :: stdevglob_temp                       ! global 3D mean ensemble standard deviation for temperature field
  
  INTEGER :: offset_prof                  ! field offset in profile vector
  
  INTEGER, ALLOCATABLE, save :: dim_nodwet_p(:) ! number of wet nodes at depth levels
  INTEGER, ALLOCATABLE, save :: dim_nodwet_g(:)
  REAL, ALLOCATABLE, save  :: invdim_nodwet_p(:)
  REAL, ALLOCATABLE, save  :: invdim_nodwet_g(:)

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
     IF (step-step_null==0) THEN
       IF (.not.(this_is_pdaf_restart)) THEN
        WRITE (*,'(a, i7,3x,a)') 'FESOM-PDAF', step,'Analyze initial state ensemble'
        WRITE (typestr,'(a1)') 'i'
       ELSE
        WRITE (*,*) 'FESOM-PDAF: This is a PDAF restart. No initial fields are written.'
       END IF
     ELSE IF (step>0) THEN
        WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Analyze assimilated state ensemble'
        WRITE (typestr,'(a1)') 'a'
     ELSE IF (step<0) THEN
        WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Analyze forecast state ensemble'
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF ! IF (mype_filter==0)
  
  ! variables allocated and saved after forecast; and deallocated after analysis
  IF (.not. ALLOCATED(stdev_SSH_f_p)) ALLOCATE(stdev_SSH_f_p(dim_fields(id%SSH)))
  IF (.not. ALLOCATED(state_fcst)) ALLOCATE(state_fcst(dim_p,dim_ens))
  
  ! allocate monthly states at initial time; never de-allocated; monthly reset to zero
  IF (write_monthly_mean) THEN
     IF (.not. ALLOCATED(monthly_state_a))     ALLOCATE(monthly_state_a(dim_p))
     IF (.not. ALLOCATED(monthly_state_f))     ALLOCATE(monthly_state_f(dim_p))
     IF (.not. ALLOCATED(monthly_state_m))     ALLOCATE(monthly_state_m(dim_p))
     IF (.not. ALLOCATED(monthly_state_ens_a)) ALLOCATE(monthly_state_ens_a(dim_p,dim_ens))
     IF (.not. ALLOCATED(monthly_state_ens_f)) ALLOCATE(monthly_state_ens_f(dim_p,dim_ens))
  ENDIF
 
  ! allocate correction-counters at initial time; never de-allocated; reset to zero during each analysis
  IF (.not. ALLOCATED(count_lim_salt0_g))   ALLOCATE(count_lim_salt0_g  (dim_ens))
  IF (.not. ALLOCATED(count_lim_salt0_p))   ALLOCATE(count_lim_salt0_p  (dim_ens))
  IF (.not. ALLOCATED(count_lim_absvel_g))  ALLOCATE(count_lim_absvel_g (dim_ens))
  IF (.not. ALLOCATED(count_lim_absvel_p))  ALLOCATE(count_lim_absvel_p (dim_ens))
  IF (.not. ALLOCATED(count_lim_ssh_g))     ALLOCATE(count_lim_ssh_g    (dim_ens))
  IF (.not. ALLOCATED(count_lim_ssh_p))     ALLOCATE(count_lim_ssh_p    (dim_ens))
  IF (.not. ALLOCATED(count_lim_tempM2_g))  ALLOCATE(count_lim_tempM2_g (dim_ens))
  IF (.not. ALLOCATED(count_lim_tempM2_p))  ALLOCATE(count_lim_tempM2_p (dim_ens))

  ! Initialize numbers
  invdim_ens = 1.0 / REAL(dim_ens)

! ****************************
! *** Perform pre/poststep ***
! ****************************


   ! *** Correction ***
   
    IF ((step - step_null) < 0) THEN
    ! *** store forecast state fields temporarily to compare with analysis afterwards ***
    state_fcst = ens_p
    
    ELSE IF ((step - step_null) > 0) THEN
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
   DO i = offset(id% SSH)+1, offset(id% SSH)+dim_fields(id% SSH)
       
       diffm = ens_p(i,member) - state_fcst(i,member)
   
       IF (ABS(diffm) > 2.0*stdev_SSH_f_p(i-offset(id% SSH))) THEN
           ens_p(i,member) = state_fcst(i,member) + SIGN(2.0*stdev_SSH_f_p(i-offset(id% SSH)),diffm)
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
   
   tiny_N   = tiny_chl/chl2N_max      ! 0.00001/ 3.15d0   Chl2N_max [mg CHL/mmol N] Maximum CHL a : N ratio = 0.3 gCHL gN^-1
   tiny_N_d = tiny_chl/chl2N_max_d    ! 0.00001/ 4.2d0
   tiny_C   = tiny_N  /NCmax          ! NCmax   = 0.2d0   [mmol N/mmol C] Maximum cell quota of nitrogen (N:C)
   tiny_C_d = tiny_N_d/NCmax_d        ! NCmax_d = 0.2d0 
   tiny_Si  = tiny_C_d/SiCmax         ! SiCmax = 0.8d0
   tiny_R   = tiny * Redfield
   
   DO member = 1, dim_ens
      DO i = 1, dim_fields(id% PhyChl) ! (identical dimension for all 3D tracer fields)
           
           ens_p(i+ offset(id% PhyN   ),member) = max(tiny_N  ,ens_p(i+ offset(id% PhyN   ),member))
           ens_p(i+ offset(id% PhyC   ),member) = max(tiny_C  ,ens_p(i+ offset(id% PhyC   ),member))
           ens_p(i+ offset(id% PhyChl ),member) = max(tiny_chl,ens_p(i+ offset(id% PhyChl ),member))
           ens_p(i+ offset(id% PhyCalc),member) = max(tiny    ,ens_p(i+ offset(id% PhyCalc),member))
           
           ens_p(i+ offset(id% DetC   ),member) = max(tiny    ,ens_p(i+ offset(id% DetC   ),member))
           ens_p(i+ offset(id% DetN   ),member) = max(tiny    ,ens_p(i+ offset(id% DetN   ),member))
           ens_p(i+ offset(id% DetCalc),member) = max(tiny    ,ens_p(i+ offset(id% DetCalc),member))
           ens_p(i+ offset(id% DetSi  ),member) = max(tiny    ,ens_p(i+ offset(id% DetSi  ),member))
           
           ens_p(i+ offset(id% Det2C   ),member) = max(tiny    ,ens_p(i+ offset(id% Det2C   ),member))
           ens_p(i+ offset(id% Det2N   ),member) = max(tiny    ,ens_p(i+ offset(id% Det2N   ),member))
           ens_p(i+ offset(id% Det2Calc),member) = max(tiny    ,ens_p(i+ offset(id% Det2Calc),member))
           ens_p(i+ offset(id% Det2Si  ),member) = max(tiny    ,ens_p(i+ offset(id% Det2Si  ),member))
           
           ens_p(i+ offset(id% Zo1N   ),member) = max(tiny    ,ens_p(i+ offset(id% Zo1N   ),member))
           ens_p(i+ offset(id% Zo1C   ),member) = max(tiny_R  ,ens_p(i+ offset(id% Zo1C   ),member))
           
           ens_p(i+ offset(id% DOC    ),member) = max(tiny    ,ens_p(i+ offset(id% DOC    ),member))
           ens_p(i+ offset(id% DON    ),member) = max(tiny    ,ens_p(i+ offset(id% DON    ),member))
           
           ens_p(i+ offset(id% DiaN   ),member) = max(tiny_N_d,ens_p(i+ offset(id% DiaN   ),member))
           ens_p(i+ offset(id% DiaC   ),member) = max(tiny_C_d,ens_p(i+ offset(id% DiaC   ),member))
           ens_p(i+ offset(id% DiaChl ),member) = max(tiny_chl,ens_p(i+ offset(id% DiaChl ),member))
           ens_p(i+ offset(id% DiaSi  ),member) = max(tiny_Si ,ens_p(i+ offset(id% DiaSi  ),member))
           
           ens_p(i+ offset(id% O2     ),member) = max(tiny    ,ens_p(i+ offset(id% O2     ),member))
           
           ens_p(i+ offset(id% Zo2N   ),member) = max(tiny    ,ens_p(i+ offset(id% Zo2N   ),member))
           ens_p(i+ offset(id% Zo2C   ),member) = max(tiny_R  ,ens_p(i+ offset(id% Zo2C   ),member))
           
           ens_p(i+ offset(id% DIN    ),member) = max(tiny*1e-3,ens_p(i+ offset(id% DIN    ),member))
           ens_p(i+ offset(id% DIC    ),member) = max(tiny*1e-3,ens_p(i+ offset(id% DIC    ),member))
           ens_p(i+ offset(id% Alk    ),member) = max(tiny*1e-3,ens_p(i+ offset(id% Alk    ),member))
           
      END DO
   END DO

   END IF ! Corrections


  ! *** Compute mean state
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

  IF ((    assim_o_sst &
      .OR. assim_o_sss &
      .OR. assim_o_sss_cci &
      .OR. assim_o_chl_cci) &
      .AND. step<0) THEN
  ! (forecast phase)

   IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- save ensemble mean forecast for observation exclusion'

     ! sea-ice concentration
     IF (    sst_exclude_ice &
        .OR. sss_exclude_ice &
        .OR. sss_cci_exclude_ice &
        .OR. chl_cci_exclude_ice) THEN 
        IF (ALLOCATED(mean_ice_p)) DEALLOCATE(mean_ice_p)
        ALLOCATE (mean_ice_p(dim_fields(id% a_ice)))
        mean_ice_p = state_p(offset(id% a_ice)+ 1 : &
                             offset(id% a_ice)+ dim_fields(id% a_ice))
     END IF

     ! SST
     IF (sst_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_sst_p)) DEALLOCATE(mean_sst_p)
        ALLOCATE (mean_sst_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sst_p(i) = state_p(offset(id% temp) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! SSS SMOS
     IF (sss_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_sss_p)) DEALLOCATE(mean_sss_p)
        ALLOCATE (mean_sss_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_p(i) = state_p(offset(id% salt) + (i-1) * (nlmax) + 1)
        END DO
     END IF

     ! SSS CCI
     IF (sss_cci_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_sss_cci_p)) DEALLOCATE(mean_sss_cci_p)
        ALLOCATE (mean_sss_cci_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_cci_p(i) = state_p(offset(id% salt) + (i-1) * (nlmax) + 1)
        END DO
     END IF
     
     ! chl CCI
     IF (chl_cci_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_chl_cci_p)) DEALLOCATE(mean_chl_cci_p)
        ALLOCATE (mean_chl_cci_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_chl_cci_p(i) = state_p(offset(id% PhyChl) + (i-1) * (nlmax) + 1) &
                            + state_p(offset(id% DiaChl) + (i-1) * (nlmax) + 1)
        END DO
     END IF

  END IF

  ! 3D temperature field
  IF (        (assim_o_en4_t .OR. assim_o_en4_s) &
      .AND. &
              (step < 0) &
      .AND. &
              (prof_exclude_diff > 0.0)) THEN

     ! Store mean temperature for profile assimilation
     IF (ALLOCATED(mean_temp_p)) DEALLOCATE(mean_temp_p)
     ALLOCATE (mean_temp_p(myDim_nod2D))
     mean_temp_p = state_p(offset(id%temp)+1 : offset(id%temp)+dim_fields(id%temp))

  END IF


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
  ALLOCATE(stdev_p(dim_p))
  stdev_p(:) = 0.0
  DO member=1, dim_ens
     DO i=1, dim_p
        stdev_p(i) = stdev_p(i) & 
                + ((ens_p(i,member) - state_p(i)) * (ens_p(i,member) - state_p(i)))
     ENDDO ! j=1, dim_p
  ENDDO ! member=1, dim_ens
  stdev_p = SQRT(invdim_ens * stdev_p)
  where (stdev_p < 1.0e-10) stdev_p = 0.0 ! avoiding precision errors
  
  IF ((step - step_null) < 0) then ! if forecast: STD of SSH is saved and used at next analysis step
     stdev_SSH_f_p = stdev_p( offset(id%SSH)+1 : offset(id%SSH)+dim_fields(id%SSH) )
  endif
  
  
  IF (step-step_null==0) THEN
     ! Get wet-node count from model topography at initial step
     IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- get topography'
      
     ALLOCATE(dim_nodwet_p(nlmax))
     ALLOCATE(dim_nodwet_g(nlmax))
     ALLOCATE(invdim_nodwet_p(nlmax))
     ALLOCATE(invdim_nodwet_g(nlmax))
     
     dim_nodwet_p(:)=0
     DO i = 1, myDim_nod2D
        dim_nodwet_p(1:mesh_fesom%nlevels_nod2D(i)-1) = dim_nodwet_p(1:mesh_fesom%nlevels_nod2D(i)-1)+1
     ENDDO ! myDim_nod2D
     
     IF (write_debug) THEN
        DO k=1,nlmax
        write(*,'(1x,a14,1x,a20,1x,i2,1x,i8)') 'FESOM-PDAF', 'dim_nodwet_p on ', k, dim_nodwet_p(k)
        ENDDO
     ENDIF ! write_debug
     
     CALL MPI_Allreduce (dim_nodwet_p, dim_nodwet_g, nlmax, &
                         MPI_INTEGER, MPI_SUM, MPI_COMM_FESOM, &
                         MPIerr)
                         
     IF (write_debug) THEN
        DO k=1,nlmax
        write(*,'(1x,a14,1x,a20,1x,i2,1x,i8)') 'FESOM-PDAF', 'dim_nodwet_g on ', k, dim_nodwet_g(k)
        ENDDO
     ENDIF ! write_debug
     
     invdim_nodwet_p = 1.0/REAL(dim_nodwet_p)
     invdim_nodwet_g = 1.0/REAL(dim_nodwet_g)
     
  ENDIF ! step-step_null==0
  
  ! Compute pe-local surface mean of ensemble STD for each field
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation surface mean'
  rmse_surf_p = 0.0
  rmse_surf_g = 0.0
  
  DO j=1,nfields
     DO i = 1, myDim_nod2D
           ! STD:
           IF (sfields(j)%ndims == 2) THEN
           ! surface of 3D fields
           rmse_surf_p(j) =  rmse_surf_p(j) &
                               +  stdev_p( offset(j) + (i-1)*(nlmax) + 1 )
           ELSE
           ! surface fields
           rmse_surf_p(j) =  rmse_surf_p(j) &
                               +  stdev_p( offset(j) + i-1)
           ENDIF
     ENDDO ! i, myDim_nod2D
  ENDDO ! j, nfields
  where (rmse_surf_p < 1.0e-11) rmse_surf_p = 0.0 ! avoiding precision errors
  
  CALL MPI_Allreduce (rmse_surf_p, rmse_surf_g, nfields, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  rmse_surf_g = rmse_surf_g / REAL(mesh_fesom%nod2D)
  
  ! Compute pe-local horizontal mean of absolute state, and mean of ensemble STD, for 3D-fields
  IF (mype_filter == 0) &
      WRITE(*, *) 'FESOM-PDAF', '--- compute ensemble standard deviation profiles'
  ALLOCATE(stdevprof_p(nfields_3D*nlmax))
  ALLOCATE(stdevprof_g(nfields_3D*nlmax))
  ALLOCATE(meanprof_p (nfields_3D*nlmax))
  ALLOCATE(meanprof_g (nfields_3D*nlmax))
  
  IF (write_debug) THEN
       ! print mean and STD profiles on PE
       fileID_debug=30
       WRITE(day_string, '(i3.3)') daynew
       open(unit=fileID_debug, file='prof_p_'//day_string//'_'//typestr//'.txt', status='unknown')
  ENDIF ! write_debug
  
  stdevprof_p(:) = 0.0
  meanprof_p (:) = 0.0
  
  DO j=1,nfields_3D
     DO i = 1, myDim_nod2D
        DO k = 1, mesh_fesom%nlevels_nod2D(i)-1
           offset_prof = (j-1)*nlmax
           ! mean:
           meanprof_p (offset_prof + k) =  meanprof_p (offset_prof + k) &
                                        +  abs(state_p( (i-1)*(nlmax) + k + offset(ids_3D(j)) ))
           ! STD:
           stdevprof_p(offset_prof + k) =  stdevprof_p (offset_prof + k) &
                                        +  stdev_p( (i-1)*(nlmax) + k + offset(ids_3D(j)) )
        ENDDO ! k, nlevels
     ENDDO ! i, myDim_nod2D
  ENDDO ! j, nfields_3D
  where (stdevprof_p < 1.0e-11) stdevprof_p = 0.0 ! avoiding precision errors
  
  ! debugging mean of absolute values and STD profiles:
  IF (write_debug) THEN
    s = 1
    DO j=1,nfields_3D
       DO k=1,nlmax
         write(fileID_debug,'(a10,1x,i8,1x,G15.6,1x,G15.6)') &
                  sfields(ids_3D(j))%variable, s, &
                  meanprof_p(s) *invdim_nodwet_p(k), &
                  stdevprof_p(s)*invdim_nodwet_p(k)
         s = s+1
       ENDDO ! k, nlmax
    ENDDO ! j, nfields_3D
    close(fileID_debug)
  ENDIF ! write_debug
  
  ! Reduce to global horizontal mean of absolute state, and mean of ensemble STD, for 3D-fields
  CALL MPI_Allreduce (meanprof_p,  meanprof_g,  nfields_3D*nlmax, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
  CALL MPI_Allreduce (stdevprof_p, stdevprof_g, nfields_3D*nlmax, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)
       
  IF (write_debug) THEN
       ! print mean and STD profiles on PE
       fileID_debug=40
       WRITE(day_string, '(i3.3)') daynew
       open(unit=fileID_debug, file='prof_g_'//day_string//'_'//typestr//'.txt', status='unknown')
  ENDIF ! write_debug
  
  s = 1
  DO j=1,nfields_3D
    DO k=1,nlmax
         meanprof_g (s) = meanprof_g (s)*invdim_nodwet_g(k)
         stdevprof_g(s) = stdevprof_g(s)*invdim_nodwet_g(k)
         
         if (write_debug) write(fileID_debug,'(a10,1x,i8,1x,G15.6,1x,G15.6)') &
                  sfields(ids_3D(j))%variable, s, &
                  meanprof_g(s), &
                  stdevprof_g(s)
         s = s+1
     ENDDO ! k,nlmax
  ENDDO ! j,nfields_3D
  if (write_debug) close(fileID_debug)
  
  ! Global 3D mean for temperature field; temperature is 4th 3D-field (u,v,w,T,S,...)  
  stdevglob_temp = SUM( stdevprof_g((4-1)*nlmax+1 : (4-1)*nlmax+nlmax) )&
                   / REAL(nlmax)

  ! Display RMS errors
  IF (mype_filter==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'FESOM-PDAF', 'RMS error according to sampled covariance:'
     WRITE (*,'(a,7x,    a14, a14,a14,a14, a14,  a14,            /a, 10x,84a)') &
          'FESOM-PDAF', 'SSH','U','V','W','temp','salt', &
          'FESOM-PDAF', ('-',i=1,84)
     WRITE (*,'(a,10x,  6es14.4, 1x,a13,a1,/a, 10x,84a)')  &
          'FESOM-PDAF', rmse_surf_g(id% SSH), &
                        rmse_surf_g(id% u), &
                        rmse_surf_g(id% v), &
                        rmse_surf_g(id% w), &
                        rmse_surf_g(id% temp), &
                        rmse_surf_g(id% salt), &
                       'surface RMSe-', typestr, 'FESOM-PDAF', ('-',i=1,84)
     WRITE (*,'(a,10x,  6es14.4, 1x,a13,a1,/a, 10x,84a)')  &
          'FESOM-PDAF', 0.0, &
                        stdevprof_g((1-1)*nlmax+20), &
                        stdevprof_g((2-1)*nlmax+20), &
                        stdevprof_g((3-1)*nlmax+20), &
                        stdevprof_g((4-1)*nlmax+20), &
                        stdevprof_g((5-1)*nlmax+20), &
                       '450m    RMSe-', typestr, 'FESOM-PDAF', ('-',i=1,84)
  END IF
  
! *******************************
! *** Adapt forgetting factor ***
! *******************************
!
! Forgetting factor increases after the start of the assimilation (days_since_DAstart).
! In case of model restarts:
!    -  days_since_DAstart is set by slurm-job-script
!    -  current forgetting factor and "longterm" temperature ensemble spread is read from atmos-perturbation file
  
  IF (step<0) THEN
  ! forecast phase   
     
     IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' days_since_DAstart ', days_since_DAstart
     
     ! reset forgetting factor:
     IF     (days_since_DAstart==  1) THEN ! set value for 1st half-month of assimilation
       forget=0.95
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
       
     ELSEIF (days_since_DAstart== 16) THEN ! set value for 2nd half of 1st month of assimilation
       forget=0.96
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
       
     ELSEIF (days_since_DAstart== 32) THEN ! set value for 2nd month of assimilation
       forget=0.97
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Resetting forget to ', forget, ' at day ', days_since_DAstart
       
     ELSEIF (days_since_DAstart== 60) THEN ! set value for 3rd month of assimilation
       forget=0.98
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' resetting forget to ', forget, ' at day ', days_since_DAstart
       
     ELSEIF (days_since_DAstart== 90) THEN ! set value for 4th-17th months of assimilation
       forget=0.99
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' resetting forget to ', forget, ' at day ', days_since_DAstart
     
     ! during month 17, save temperature ensemble spread and from now on, use this RMSE at target value:
     ! (note: it is written with atmospheric perturbation to restart files)
     ELSEIF ((days_since_DAstart >= 485) .and. (days_since_DAstart <= 516)) THEN
       stable_rmse = stable_rmse + stdevglob_temp/31
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' Saving ensemble spread to adapt forget at day ', days_since_DAstart
       
     ! after month 17, the forgetting factor is reset whenever the ensemble spread becomes larger or smaller than target value:
     ELSEIF ((days_since_DAstart >= 516) .and. (stdevglob_temp > stable_rmse+0.0025)) THEN
       forget=1.00
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' resetting forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' stable RMSE: ', stable_rmse, &
                                      ' new forget: ',  forget
       
     ELSEIF ((days_since_DAstart >= 516) .and. (stdevglob_temp < stable_rmse-0.0025)) THEN
       forget=0.99
       CALL PDAF_reset_forget(forget)
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' resetting forget ',&
                                      'Current RMSE: ', stdevglob_temp,&
                                      ' target RMSE: ', stable_rmse,&
                                      ' new forget: ', forget
       
     ELSE
       IF (mype_filter==0) write(*,*) 'FESOM-PDAF', ' not resetting forget ',&
                                      'Current RMSE: ', stdevglob_temp, &
                                      ' target RMSE: ', stable_rmse,&
                                      'forget', forget
     
     ENDIF ! (whether to reset forgetting factor)
     
     ! count upwards during daily forecast phase:
     days_since_DAstart=days_since_DAstart+1
     
  ENDIF ! (forecast phase)


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

  IF (loctype==1 .AND. ((step - step_null) > 0)) THEN

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

  IF ((step - step_null) > 0) THEN
     timemean = timemean + state_p / delt_obs_ocn
  ENDIF

! *****************************
! *** Compute monthly means ***
! *****************************

  debugging_monthlymean = .false.

  IF (write_monthly_mean) THEN
  
      now_to_write_monthly = .FALSE.
      call check_fleapyr(yearold, fleap)
      
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'dayold: ', dayold
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'daynew: ', daynew
      
      ! set "now_to_write_monthly" at last day of month
      find_month: DO month_iter=1,12
        if (dayold <= endday_of_month_in_year(fleap,month_iter)) THEN
            whichmonth = month_iter
            if (dayold == endday_of_month_in_year(fleap,month_iter)) THEN
               now_to_write_monthly = .TRUE.
            endif
            exit find_month
        endif
      ENDDO find_month
      
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'whichmonth: ', whichmonth
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'now_to_write_monthly: ', now_to_write_monthly

      ! reset monthly_state to zero at first day of month
      monthly_state_p_0: DO month_iter=1,12
        if (dayold == startday_of_month_in_year(fleap,month_iter)) THEN
          IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'resetting monthly mean to zero.'
          IF ((step - step_null) > 0) THEN
          ! *** assimilated state fields ***
            monthly_state_a= 0.0D0
            monthly_state_m= 0.0D0
            IF (write_ens) monthly_state_ens_a= 0.0D0
          ELSE IF ((step - step_null) < 0) THEN
          ! *** forecasted state fields ***
            monthly_state_f= 0.0D0
            IF (write_ens) monthly_state_ens_f= 0.0D0
          END IF
          exit monthly_state_p_0
        endif
      ENDDO monthly_state_p_0
      
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'num_day_in_month: ', num_day_in_month(fleap,whichmonth)
      
      ! include state into monthly mean
      IF ((step - step_null) > 0) THEN
      ! *** analyzed state fields ***
        IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'adding to monthly mean.'
        monthly_state_a = monthly_state_a + state_p /num_day_in_month(fleap,whichmonth)
        monthly_state_m = monthly_state_m + timemean/num_day_in_month(fleap,whichmonth)
        IF (write_ens) monthly_state_ens_a = monthly_state_ens_a + ens_p/num_day_in_month(fleap,whichmonth)
      ELSE IF ((step - step_null) < 0) THEN
      ! *** forecasted state fields ***
        IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'adding to monthly mean.'
        monthly_state_f = monthly_state_f + state_p/num_day_in_month(fleap,whichmonth)
        IF (write_ens) monthly_state_ens_f = monthly_state_ens_f + ens_p/num_day_in_month(fleap,whichmonth)
      END IF

   ENDIF
   
! **************************
! *** Write output files ***
! **************************

! *** write initial state fields ***
IF ((step - step_null)==0 .and. ( .not. this_is_pdaf_restart)) THEN
      CALL netCDF_out('i',daynew,step,state_p,ens_p,rmse_surf_g,forget)
      CALL netCDF_STD_out('i',daynew,step,stdev_p,stdevprof_g,meanprof_g)
ENDIF

IF (.not. write_monthly_mean) THEN
  ! daily output
  write_pos_da = daynew
  IF ((step - step_null) < 0) THEN
        ! *** write forecast state fields ***
        CALL netCDF_out('f',daynew,step,state_p,ens_p,rmse_surf_g,forget)
        CALL netCDF_STD_out('f',daynew,step,stdev_p,stdevprof_g,meanprof_g)
  ELSE IF ((step - step_null) > 0) THEN
        ! *** write analysis state fields ***
        CALL netCDF_out('a',daynew,step,state_p, ens_p,rmse_surf_g,forget)
        CALL netCDF_out('m',daynew,step,timemean,ens_p,rmse_surf_g,forget) ! Warning: ensemble data for day-averages ('m') not yet available.
        CALL netCDF_STD_out('a',daynew,step,stdev_p,stdevprof_g,meanprof_g)
  END IF

ELSEIF (write_monthly_mean .and. now_to_write_monthly) THEN
  ! monthly output
  write_pos_da = whichmonth
  IF ((step - step_null) < 0) THEN
        ! *** write forecast state fields ***
!~         CALL netCDF_out('f',write_pos_da,step,monthly_state_f,monthly_state_ens_f,rmse_surf_g,forget)  ! writes a monthly mean
!~         CALL netCDF_STD_out('f',write_pos_da,step,stdev_p,stdevprof_g,meanprof_g)               ! writes a monthly snapshot
  ELSE IF ((step - step_null) > 0) THEN
        ! *** write analysis state fields ***
!~         CALL netCDF_out('a',write_pos_da,step,monthly_state_a,monthly_state_ens_a,rmse_surf_g,forget) ! writes a monthly mean
!~         CALL netCDF_out('m',write_pos_da,step,monthly_state_m,monthly_state_ens_a,rmse_surf_g,forget) ! Warning: ensemble data for day-averages ('m') not yet available.
!~         CALL netCDF_STD_out('a',write_pos_da,step,stdev_p,stdevprof_g,meanprof_g)              ! writes a monthly snapshot
  END IF
END IF ! write_monthly_mean

! ********************
! *** finishing up ***
! ********************

  IF ((step - step_null) >= 0) DEALLOCATE(stdev_SSH_f_p,state_fcst)     ! if forecast, keep stdev_SSH_f_p; it is used at next analysis step.
  DEALLOCATE(stdevprof_p, stdevprof_g, meanprof_p, meanprof_g,stdev_p)

END SUBROUTINE prepoststep_pdaf
