!!  Routine to call PDAF for analysis step
!!
!! This routine is called during the model integrations at each time 
!! step. It calls the filter-specific assimilation routine of PDAF 
!! (PDAF_assimilate_X), which checks whether the forecast phase is
!! completed. If so, the analysis step is computed inside PDAF
!!
!! __Revision history:__
!! * 2013-08 - Lars Nerger - Initial code for AWI-CM

SUBROUTINE assimilate_pdaf(istep)

  USE pdaf_interfaces_module, &   ! Interface definitions to PDAF core routines
       ONLY: PDAFomi_assimilate_local, PDAFomi_assimilate_global, &
       PDAFomi_assimilate_lenkf, PDAFomi_generate_obs, PDAF_get_localfilter
  USE PDAF_mod_filter, &
       ONLY: cnt_steps, assim_flag
  USE mod_parallel_pdaf, &        ! Parallelization variables
       ONLY: mype_world, abort_parallel, task_id, mype_submodel, &
       COMM_COUPLE, filterpe
  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: filtertype, istep_asml, step_null, timemean, &
       dim_state_p, delt_obs_ocn, dim_ens, timemean_s, &
       monthly_state_sm, monthly_state_m, &
       compute_monthly_mm, compute_monthly_sm
  USE mod_nc_out_variables, &
       ONLY: w_mm, w_sm, w_dayensm, w_monensm
  USE mod_nc_out_routines, &
       ONLY: netCDF_out
  USE g_clock, &
       ONLY: timenew, daynew, yearnew, month, &
       num_day_in_month, fleapyear
  USE g_events, &
       ONLY: daily_event, monthly_event

  IMPLICIT NONE
  include 'mpif.h'

! *** Arguments ***
  INTEGER, INTENT(in) :: istep       !< current time step

! *** Local variables ***
  INTEGER :: status_pdaf             ! PDAF status flag
  INTEGER :: localfilter             ! Flag for domain-localized filter (1=true)
  REAL, ALLOCATABLE :: state_p(:)    ! Ensemble member / mean state
  REAL, ALLOCATABLE :: stdev_p(:)    ! Standard deviation
  REAL, ALLOCATABLE :: ensm_p(:)     ! Ensemble mean state
  INTEGER :: mpierror
  real :: invdim_ens                 ! Inverse ensemble size
  real :: weights
  
  logical :: IsLastStepDay
  logical :: IsLastStepMonth
  
  INTEGER, parameter :: int0 = 0

  ! External subroutines
  EXTERNAL :: collect_state_pdaf, &  ! Routine to collect a state vector from model fields
       distribute_state_pdaf, &      ! Routine to distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       prepoststep_pdaf              ! User supplied pre/poststep routine
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi        ! Get dimension of obs. vector for local analysis domain
  ! Subroutines used for generating observations
  EXTERNAL :: get_obs_f_pdaf         ! Get vector of synthetic observations from PDAF

! *********************************
! *** Call assimilation routine ***
! *********************************

  call daily_event  (IsLastStepDay,  1)
  call monthly_event(IsLastStepMonth,1)

  istep_asml = istep + step_null  ! istep:       starting at 1 at each model (re)start
                                  ! istep_asml:  starting at 1 at beginning of each year

  if(mype_submodel==0 .and. task_id==1) write (*,'(a,1x,a,1x,a,1x,i5,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
          'FESOM-PDAF','assimilate_pdaf','istep', istep, 'istep_asml', istep_asml, 'day', daynew, 'time', FLOOR(timenew/3600.0),'h',INT(MOD(timenew,3600.0)/60.0),'min'

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine for global or local filter
  IF (localfilter==1) THEN
     CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
          init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
          init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
          next_observation_pdaf, status_pdaf)
  ELSE
     IF (filtertype==11) THEN
!~         ! Observation generation has its own OMI interface routine
!~         CALL PDAFomi_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
!~              init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_f_pdaf, &
!~              prepoststep_pdaf, next_observation_pdaf, status_pdaf)
     ELSE
        ! All global filters except LEnKF
        CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
             next_observation_pdaf, status_pdaf)
     END IF
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF_put_state - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF
  
  ! *********************************
  ! *** Compute daily mean        ***
  ! *********************************
  ! daily means ("m"-state) are averaged over 1 analysis step followed by step_per_day-minus-1 model forecast steps

  ! note: computing ensemble mean of state fields at every step is less efficient,
  ! but required to compute ensemble standard deviation at every step
     
  IF (w_mm .or. w_sm) THEN
     IF ( .not. ALLOCATED(state_p))            ALLOCATE(state_p(dim_state_p))
     IF ( .not. ALLOCATED(ensm_p ))            ALLOCATE(ensm_p (dim_state_p))
     IF ( .not. ALLOCATED(stdev_p) .and. w_sm) ALLOCATE(stdev_p(dim_state_p))
     
     ! *** in between assimilation steps, add forecast steps to m-fields ***
     ! note: assimilation step is at first time step of day
     IF (assim_flag == 0) THEN
        ! collect instantenous model data
        CALL collect_state_pdaf(dim_state_p, state_p)
        ! compute ensemble mean
        invdim_ens = 1.0 / REAL(dim_ens)
        CALL MPI_ALLREDUCE((state_p*invdim_ens),ensm_p,dim_state_p,MPI_DOUBLE_PRECISION,MPI_SUM,COMM_COUPLE,mpierror)
        ! compute ensemble mean of squared deviations on filterpe
        IF (w_sm) CALL MPI_REDUCE(((ensm_p-state_p)*(ensm_p-state_p)*invdim_ens),stdev_p,dim_state_p,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_COUPLE,mpierror)
        ! add to daily mean on filterpe
        IF (filterpe) then
           timemean    = timemean    + ensm_p / delt_obs_ocn
           ! add daily mean of standard deviation
           IF (w_sm) stdev_p     = SQRT(invdim_ens * stdev_p)
           IF (w_sm) timemean_s  = timemean_s  + stdev_p / delt_obs_ocn
        ENDIF ! filterpe
     ENDIF ! assim_flag
     
     ! *** compute monthly means ***
     IF (filterpe) THEN
     IF (IsLastStepDay) THEN
        ! add m-fields to monthly mean
        IF (compute_monthly_mm) monthly_state_m  = monthly_state_m  + timemean
        IF (compute_monthly_sm) monthly_state_sm = monthly_state_sm + timemean_s
     ENDIF
     IF (IsLastStepMonth) THEN
        ! compute monthly mean
        weights = 1.0/REAL(num_day_in_month(fleapyear,month))
        IF (compute_monthly_mm) monthly_state_m  = monthly_state_m  * weights
        IF (compute_monthly_sm) monthly_state_sm = monthly_state_sm * weights
     ENDIF
     
     ! *** write output and reset to zero ***
     IF (IsLastStepMonth) THEN
        ! monthly and daily output
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('mm',timemean  , int0, IsLastStepMonth, m_state_p=monthly_state_m )
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('sm',timemean_s, int0, IsLastStepMonth, m_state_p=monthly_state_sm)
        ! reset monthly and daily
        timemean = 0
        if (w_sm) timemean_s = 0
        if (compute_monthly_mm) monthly_state_m = 0
        if (compute_monthly_sm) monthly_state_sm = 0
     ELSEIF (IsLastStepDay) THEN
        ! daily output
        IF (w_dayensm)            CALL netCDF_out('mm',timemean,   int0, IsLastStepMonth)
        IF (w_dayensm .and. w_sm) CALL netCDF_out('sm',timemean_s, int0, IsLastStepMonth)
        ! reset daily
        timemean = 0
        if (w_sm) timemean_s = 0
     ENDIF
     ENDIF ! filterpe
  ENDIF ! w_mm

END SUBROUTINE assimilate_pdaf
