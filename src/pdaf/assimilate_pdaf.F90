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
       compute_monthly_mm, compute_monthly_sm, &
       assimilatePHY, assimilateBGC, &
       cda_phy, cda_bio, &
       type_trans
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
  INTEGER, INTENT(in) :: istep       ! current time step of model main loop

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
  
  ! work-around to update random seed
  REAL, ALLOCATABLE  :: dummyrndmat(:,:)
  INTEGER            :: nsteps_dummy
  INTEGER            :: next_assim_flag
  REAL               :: time_dummy

  ! External subroutines
  EXTERNAL :: collect_state_pdaf, &  ! Routine to collect a state vector from model fields
       distribute_state_pdaf, &      ! Routine to distribute a state vector to model fields
       next_observation_pdaf, &      ! Provide time step of next observation
       next_observation_pdaf_ncalls2_1, &      
       next_observation_pdaf_ncalls2_2, &      
       prepoststep_pdaf, &            ! User supplied pre/poststep routine
       prestep_pdaf, poststep_pdaf
  ! Localization of state vector
  EXTERNAL :: init_n_domains_pdaf, & ! Provide number of local analysis domains
       init_dim_l_pdaf, &            ! Initialize state dimension for local analysis domain
       init_dim_l_pdaf_PHY, &
       init_dim_l_pdaf_BGC, &
       g2l_state_pdaf, &             ! Get state on local analysis domain from global state
       l2g_state_pdaf                ! Update global state from state on local analysis domain
  ! Interface to PDAF-OMI for local and global filters
  EXTERNAL :: &
       init_dim_obs_pdafomi, &       ! Get dimension of full obs. vector for PE-local domain
       obs_op_pdafomi, &             ! Obs. operator for full obs. vector for PE-local domain
       init_dim_obs_l_pdafomi, &     ! Get dimension of obs. vector for local analysis domain
       init_dim_obs_pdafomi_PHY, &   !  """
       obs_op_pdafomi_PHY, &         !  """
       init_dim_obs_l_pdafomi_PHY, & !  """
       init_dim_obs_pdafomi_BGC, &   !  """
       obs_op_pdafomi_BGC, &         !  """
       init_dim_obs_l_pdafomi_BGC    !  """
  ! Subroutines used for generating observations
  EXTERNAL :: get_obs_f_pdaf         ! Get vector of synthetic observations from PDAF
  
  ! Variables for debugging:
  LOGICAL :: simplify_assim = .false.
  
  ! Possibility to reduce functionality for debugging
  simplify_assim = .true.

! *********************************
! *** Call assimilation routine ***
! *********************************

  ! init
  call daily_event  (IsLastStepDay,  1)
  call monthly_event(IsLastStepMonth,1)
  
  if (filterpe) allocate(dummyrndmat(dim_ens,dim_ens-1))

  ! istep:       Fesom's step:
  !              - starts at 1 at each model (re)start
  ! istep_asml:  imitates PDAF's step:
  !              - starts at 1 at beginning of each calendar year
  !              - starts at step_null at each (re)start

  if(mype_submodel==0 .and. task_id==1) write (*,'(a,1x,a,1x,a,1x,i5,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
          'FESOM-PDAF','assimilate_pdaf','istep', istep, 'istep_asml', istep_asml, 'day', daynew, 'time', FLOOR(timenew/3600.0),'h',INT(MOD(timenew,3600.0)/60.0),'min'

  ! Check  whether the filter is domain-localized
  CALL PDAF_get_localfilter(localfilter)

  ! Call assimilate routine
  IF (localfilter==1) THEN
     
     ! One or two calls to assimilation routine?
     IF ((assimilateBGC) .and. (assimilatePHY) .and. (trim(cda_phy)=='weak') .and. (trim(cda_bio)=='weak')) THEN
     ! Two consecutive calls for weakly coupled assimilation of PHY and BGC observations
     
        ! PHY assimilation (1)
        istep_asml = istep_asml + 1
        CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi_PHY, obs_op_pdafomi_PHY, prestep_pdaf, init_n_domains_pdaf, &
             init_dim_l_pdaf_PHY, init_dim_obs_l_pdafomi_PHY, g2l_state_pdaf, l2g_state_pdaf, &
             next_observation_pdaf_ncalls2_1, status_pdaf)
        ! BGC assimilation (2)
        istep_asml = istep_asml + 1
        CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
             init_dim_obs_pdafomi_BGC, obs_op_pdafomi_BGC, poststep_pdaf, init_n_domains_pdaf, &
             init_dim_l_pdaf_BGC, init_dim_obs_l_pdafomi_BGC, g2l_state_pdaf, l2g_state_pdaf, &
             next_observation_pdaf_ncalls2_2, status_pdaf)
     
     ELSE
     ! Assimilation routine is called once
        istep_asml = istep_asml + 1
        
        ! Select specific user supplied routines:   
        IF ((.not. assimilateBGC) .and. (assimilatePHY) .and. (trim(cda_phy)=='weak')) THEN
        ! One call for only-PHY assimilation
           
           CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
                init_dim_obs_pdafomi_PHY, obs_op_pdafomi_PHY, prepoststep_pdaf, init_n_domains_pdaf, &
                init_dim_l_pdaf_PHY, init_dim_obs_l_pdafomi_PHY, g2l_state_pdaf, l2g_state_pdaf, &
                next_observation_pdaf, status_pdaf)
           ! for consistency: work-around to update random seed after
           if (filterpe .and. (assim_flag==1)) CALL PDAF_seik_omega(dim_ens-1, dummyrndmat, type_trans, 1)
           
        ELSEIF ((assimilateBGC) .and. (.not. assimilatePHY) .and. (trim(cda_bio)=='weak')) THEN
        ! One call for only-BGC assimilation
           
           ! for consistency: work-around to update random seed before
           if (filterpe) then
              call next_observation_pdaf(istep_asml,nsteps_dummy,next_assim_flag,time_dummy)
              if (next_assim_flag==1) CALL PDAF_seik_omega(dim_ens-1, dummyrndmat, type_trans, 1)
           endif
           
           CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
                init_dim_obs_pdafomi_BGC, obs_op_pdafomi_BGC, prepoststep_pdaf, init_n_domains_pdaf, &
                init_dim_l_pdaf_BGC, init_dim_obs_l_pdafomi_BGC, g2l_state_pdaf, l2g_state_pdaf, &
                next_observation_pdaf, status_pdaf)
        
        ELSEIF ((assimilateBGC) .and. (assimilatePHY) .and. (trim(cda_phy)=='strong') .and. (trim(cda_bio)=='strong')) THEN
        ! Combined call for strongly coupled assimilation of PHY and BGC observations
           CALL PDAFomi_assimilate_local(collect_state_pdaf, distribute_state_pdaf, &
                init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, init_n_domains_pdaf, &
                init_dim_l_pdaf, init_dim_obs_l_pdafomi, g2l_state_pdaf, l2g_state_pdaf, &
                next_observation_pdaf, status_pdaf)
           ! for consistency: work-around to update random seed twice in total
           if (filterpe .and. (assim_flag==1)) CALL PDAF_seik_omega(dim_ens-1, dummyrndmat, type_trans, 1)
                
        ENDIF ! case-specific user supplied routines
     ENDIF ! one or two calls
     
  ELSE
!    IF (filtertype==11) THEN
!          ! Observation generation has its own OMI interface routine
!          CALL PDAFomi_generate_obs(collect_state_pdaf, distribute_state_pdaf, &
!               init_dim_obs_pdafomi, obs_op_pdafomi, get_obs_f_pdaf, &
!               prepoststep_pdaf, next_observation_pdaf, status_pdaf)
!    ELSE
!       ! All global filters except LEnKF
!       CALL PDAFomi_assimilate_global(collect_state_pdaf, distribute_state_pdaf, &
!            init_dim_obs_pdafomi, obs_op_pdafomi, prepoststep_pdaf, &
!            next_observation_pdaf, status_pdaf)
!    END IF
     WRITE (*,'(/a,i4,a1/)') &
          ' This code implementation is for LESTKF - stopping! (PE ', mype_world,')'
     CALL  abort_parallel()
  END IF

  ! Check for errors during execution of PDAF

  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a33,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in PDAF - stopping! (PE ', mype_world,')'
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
     
     IF (.not. simplify_assim) THEN ! simplify_assim-01
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
     ENDIF ! simplify_assim-01
     
     IF (.not. simplify_assim) THEN ! simplify_assim-02
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
     ENDIF ! simplify_assim-02
     
     ! *** write output and reset to zero ***
     IF (IsLastStepMonth) THEN
        IF (.not. simplify_assim) THEN ! simplify_assim-03
        ! monthly and daily output
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('mm',timemean  , int0, IsLastStepMonth, m_state_p=monthly_state_m )
        IF (w_dayensm .or. w_monensm) CALL netCDF_out('sm',timemean_s, int0, IsLastStepMonth, m_state_p=monthly_state_sm)
        ENDIF ! simplify_assim-03
        IF (.not. simplify_assim) THEN ! simplify_assim-04
        ! reset monthly and daily
        timemean = 0
        if (w_sm) timemean_s = 0
        if (compute_monthly_mm) monthly_state_m = 0
        if (compute_monthly_sm) monthly_state_sm = 0
        ENDIF ! simplify_assim-04
     ELSEIF (IsLastStepDay) THEN
        IF (.not. simplify_assim) THEN ! simplify_assim-05
        ! daily output
        IF (w_dayensm)            CALL netCDF_out('mm',timemean,   int0, IsLastStepMonth)
        IF (w_dayensm .and. w_sm) CALL netCDF_out('sm',timemean_s, int0, IsLastStepMonth)
        ENDIF ! simplify_assim-05
        IF (.not. simplify_assim) THEN ! simplify_assim-06
        ! reset daily
        timemean = 0
        if (w_sm) timemean_s = 0
        ENDIF ! simplify_assim-06
     ENDIF
     ENDIF ! filterpe
  ENDIF ! w_mm
  
  ! clean up
  if (filterpe) deallocate(dummyrndmat)
  
END SUBROUTINE assimilate_pdaf
