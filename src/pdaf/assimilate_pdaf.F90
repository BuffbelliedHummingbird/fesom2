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
       dim_state_p, delt_obs_ocn, dim_ens
  USE g_clock, &
       ONLY: timenew

  IMPLICIT NONE
  include 'mpif.h'

! *** Arguments ***
  INTEGER, INTENT(in) :: istep       !< current time step

! *** Local variables ***
  INTEGER :: status_pdaf             ! PDAF status flag
  INTEGER :: localfilter             ! Flag for domain-localized filter (1=true)
  REAL, ALLOCATABLE :: state_p(:)    ! Ensemble member / mean state
  INTEGER :: mpierror

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

  istep_asml = istep + step_null

  if(mype_submodel==0 .and. task_id==1) write (*,'(a,i1.1,a,i,a,i,a,i)') &
          'FESOM ',task_id,' ',mype_submodel,' assimilate_pdaf, istep', istep, '  istep_asml', istep_asml

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
  
  ! Compute daily mean
  
  IF ( .not. ALLOCATED(state_p)) ALLOCATE(state_p(dim_state_p))
  
  ! 1. add snapshots to daily mean in between assimilation steps:
  IF (assim_flag == 0) THEN

        ! collect snapshot:
        CALL collect_state_pdaf(dim_state_p, state_p)

        ! add to daily mean:
        timemean = timemean + state_p / dim_ens / delt_obs_ocn

        ! compute ensemble mean on filter PE before assimilation step:
        IF (MOD(istep,delt_obs_ocn)==delt_obs_ocn-1) THEN
           IF (filterpe) THEN
              CALL MPI_REDUCE(MPI_IN_PLACE,timemean,dim_state_p,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_COUPLE,mpierror)
           ELSE
              CALL MPI_REDUCE(timemean,timemean,dim_state_p,MPI_DOUBLE_PRECISION,MPI_SUM,0,COMM_COUPLE,mpierror)
           ENDIF
        ENDIF

  ! 2. reset daily mean to zero after assimilation step:
  ELSEIF (assim_flag == 1) THEN
    timemean = 0.0
  ENDIF
  

END SUBROUTINE assimilate_pdaf
