!! This file provides interface routines between the call-back routines
!! of PDAF and the observation-specific routines in PDAF-OMI. This structure
!! collects all calls to observation-specific routines in this single file
!! to make it easier to find the routines that need to be adapted.
!!
!! The routines here are mainly pure pass-through routines. Thus they
!! simply call one of the routines from PDAF-OMI. Partly some addtional
!! variable is required, e.g. to specify the offset of an observation
!! in the observation vector containing all observation types. These
!! cases are described in the routines.
!!
!! **Adding an observation type:**
!!   When adding an observation type, one has to add one module
!!   obs_TYPE_pdafomi (based on the template obs_TYPE_pdafomi_TEMPLATE.F90).
!!   In addition one has to add a call to the different routines include
!!   in this file. It is recommended to keep the order of the calls
!!   consistent over all files. 
!! 
!! __Revision history:__
!! * 2019-12 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!

!-------------------------------------------------------------------------------

!> Call-back routine for init_dim_obs
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_TYPE.
!!
SUBROUTINE init_dim_obs_pdafomi(step, dim_obs)

  ! Include functions for different observations
  USE mod_assim_pdaf, &
       ONLY:  proffiles_o, start_year_o, end_year_o, mype_debug, node_debug
  USE mod_parallel_pdaf, &
       ONLY: abort_parallel
  USE obs_sst_pdafomi, &
       ONLY: assim_o_sst, init_dim_obs_sst
  USE obs_sss_smos_pdafomi, &
       ONLY: assim_o_sss, init_dim_obs_sss
  USE obs_sss_cci_pdafomi, &
       ONLY: assim_o_sss_cci, init_dim_obs_sss_cci
  USE obs_ssh_cmems_pdafomi, &
       ONLY: assim_o_ssh, init_dim_obs_ssh
  USE obs_TSprof_EN4_pdafomi, &
       ONLY: assim_o_en4_t, assim_o_en4_s, init_dim_obs_prof
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs   !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_sst     ! Full number of SST observations
  INTEGER :: dim_obs_sss_cci ! Full number of SSS (CCI) observations
  INTEGER :: dim_obs_sss     ! Full number of SSS (SMOS) observations
  INTEGER :: dim_obs_ssh     ! Full number of SSH observations
  INTEGER :: dim_obs_prof    ! Full number of subsurface profile observations
  INTEGER :: dim_obs_en4ana  ! Full number of EN4 analysis profile observations


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_sst = 0
  dim_obs_sss = 0
  dim_obs_sss_cci = 0
  dim_obs_ssh = 0
  dim_obs_prof = 0
  dim_obs_en4ana = 0

  ! Call observation specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  
  ! No domain_p, thus no debugging call.
  
  IF (assim_o_sst)     CALL init_dim_obs_sst(step, dim_obs_sst)
  IF (assim_o_sss)     CALL init_dim_obs_sss(step, dim_obs_sss)
  IF (assim_o_sss_cci) CALL init_dim_obs_sss_cci(step, dim_obs_sss_cci)
  IF (assim_o_ssh)     CALL init_dim_obs_ssh(step, dim_obs_ssh)
  IF (assim_o_en4_t .OR. assim_o_en4_s) CALL init_dim_obs_prof(step, dim_obs_prof)

  dim_obs = dim_obs_sst + dim_obs_sss + dim_obs_sss_cci + dim_obs_ssh + dim_obs_prof + dim_obs_en4ana


  ! *** Generate profile observation files ***
  IF (proffiles_o == 1) THEN
     ! Generate distributed files
     CALL init_dim_obs_f_proffile_pdaf(start_year_o,end_year_o)
          
  ELSEIF (proffiles_o == 2) THEN
     ! Generate one global file
     WRITE(*,*) 'Generation of global file from EN4 raw profile data ', &
                'not yet implemented in FESOM Version 2.0 - stopping!'
     CALL abort_parallel
!~      CALL init_dim_obs_f_proffile_g_pdaf(step,dim_obs_prof)
  END IF

END SUBROUTINE init_dim_obs_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for obs_op
!!
!! This routine calls the observation-specific
!! routines obs_op_TYPE.
!!
SUBROUTINE obs_op_pdafomi(step, dim_p, dim_obs, state_p, ostate)

  ! Include functions for different observations
  USE obs_sst_pdafomi, ONLY: obs_op_sst
  USE obs_sss_smos_pdafomi, ONLY: obs_op_sss
  USE obs_sss_cci_pdafomi, ONLY: obs_op_sss_cci
  USE obs_ssh_cmems_pdafomi, ONLY: obs_op_ssh
  USE obs_TSprof_EN4_pdafomi, ONLY: obs_op_prof, thisobs
  USE mod_parallel_pdaf, ONLY: mype_filter

  USE mod_assim_pdaf, ONLY: mype_debug, node_debug
  USE PDAFomi, ONLY: PDAFomi_set_debug_flag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs              !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate(dim_obs)      !< PE-local full observed state

! No domain_p, thus no debugging call.

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  CALL obs_op_sst    (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_sss    (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_sss_cci(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_ssh    (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_prof   (dim_p, dim_obs, state_p, ostate)

END SUBROUTINE obs_op_pdafomi



!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs, dim_obs_l)

  ! Include observation types:
  USE obs_sst_pdafomi, ONLY: init_dim_obs_l_sst
  USE obs_sss_smos_pdafomi, ONLY: init_dim_obs_l_sss
  USE obs_sss_cci_pdafomi, ONLY: init_dim_obs_l_sss_cci
  USE obs_ssh_cmems_pdafomi, ONLY: init_dim_obs_l_ssh
  USE obs_TSprof_EN4_pdafomi, ONLY: init_dim_obs_l_prof

  ! General modules:
  USE PDAFomi, ONLY: PDAFomi_set_debug_flag
  USE mod_parallel_pdaf, ONLY: mype_filter
  USE g_parsup, ONLY: myList_nod2D
  USE mod_assim_pdaf, ONLY: debug_id_nod2, mype_debug, node_debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector
   
!~    ! Debugging:
!~    IF (mype_filter==mype_debug .AND. domain_p==node_debug) THEN
!~    CALL PDAFomi_set_debug_flag(domain_p)
!~    ELSE
!~    CALL PDAFomi_set_debug_flag(0)
!~    ENDIF


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

   CALL init_dim_obs_l_sst    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_sss    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_sss_cci(domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_ssh    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_prof   (domain_p, step, dim_obs, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdafomi
