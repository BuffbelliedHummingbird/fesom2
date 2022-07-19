!$Id: l2g_state_pdaf.F90 2466 2021-02-25 12:46:34Z lnerger $
!>  Routine to initialize full state from local analysis
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called during the loop over all
!! local analysis domains in the domain local filters
!! after the analysis and ensemble transformation 
!! on a single local analysis domain. It has to 
!! initialize elements of the PE-local full state 
!! vector from the provided analysis state vector 
!! on the local analysis domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE l2g_state_pdaf(step, domain, dim_l, state_l, dim_p, state_p)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: id_lstate_in_pstate, ens_member_debug
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, mype_model
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain         !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< Process-local full state dimension
  REAL, INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL, INTENT(inout) :: state_p(dim_p) !< Process-local full state vector 
  
! *** Local variables *** 
  INTEGER :: i                    !< Counter
  INTEGER :: memberid             !< Ensemble member
  CHARACTER(LEN=17) :: filename   !< Filename for debugging output
  
! **************************************************
! *** Initialize elements of global state vector ***
! **************************************************

  DO i = 1, dim_l
     state_p(id_lstate_in_pstate(i)) = state_l(i)
  END DO
  
!~   IF ((mype_model==55) .AND. (domain==669)) THEN
!~         ens_member_debug = ens_member_debug + 1
!~         CALL PDAF_get_memberid(memberid)
!~         write (*,*) 'l2g Frauke-debug', memberid
!~         write (filename, "(A11,I2,A4)") "l2g_state_l", ens_member_debug, ".dat"
!~         open (20+ens_member_debug, file = TRIM(filename))
!~         write(20+ens_member_debug,*) state_l
!~         close(20+ens_member_debug) 
!~   END IF

END SUBROUTINE l2g_state_pdaf
