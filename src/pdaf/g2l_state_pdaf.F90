!$Id: g2l_state_pdaf.F90 2466 2021-02-25 12:46:34Z lnerger $
!> Routine to restrict a model state to a local analysis domain
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called during the loop over all
!! local analysis domains in the domain local filters
!! before the analysis on a single local analysis 
!! domain. It has to project the full process-local 
!! model state onto the current local analysis 
!! domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE g2l_state_pdaf(step, domain_p, dim_p, state_p, dim_l, state_l)

  USE mod_assim_pdaf, &           ! Variables for assimilation
       ONLY: id_lstate_in_pstate, ens_member_debug
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, mype_model
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  REAL, INTENT(in)    :: state_p(dim_p) !< PE-local full state vector 
  REAL, INTENT(out)   :: state_l(dim_l) !< State vector on local analysis domain
  
! *** Local variables *** 
  INTEGER :: i                    !< Counter
  INTEGER :: memberid             !< Ensemble member
  CHARACTER(LEN=17) :: filename   !< Filename for debugging output

! *************************************
! *** Initialize local state vector ***
! *************************************
 
  DO i = 1, dim_l
     state_l(i) = state_p(id_lstate_in_pstate(i))
  ENDDO
  
!~   IF ((mype_model==55) .AND. (domain_p==669)) THEN
!~         ens_member_debug = ens_member_debug + 1
!~         CALL PDAF_get_memberid(memberid)
!~         write (*,*) 'g2l Frauke-debug', memberid
!~         write (filename, "(A11,I2,A4)") "g2l_state_l", ens_member_debug, ".dat"
!~         open (ens_member_debug, file = TRIM(filename))
!~         write(ens_member_debug,*) state_l
!~         close(ens_member_debug) 
!~   END IF
  
  

END SUBROUTINE g2l_state_pdaf
