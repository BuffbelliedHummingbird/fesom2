!$Id: init_n_domains_pdaf.F90 2466 2021-02-25 12:46:34Z lnerger $
!>  Routine to set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! The routine is called in PDAF_X_update 
!! or PDAF_assimialte_X at the beginning of the 
!! analysis step before the loop through all local
!! analysis domains. It has to set the number of 
!! local analysis domains for the process-local domain.
!!
!! The routine is called by all filter processes.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code for AWI-CM
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)

  USE g_parsup, &
      ONLY: myDim_nod2D, &      ! Process-local number of vertices
            myDim_elem2D        ! Process-local number of elements
  USE mod_parallel_pdaf, &
      ONLY: mype_filter
  USE PDAFomi, &
      ONLY: PDAFomi_set_debug_flag

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        ! Current time step
  INTEGER, INTENT(out) :: n_domains_p ! Process-local number of analysis domains


! ************************************
! *** Initialize number of domains ***
! ************************************

  n_domains_p = myDim_nod2D

END SUBROUTINE init_n_domains_pdaf
