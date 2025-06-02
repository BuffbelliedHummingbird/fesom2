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
       
  USE obs_chl_cci_pdafomi, &
       ONLY: assim_o_chl_cci, init_dim_obs_chl_cci
  USE obs_DIC_glodap_pdafomi, &
       ONLY: assim_o_DIC_glodap, init_dim_obs_DIC_glodap
  USE obs_Alk_glodap_pdafomi, &
       ONLY: assim_o_Alk_glodap, init_dim_obs_Alk_glodap
  USE obs_pco2_SOCAT_pdafomi, &
       ONLY: assim_o_pCO2_SOCAT, init_dim_obs_pCO2_SOCAT
  USE obs_O2_comf_pdafomi, &
       ONLY: assim_o_O2_comf, init_dim_obs_O2_comf
  USE obs_N_comf_pdafomi, &
       ONLY: assim_o_N_comf, init_dim_obs_N_comf
  USE obs_O2_argo_pdafomi, &
       ONLY: assim_o_O2_argo, init_dim_obs_O2_argo
  USE obs_N_argo_pdafomi, &
       ONLY: assim_o_N_argo, init_dim_obs_N_argo
  USE obs_O2_merged_pdafomi, &
       ONLY: assim_o_O2_merged, init_dim_obs_O2_merged
  USE obs_n_merged_pdafomi, &
       ONLY: assim_o_n_merged, init_dim_obs_n_merged
              
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
  
  INTEGER :: dim_obs_chl_cci ! Full number of biogeochem.-TYPE observations
  INTEGER :: dim_obs_DIC_glodap
  INTEGER :: dim_obs_Alk_glodap
  INTEGER :: dim_obs_pCO2_socat
  INTEGER :: dim_obs_O2_comf
  INTEGER :: dim_obs_N_comf
  INTEGER :: dim_obs_O2_argo
  INTEGER :: dim_obs_N_argo
  INTEGER :: dim_obs_O2_merged
  INTEGER :: dim_obs_n_merged

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
  dim_obs_chl_cci = 0
  dim_obs_DIC_glodap = 0
  dim_obs_Alk_glodap = 0
  dim_obs_pCO2_socat = 0
  dim_obs_O2_comf = 0
  dim_obs_N_comf = 0
  dim_obs_O2_argo = 0
  dim_obs_N_argo = 0
  dim_obs_O2_merged = 0
  dim_obs_n_merged = 0


  ! Call observation specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  
  ! No domain_p, thus no debugging call from here
  
  IF (assim_o_sst)     CALL init_dim_obs_sst(step, dim_obs_sst)
  IF (assim_o_sss)     CALL init_dim_obs_sss(step, dim_obs_sss)
  IF (assim_o_sss_cci) CALL init_dim_obs_sss_cci(step, dim_obs_sss_cci)
  IF (assim_o_ssh)     CALL init_dim_obs_ssh(step, dim_obs_ssh)
  IF (assim_o_en4_t .OR. assim_o_en4_s) CALL init_dim_obs_prof(step, dim_obs_prof)
  
  IF (assim_o_chl_cci)      CALL init_dim_obs_chl_cci(step, dim_obs_chl_cci)
  IF (assim_o_DIC_glodap)   CALL init_dim_obs_DIC_glodap(step, dim_obs_DIC_glodap)
  IF (assim_o_Alk_glodap)   CALL init_dim_obs_Alk_glodap(step, dim_obs_Alk_glodap)
  IF (assim_o_pCO2_SOCAT)   CALL init_dim_obs_pCO2_SOCAT(step, dim_obs_pCO2_socat)
  IF (assim_o_O2_comf)      CALL init_dim_obs_O2_comf(step, dim_obs_O2_comf)
  IF (assim_o_N_comf)       CALL init_dim_obs_N_comf(step, dim_obs_N_comf)
  IF (assim_o_O2_argo)      CALL init_dim_obs_O2_argo(step, dim_obs_O2_argo)
  IF (assim_o_N_argo)       CALL init_dim_obs_N_argo(step, dim_obs_N_argo)
  IF (assim_o_O2_merged)    CALL init_dim_obs_O2_merged(step, dim_obs_O2_merged)
  IF (assim_o_n_merged)     CALL init_dim_obs_n_merged(step, dim_obs_n_merged)

  dim_obs =   dim_obs_sst + dim_obs_sss + dim_obs_sss_cci + dim_obs_ssh + dim_obs_prof + dim_obs_en4ana &
            + dim_obs_chl_cci + dim_obs_DIC_glodap + dim_obs_Alk_glodap + dim_obs_pCO2_socat &
            + dim_obs_O2_comf + dim_obs_N_comf + dim_obs_o2_argo + dim_obs_N_argo &
            + dim_obs_O2_merged + dim_obs_n_merged


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
  USE obs_TSprof_EN4_pdafomi, ONLY: obs_op_prof
  
  USE obs_chl_cci_pdafomi, ONLY: obs_op_chl_cci
  USE obs_DIC_glodap_pdafomi, ONLY: obs_op_DIC_glodap
  USE obs_Alk_glodap_pdafomi, ONLY: obs_op_Alk_glodap
  USE obs_pco2_SOCAT_pdafomi, ONLY: obs_op_pCO2_SOCAT
  USE obs_O2_comf_pdafomi, ONLY: obs_op_O2_comf
  USE obs_N_comf_pdafomi, ONLY: obs_op_N_comf
  USE obs_O2_argo_pdafomi, ONLY: obs_op_o2_argo
  USE obs_N_argo_pdafomi, ONLY: obs_op_N_argo
  USE obs_O2_merged_pdafomi, ONLY: obs_op_o2_merged
  USE obs_n_merged_pdafomi, ONLY: obs_op_n_merged
  
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

! No domain_p, thus no debugging call from here

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
  
  CALL obs_op_chl_cci   (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_DIC_glodap(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_Alk_glodap(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_pCO2_SOCAT(dim_p, dim_obs, state_p, ostate)
  CALL obs_op_O2_comf   (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_N_comf    (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_O2_argo   (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_N_argo    (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_O2_merged (dim_p, dim_obs, state_p, ostate)
  CALL obs_op_n_merged  (dim_p, dim_obs, state_p, ostate)

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
  
  USE obs_chl_cci_pdafomi, ONLY: init_dim_obs_l_chl_cci
  USE obs_DIC_glodap_pdafomi, ONLY: init_dim_obs_l_DIC_glodap
  USE obs_Alk_glodap_pdafomi, ONLY: init_dim_obs_l_Alk_glodap
  USE obs_pco2_SOCAT_pdafomi, ONLY: init_dim_obs_l_pCO2_SOCAT
  USE obs_o2_comf_pdafomi, ONLY: init_dim_obs_l_o2_comf
  USE obs_n_comf_pdafomi, ONLY: init_dim_obs_l_n_comf
  USE obs_o2_argo_pdafomi, ONLY: init_dim_obs_l_o2_argo
  USE obs_n_argo_pdafomi, ONLY: init_dim_obs_l_n_argo
  USE obs_o2_merged_pdafomi, ONLY: init_dim_obs_l_o2_merged
  USE obs_n_merged_pdafomi, ONLY: init_dim_obs_l_n_merged

  ! General modules:
  USE PDAFomi, ONLY: PDAFomi_set_debug_flag
  USE mod_parallel_pdaf, ONLY: mype_filter
  USE g_parsup, ONLY: myList_nod2D, myDim_nod2D
  USE mod_assim_pdaf, ONLY: debug_id_nod2, mype_debug, node_debug

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs    !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector
   
   ! Debugging:
!  IF (mype_filter==mype_debug .AND. modulo(domain_p,myDim_nod2D)==node_debug) THEN
!   IF ( (mype_filter==57 .AND. modulo(domain_p,myDim_nod2D)==23)  .OR. &
!        (mype_filter==56 .AND. modulo(domain_p,myDim_nod2D)==129) .OR. &
!        (mype_filter==56 .AND. modulo(domain_p,myDim_nod2D)==802) .OR. &
!        (mype_filter==56 .AND. modulo(domain_p,myDim_nod2D)==880)      &
!        ) THEN
!    CALL PDAFomi_set_debug_flag(domain_p)
!    ELSE
!   CALL PDAFomi_set_debug_flag(0)
!   ENDIF



! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

   CALL init_dim_obs_l_sst    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_sss    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_sss_cci(domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_ssh    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_prof   (domain_p, step, dim_obs, dim_obs_l)
   
   CALL init_dim_obs_l_chl_cci   (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_DIC_glodap(domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_Alk_glodap(domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_pCO2_SOCAT(domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_o2_comf   (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_n_comf    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_o2_argo   (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_n_argo    (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_o2_merged (domain_p, step, dim_obs, dim_obs_l)
   CALL init_dim_obs_l_n_merged  (domain_p, step, dim_obs, dim_obs_l)

END SUBROUTINE init_dim_obs_l_pdafomi
