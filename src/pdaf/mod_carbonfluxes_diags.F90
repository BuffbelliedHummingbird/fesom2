MODULE mod_carbon_fluxes_diags


IMPLICIT NONE
 
   logical            :: cfdiags_debug ! whether to write debug output for tracer mass conservation
   logical            :: cfconc        ! whether to write concentration or mass fields
   character(len=100) :: vname_cfdiags ! variable name for debugging output


!  ********************************************
!  *** Diagnostic variables (instantaneous) ***
!  ********************************************
!  dimension: (nlay, myDim_nod2D)

!  _______________________________________________
!  *** Transport of carbon [mmol C / sec / m2] ***
!  positive: upward/eastward/northward ___________
   
!  alkalinity
   real, allocatable :: t_u_alk            (:,:) ! advective transport
   real, allocatable :: t_v_alk            (:,:)
   real, allocatable :: t_w_alk            (:,:)

!  DIC
   real, allocatable :: t_u_dic            (:,:) ! advective transport
   real, allocatable :: t_v_dic            (:,:)
   real, allocatable :: t_w_dic            (:,:)
   
!  Living biomass (Phy, Dia, Het, Zoo)
   real, allocatable :: t_u_livingmatter   (:,:) ! advective transport
   real, allocatable :: t_v_livingmatter   (:,:)
   real, allocatable :: t_w_livingmatter   (:,:)
   real, allocatable :: t_sink_livingmatter(:,:) ! sinking in the water column

!  Dead biomass (Detritus and DOC)
   real, allocatable :: t_u_deadmatter     (:,:) ! advective transport of detritus and DOC
   real, allocatable :: t_v_deadmatter     (:,:)
   real, allocatable :: t_w_deadmatter     (:,:)
   real, allocatable :: t_export           (:,:) ! sinking of detritus in the water column

!  _____________________________________________________   
!  *** Sources / sinks of carbon [mmol C / m3 / sec] ***

!  alkalinity
   real, allocatable :: s_hor_alk           (:,:)  ! advection
   real, allocatable :: s_ver_alk           (:,:)
   real, allocatable :: s_diffV_alk         (:,:)  ! vertical diffusion
   real, allocatable :: s_diffH_alk         (:,:)  ! horizontal diffusion
   real, allocatable :: s_bio_alk           (:,:)  ! remineralization and calcification

!  DIC
   real, allocatable :: s_hor_dic           (:,:)  ! advection
   real, allocatable :: s_ver_dic           (:,:)
   real, allocatable :: s_diffV_dic         (:,:)  ! diffusion
   real, allocatable :: s_diffH_dic         (:,:)
   real, allocatable :: s_bio_dic           (:,:)  ! source of DIC from deadmatter (remineralization of DOC; calcite dissolution)

!  Living biomass (Phy, Dia, Het, Zoo)
   real, allocatable :: s_hor_livingmatter  (:,:)  ! advection
   real, allocatable :: s_ver_livingmatter  (:,:)
   real, allocatable :: s_diffV_livingmatter(:,:)  ! diffusion
   real, allocatable :: s_diffH_livingmatter(:,:)
   real, allocatable :: s_sink_livingmatter (:,:)  ! sinking in the water column
   real, allocatable :: s_bio_livingmatter  (:,:)  ! net source of alive biomass from DIC (photosynthesis - respiration; calcification)

!  Dead organic carbon (Detritus and DOC)
   real, allocatable :: s_hor_deadmatter    (:,:)  ! advection of detritus and DOC
   real, allocatable :: s_ver_deadmatter    (:,:)
   real, allocatable :: s_diffV_deadmatter  (:,:)  ! diffusion 
   real, allocatable :: s_diffH_deadmatter  (:,:)
   real, allocatable :: s_export            (:,:)  ! sinking of detritus in the water column
   real, allocatable :: s_bio_deadmatter    (:,:)  ! net source of dead organic C from living biomass (aggregation, mortality, excretion, sloppy feeding, etc.)
   
!  Domain-local arrays required for REcoM fluxes
   real, allocatable :: s_bio_dicLoc          (:)
   real, allocatable :: s_bio_livingmatterLoc (:)
   real, allocatable :: s_bio_deadmatterLoc   (:)
   real, allocatable :: s_bio_alkLoc          (:)
   
!  ********************************************
!  *** Diagnostic variables (monthly means) ***
!  ********************************************
!  dimension: (nlay, myDim_nod2D)

   real, allocatable :: tM_u_alk             (:,:)
   real, allocatable :: tM_v_alk             (:,:)
   real, allocatable :: tM_w_alk             (:,:)
   real, allocatable :: tM_u_dic             (:,:)
   real, allocatable :: tM_v_dic             (:,:)
   real, allocatable :: tM_w_dic             (:,:)
   real, allocatable :: tM_u_livingmatter    (:,:)
   real, allocatable :: tM_v_livingmatter    (:,:)
   real, allocatable :: tM_w_livingmatter    (:,:)
   real, allocatable :: tM_sink_livingmatter (:,:)
   real, allocatable :: tM_u_deadmatter      (:,:)
   real, allocatable :: tM_v_deadmatter      (:,:)
   real, allocatable :: tM_w_deadmatter      (:,:)
   real, allocatable :: tM_export            (:,:)
   real, allocatable :: sM_hor_alk           (:,:)
   real, allocatable :: sM_ver_alk           (:,:)
   real, allocatable :: sM_diffV_alk         (:,:)
   real, allocatable :: sM_diffH_alk         (:,:)
   real, allocatable :: sM_bio_alk           (:,:)
   real, allocatable :: sM_hor_dic           (:,:)
   real, allocatable :: sM_ver_dic           (:,:)
   real, allocatable :: sM_diffV_dic         (:,:)
   real, allocatable :: sM_diffH_dic         (:,:)
   real, allocatable :: sM_bio_dic           (:,:)
   real, allocatable :: sM_hor_livingmatter  (:,:)
   real, allocatable :: sM_ver_livingmatter  (:,:)
   real, allocatable :: sM_diffV_livingmatter(:,:)
   real, allocatable :: sM_diffH_livingmatter(:,:)
   real, allocatable :: sM_sink_livingmatter (:,:)
   real, allocatable :: sM_bio_livingmatter  (:,:)
   real, allocatable :: sM_hor_deadmatter    (:,:)
   real, allocatable :: sM_ver_deadmatter    (:,:)
   real, allocatable :: sM_diffV_deadmatter  (:,:)
   real, allocatable :: sM_diffH_deadmatter  (:,:)
   real, allocatable :: sM_export            (:,:)
   real, allocatable :: sM_bio_deadmatter    (:,:)
   
CONTAINS

! ***********************************************
! ***                                         ***
! ***   init_carbonfluxes_diags_arrays        ***
! ***                                         ***
! ***********************************************
SUBROUTINE init_carbonfluxes_diags_arrays()
USE g_PARSUP, ONLY: myDim_nod2D, eDim_nod2D
USE mod_assim_pdaf, ONLY: mesh_fesom

implicit none
integer :: nl

   ! To write debugging output:
   cfdiags_debug = .False.
   
   ! True:   Write out carbon concentration fields
   ! False:  Write out carbon mass fields
   cfconc = .False.

   nl = mesh_fesom%nl

   allocate(t_u_alk              (nl-2,myDim_nod2D))
   allocate(t_v_alk              (nl-2,myDim_nod2D))
   allocate(t_w_alk              (nl-2,myDim_nod2D))
   allocate(t_u_dic              (nl-2,myDim_nod2D))
   allocate(t_v_dic              (nl-2,myDim_nod2D))
   allocate(t_w_dic              (nl-2,myDim_nod2D))
   allocate(t_u_livingmatter     (nl-2,myDim_nod2D))
   allocate(t_v_livingmatter     (nl-2,myDim_nod2D))
   allocate(t_w_livingmatter     (nl-2,myDim_nod2D))
   allocate(t_sink_livingmatter  (nl-2,myDim_nod2D))
   allocate(t_u_deadmatter       (nl-2,myDim_nod2D))
   allocate(t_v_deadmatter       (nl-2,myDim_nod2D))
   allocate(t_w_deadmatter       (nl-2,myDim_nod2D))
   allocate(t_export             (nl-2,myDim_nod2D))
   
   allocate(s_hor_alk            (nl-2,myDim_nod2D))
   allocate(s_ver_alk            (nl-2,myDim_nod2D))
   allocate(s_diffV_alk          (nl-2,myDim_nod2D))
   allocate(s_diffH_alk          (nl-2,myDim_nod2D+eDim_nod2D)) ! edges required for horizontal diffusion
   allocate(s_bio_alk            (nl-2,myDim_nod2D))
   allocate(s_hor_dic            (nl-2,myDim_nod2D))
   allocate(s_ver_dic            (nl-2,myDim_nod2D))
   allocate(s_diffV_dic          (nl-2,myDim_nod2D))
   allocate(s_diffH_dic          (nl-2,myDim_nod2D+eDim_nod2D))
   allocate(s_bio_dic            (nl-2,myDim_nod2D))
   allocate(s_hor_livingmatter   (nl-2,myDim_nod2D))
   allocate(s_ver_livingmatter   (nl-2,myDim_nod2D))
   allocate(s_diffV_livingmatter (nl-2,myDim_nod2D))
   allocate(s_diffH_livingmatter (nl-2,myDim_nod2D+eDim_nod2D))
   allocate(s_sink_livingmatter  (nl-2,myDim_nod2D))
   allocate(s_bio_livingmatter   (nl-2,myDim_nod2D))
   allocate(s_hor_deadmatter     (nl-2,myDim_nod2D))
   allocate(s_ver_deadmatter     (nl-2,myDim_nod2D))
   allocate(s_diffV_deadmatter   (nl-2,myDim_nod2D))
   allocate(s_diffH_deadmatter   (nl-2,myDim_nod2D+eDim_nod2D))
   allocate(s_export             (nl-2,myDim_nod2D))
   allocate(s_bio_deadmatter     (nl-2,myDim_nod2D))

   allocate(tM_u_alk             (nl-2,myDim_nod2D))
   allocate(tM_v_alk             (nl-2,myDim_nod2D))
   allocate(tM_w_alk             (nl-2,myDim_nod2D))
   allocate(tM_u_dic             (nl-2,myDim_nod2D))
   allocate(tM_v_dic             (nl-2,myDim_nod2D))
   allocate(tM_w_dic             (nl-2,myDim_nod2D))
   allocate(tM_u_livingmatter    (nl-2,myDim_nod2D))
   allocate(tM_v_livingmatter    (nl-2,myDim_nod2D))
   allocate(tM_w_livingmatter    (nl-2,myDim_nod2D))
   allocate(tM_sink_livingmatter (nl-2,myDim_nod2D))
   allocate(tM_u_deadmatter      (nl-2,myDim_nod2D))
   allocate(tM_v_deadmatter      (nl-2,myDim_nod2D))
   allocate(tM_w_deadmatter      (nl-2,myDim_nod2D))
   allocate(tM_export            (nl-2,myDim_nod2D))
   allocate(sM_hor_alk           (nl-2,myDim_nod2D))
   allocate(sM_ver_alk           (nl-2,myDim_nod2D))
   allocate(sM_diffV_alk         (nl-2,myDim_nod2D))
   allocate(sM_diffH_alk         (nl-2,myDim_nod2D))
   allocate(sM_bio_alk           (nl-2,myDim_nod2D))
   allocate(sM_hor_dic           (nl-2,myDim_nod2D))
   allocate(sM_ver_dic           (nl-2,myDim_nod2D))
   allocate(sM_diffV_dic         (nl-2,myDim_nod2D))
   allocate(sM_diffH_dic         (nl-2,myDim_nod2D))
   allocate(sM_bio_dic           (nl-2,myDim_nod2D))
   allocate(sM_hor_livingmatter  (nl-2,myDim_nod2D))
   allocate(sM_ver_livingmatter  (nl-2,myDim_nod2D))
   allocate(sM_diffV_livingmatter(nl-2,myDim_nod2D))
   allocate(sM_diffH_livingmatter(nl-2,myDim_nod2D))
   allocate(sM_sink_livingmatter (nl-2,myDim_nod2D))
   allocate(sM_bio_livingmatter  (nl-2,myDim_nod2D))
   allocate(sM_hor_deadmatter    (nl-2,myDim_nod2D))
   allocate(sM_ver_deadmatter    (nl-2,myDim_nod2D))
   allocate(sM_diffV_deadmatter  (nl-2,myDim_nod2D))
   allocate(sM_diffH_deadmatter  (nl-2,myDim_nod2D))
   allocate(sM_export            (nl-2,myDim_nod2D))
   allocate(sM_bio_deadmatter    (nl-2,myDim_nod2D))
   
   t_u_alk              =0.0
   t_v_alk              =0.0
   t_w_alk              =0.0
   t_u_dic              =0.0
   t_v_dic              =0.0
   t_w_dic              =0.0
   t_u_livingmatter     =0.0
   t_v_livingmatter     =0.0
   t_w_livingmatter     =0.0
   t_sink_livingmatter  =0.0
   t_u_deadmatter       =0.0
   t_v_deadmatter       =0.0
   t_w_deadmatter       =0.0
   t_export             =0.0
   
   s_hor_alk            =0.0
   s_ver_alk            =0.0
   s_diffV_alk          =0.0
   s_diffH_alk          =0.0
   s_bio_alk            =0.0
   s_hor_dic            =0.0
   s_ver_dic            =0.0
   s_diffV_dic          =0.0
   s_diffH_dic          =0.0
   s_bio_dic            =0.0
   s_hor_livingmatter   =0.0
   s_ver_livingmatter   =0.0
   s_diffV_livingmatter =0.0
   s_diffH_livingmatter =0.0
   s_sink_livingmatter  =0.0
   s_bio_livingmatter   =0.0
   s_hor_deadmatter     =0.0
   s_ver_deadmatter     =0.0
   s_diffV_deadmatter   =0.0
   s_diffH_deadmatter   =0.0
   s_export             =0.0
   s_bio_deadmatter     =0.0

   tM_u_alk             =0.0
   tM_v_alk             =0.0
   tM_w_alk             =0.0
   tM_u_dic             =0.0
   tM_v_dic             =0.0
   tM_w_dic             =0.0
   tM_u_livingmatter    =0.0
   tM_v_livingmatter    =0.0
   tM_w_livingmatter    =0.0
   tM_sink_livingmatter =0.0
   tM_u_deadmatter      =0.0
   tM_v_deadmatter      =0.0
   tM_w_deadmatter      =0.0
   
   tM_export            =0.0
   sM_hor_alk           =0.0
   sM_ver_alk           =0.0
   sM_diffV_alk         =0.0
   sM_diffH_alk         =0.0
   sM_bio_alk           =0.0
   sM_hor_dic           =0.0
   sM_ver_dic           =0.0
   sM_diffV_dic         =0.0
   sM_diffH_dic         =0.0
   sM_bio_dic           =0.0
   sM_hor_livingmatter  =0.0
   sM_ver_livingmatter  =0.0
   sM_diffV_livingmatter=0.0
   sM_diffH_livingmatter=0.0
   sM_sink_livingmatter =0.0
   sM_bio_livingmatter  =0.0
   sM_hor_deadmatter    =0.0
   sM_ver_deadmatter    =0.0
   sM_diffV_deadmatter  =0.0
   sM_diffH_deadmatter  =0.0
   sM_export            =0.0
   sM_bio_deadmatter    =0.0

END SUBROUTINE init_carbonfluxes_diags_arrays

! ***********************************************
! ***                                         ***
! ***   carbonfluxes_diags_output_monthly     ***
! ***                                         ***
! ***********************************************
! Compute monthly means ensemble means and writes output.

SUBROUTINE carbonfluxes_diags_output_monthly(mstep)

      USE g_clock, &
        ONLY: month,num_day_in_month,fleapyear,cyearnew
      USE g_config, &
        ONLY: step_per_day
      USE g_parsup, &
        ONLY: myDim_nod2D
      USE mod_assim_pdaf, &
        ONLY: nlmax, dim_ens
      USE mod_parallel_pdaf, &
        ONLY: mype_world, abort_parallel, task_id, mype_submodel, &
        COMM_COUPLE, filterpe, writepe


      IMPLICIT NONE
      include 'mpif.h'
      
      ! arguments
      INTEGER, INTENT(in) :: mstep
      ! local variables
      LOGICAL :: now_to_write_monthly
      INTEGER :: month_iter, whichmonth,fleap
      REAL    :: weights
      INTEGER :: mpierror
      character(len=200) :: filename
      
      now_to_write_monthly = .FALSE.

      ! *** set "now_to_write_monthly" at last day of month ***
      call monthly_event(now_to_write_monthly)
      
      IF (writepe) WRITE (*, '(/a, 1x, a, 1x, i3, 1x, a, 1x, l)') 'FESOM-PDAF', 'Call carbonflux diagnostics at step', mstep, 'nowtowritemonthly', now_to_write_monthly
      
      IF (.not. now_to_write_monthly) THEN
      ! ***
      ! *** add instantaneous data to monthly
      ! ***
      
      tM_u_alk              = tM_u_alk              + t_u_alk            
      tM_v_alk              = tM_v_alk              + t_v_alk            
      tM_w_alk              = tM_w_alk              + t_w_alk
      tM_u_dic              = tM_u_dic              + t_u_dic
      tM_v_dic              = tM_v_dic              + t_v_dic
      tM_w_dic              = tM_w_dic              + t_w_dic
      tM_u_livingmatter     = tM_u_livingmatter     + t_u_livingmatter
      tM_v_livingmatter     = tM_v_livingmatter     + t_v_livingmatter
      tM_w_livingmatter     = tM_w_livingmatter     + t_w_livingmatter
      tM_sink_livingmatter  = tM_sink_livingmatter  + t_sink_livingmatter
      tM_u_deadmatter       = tM_u_deadmatter       + t_u_deadmatter
      tM_v_deadmatter       = tM_v_deadmatter       + t_v_deadmatter
      tM_w_deadmatter       = tM_w_deadmatter       + t_w_deadmatter
      tM_export             = tM_export             + t_export
      
      sM_hor_alk            = sM_hor_alk            + s_hor_alk
      sM_ver_alk            = sM_ver_alk            + s_ver_alk
      sM_diffV_alk          = sM_diffV_alk          + s_diffV_alk
      sM_diffH_alk          = sM_diffH_alk          + s_diffH_alk          (:,:myDim_nod2D)
      sM_bio_alk            = sM_bio_alk            + s_bio_alk
      sM_hor_dic            = sM_hor_dic            + s_hor_dic
      sM_ver_dic            = sM_ver_dic            + s_ver_dic
      sM_diffV_dic          = sM_diffV_dic          + s_diffV_dic
      sM_diffH_dic          = sM_diffH_dic          + s_diffH_dic          (:,:myDim_nod2D)
      sM_bio_dic            = sM_bio_dic            + s_bio_dic
      sM_hor_livingmatter   = sM_hor_livingmatter   + s_hor_livingmatter
      sM_ver_livingmatter   = sM_ver_livingmatter   + s_ver_livingmatter
      sM_diffV_livingmatter = sM_diffV_livingmatter + s_diffV_livingmatter
      sM_diffH_livingmatter = sM_diffH_livingmatter + s_diffH_livingmatter (:,:myDim_nod2D)
      sM_sink_livingmatter  = sM_sink_livingmatter  + s_sink_livingmatter
      sM_bio_livingmatter   = sM_bio_livingmatter   + s_bio_livingmatter
      sM_hor_deadmatter     = sM_hor_deadmatter     + s_hor_deadmatter
      sM_ver_deadmatter     = sM_ver_deadmatter     + s_ver_deadmatter
      sM_diffV_deadmatter   = sM_diffV_deadmatter   + s_diffV_deadmatter
      sM_diffH_deadmatter   = sM_diffH_deadmatter   + s_diffH_deadmatter   (:,:myDim_nod2D)
      sM_export             = sM_export             + s_export
      sM_bio_deadmatter     = sM_bio_deadmatter     + s_bio_deadmatter
      
      ELSE
      ! ***
      ! *** compute monthly mean, compute ensemble mean, write output and reset monthly data to zero
      ! ***
      
      ! apply weighting to monthly mean
      weights = 1.0/REAL(num_day_in_month(fleapyear,month)*step_per_day)/REAL(dim_ens)
      
      tM_u_alk             = tM_u_alk              * weights
      tM_v_alk             = tM_v_alk              * weights
      tM_w_alk             = tM_w_alk              * weights
      tM_u_dic             = tM_u_dic              * weights
      tM_v_dic             = tM_v_dic              * weights
      tM_w_dic             = tM_w_dic              * weights
      tM_u_livingmatter    = tM_u_livingmatter     * weights
      tM_v_livingmatter    = tM_v_livingmatter     * weights
      tM_w_livingmatter    = tM_w_livingmatter     * weights
      tM_sink_livingmatter = tM_sink_livingmatter  * weights
      tM_u_deadmatter      = tM_u_deadmatter       * weights
      tM_v_deadmatter      = tM_v_deadmatter       * weights
      tM_w_deadmatter      = tM_w_deadmatter       * weights
      tM_export            = tM_export             * weights
      
      sM_hor_alk           = sM_hor_alk            * weights
      sM_ver_alk           = sM_ver_alk            * weights
      sM_diffV_alk         = sM_diffV_alk          * weights
      sM_diffH_alk         = sM_diffH_alk          * weights
      sM_bio_alk           = sM_bio_alk            * weights
      sM_hor_dic           = sM_hor_dic            * weights
      sM_ver_dic           = sM_ver_dic            * weights
      sM_diffV_dic         = sM_diffV_dic          * weights
      sM_diffH_dic         = sM_diffH_dic          * weights
      sM_bio_dic           = sM_bio_dic            * weights
      sM_hor_livingmatter  = sM_hor_livingmatter   * weights
      sM_ver_livingmatter  = sM_ver_livingmatter   * weights
      sM_diffV_livingmatter= sM_diffV_livingmatter * weights
      sM_diffH_livingmatter= sM_diffH_livingmatter * weights
      sM_sink_livingmatter = sM_sink_livingmatter  * weights
      sM_bio_livingmatter  = sM_bio_livingmatter   * weights
      sM_hor_deadmatter    = sM_hor_deadmatter     * weights
      sM_ver_deadmatter    = sM_ver_deadmatter     * weights
      sM_diffV_deadmatter  = sM_diffV_deadmatter   * weights
      sM_diffH_deadmatter  = sM_diffH_deadmatter   * weights
      sM_export            = sM_export             * weights
      sM_bio_deadmatter    = sM_bio_deadmatter     * weights

      
      ! compute ensemble mean
      IF (filterpe) THEN
         !                - send -       - receive -            - size -             - type -              - sum -  -receiver-   -ensemble-    -check-
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_u_alk             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_v_alk             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_w_alk             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_u_dic             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_v_dic             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_w_dic             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_u_livingmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_v_livingmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_w_livingmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_sink_livingmatter , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_u_deadmatter      , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_v_deadmatter      , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_w_deadmatter      , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  tM_export            , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_hor_alk           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_ver_alk           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffV_alk         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffH_alk         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_bio_alk           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_hor_dic           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_ver_dic           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffV_dic         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffH_dic         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_bio_dic           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_hor_livingmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_ver_livingmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffV_livingmatter, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffH_livingmatter, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_sink_livingmatter , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_bio_livingmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_hor_deadmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_ver_deadmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffV_deadmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_diffH_deadmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_export            , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( MPI_IN_PLACE,  sM_bio_deadmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM,  0          , COMM_COUPLE, mpierror)
         
         ! write output
         call write_carbonfluxes_diags_out(month)
      ELSE
         
         !                - send -               - receive -            - size -             - type -              - sum -  -receiver-   -ensemble-    -check-
         CALL MPI_REDUCE( tM_u_alk             , tM_u_alk             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_v_alk             , tM_v_alk             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_w_alk             , tM_w_alk             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_u_dic             , tM_u_dic             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_v_dic             , tM_v_dic             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_w_dic             , tM_w_dic             , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_u_livingmatter    , tM_u_livingmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_v_livingmatter    , tM_v_livingmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_w_livingmatter    , tM_w_livingmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_sink_livingmatter , tM_sink_livingmatter , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_u_deadmatter      , tM_u_deadmatter      , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_v_deadmatter      , tM_v_deadmatter      , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_w_deadmatter      , tM_w_deadmatter      , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( tM_export            , tM_export            , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_hor_alk           , sM_hor_alk           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_ver_alk           , sM_ver_alk           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffV_alk         , sM_diffV_alk         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffH_alk         , sM_diffH_alk         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_bio_alk           , sM_bio_alk           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_hor_dic           , sM_hor_dic           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_ver_dic           , sM_ver_dic           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffV_dic         , sM_diffV_dic         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffH_dic         , sM_diffH_dic         , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_bio_dic           , sM_bio_dic           , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_hor_livingmatter  , sM_hor_livingmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_ver_livingmatter  , sM_ver_livingmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffV_livingmatter, sM_diffV_livingmatter, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffH_livingmatter, sM_diffH_livingmatter, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_sink_livingmatter , sM_sink_livingmatter , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_bio_livingmatter  , sM_bio_livingmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_hor_deadmatter    , sM_hor_deadmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_ver_deadmatter    , sM_ver_deadmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffV_deadmatter  , sM_diffV_deadmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_diffH_deadmatter  , sM_diffH_deadmatter  , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_export            , sM_export            , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( sM_bio_deadmatter    , sM_bio_deadmatter    , nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         
      ENDIF ! filterpe
      
      ! reset monthly data to zero
      tM_u_alk             = 0.0
      tM_v_alk             = 0.0
      tM_w_alk             = 0.0
      tM_u_dic             = 0.0
      tM_v_dic             = 0.0
      tM_w_dic             = 0.0
      tM_u_livingmatter    = 0.0
      tM_v_livingmatter    = 0.0
      tM_w_livingmatter    = 0.0
      tM_sink_livingmatter = 0.0
      tM_u_deadmatter      = 0.0
      tM_v_deadmatter      = 0.0
      tM_w_deadmatter      = 0.0
      tM_export            = 0.0
      sM_hor_alk           = 0.0
      sM_ver_alk           = 0.0
      sM_diffV_alk         = 0.0
      sM_diffH_alk         = 0.0
      sM_bio_alk           = 0.0
      sM_hor_dic           = 0.0
      sM_ver_dic           = 0.0
      sM_diffV_dic         = 0.0
      sM_diffH_dic         = 0.0
      sM_bio_dic           = 0.0
      sM_hor_livingmatter  = 0.0
      sM_ver_livingmatter  = 0.0
      sM_diffV_livingmatter= 0.0
      sM_diffH_livingmatter= 0.0
      sM_sink_livingmatter = 0.0
      sM_bio_livingmatter  = 0.0
      sM_hor_deadmatter    = 0.0
      sM_ver_deadmatter    = 0.0
      sM_diffV_deadmatter  = 0.0
      sM_diffH_deadmatter  = 0.0
      sM_export            = 0.0
      sM_bio_deadmatter    = 0.0
      
      ENDIF ! now_to_write_monthly
  
END SUBROUTINE carbonfluxes_diags_output_monthly

! ********************************
! ***                          ***
! ***   netCDF check           ***
! ***                          ***
! ********************************
! Checks for errors during netCDF operations.

SUBROUTINE check(status)
     USE netcdf
     USE mod_parallel_pdaf, ONLY: abort_parallel
     
     IMPLICIT NONE

     ! *** Arguments ***
     integer, intent ( in) :: status   ! Reading status

     if(status /= nf90_noerr) then
        print *, trim(nf90_strerror(status))
        call abort_parallel()
     end if

END SUBROUTINE check

! ***********************************************
! ***                                         ***
! ***   init_carbonfluxes_diags_out           ***
! ***                                         ***
! ***********************************************
! Initializes netCDF output file for carbon flux diagnostics. 
SUBROUTINE init_carbonfluxes_diags_out()
      
      USE mod_assim_pdaf, &
         ONLY: mesh_fesom, nlmax, DAoutput_path
      USE mod_obs_f_pdaf, &
         ONLY: pi
      USE g_config, &
         ONLY: runid
      USE recom_config, &
         ONLY: secondsperday
      USE g_clock, &
         ONLY: cyearnew, num_day_in_month, fleapyear, yearold
      USE g_parsup, &
         ONLY: myDim_nod2D
      USE g_comm_auto, &
         ONLY: gather_nod
      USE mod_parallel_pdaf, &
         ONLY: writepe
      USE netcdf

      IMPLICIT NONE
      
      character(len=200) :: filename
      INTEGER :: fileid
      INTEGER :: dimIDs(3)
      INTEGER :: varid_nod2, varid_iter, &
                 varid_nz1, &
                 varid_lon, varid_lat
      INTEGER :: dimid_nod2, dimid_iter, &
                 dimid_nz1
      INTEGER :: ndims
      INTEGER :: varid
      INTEGER :: n
      character(40) :: att_text
      INTEGER :: firstdayofmonth(12)
      
      REAL, allocatable :: lon(:)
      REAL, allocatable :: lat(:)
      
      ! gather GEO coordinates (from all PEs)
      allocate(lon(mesh_fesom% nod2D),lat(mesh_fesom% nod2D))
      call gather_nod(mesh_fesom%geo_coord_nod2D(1, 1:myDim_nod2D), lon)
      call gather_nod(mesh_fesom%geo_coord_nod2D(2, 1:myDim_nod2D), lat)
      
      IF (writepe) THEN
      
      ! initialize file
      filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.cfx.'//cyearnew//'.nc'
      WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize carbon flux diags NetCDF file: '//TRIM(filename)
      
      ! open file
      call check(NF90_CREATE(trim(filename),NF90_NETCDF4,fileid))
      
      ! define dimensions
      call check( NF90_DEF_DIM(fileid, 'nod2', mesh_fesom%nod2d, dimID_nod2))
      call check( NF90_DEF_DIM(fileid, 'nz1',  nlmax,            dimID_nz1))
      call check( NF90_DEF_DIM(fileid, 'time', 12,               dimID_iter))
      
      ! dimension variables
      call check( nf90_def_var(fileid, 'nod2', NF90_INT,   dimID_nod2, varid_nod2))
      call check( nf90_def_var(fileid, 'nz1',  NF90_FLOAT, dimID_nz1,  varid_nz1))
      call check( nf90_def_var(fileid, 'time', NF90_INT,   dimID_iter, varid_iter))
      
      call check( nf90_def_var(fileid, 'lon', NF90_FLOAT,   dimID_nod2, varid_lon))
      call check( nf90_def_var(fileid, 'lat', NF90_FLOAT,   dimID_nod2, varid_lat))
      
      ! dimension description
      call check( nf90_put_att(fileid, varid_nod2, 'long_name', 'surface nodes'))
      call check( nf90_put_att(fileid, varid_nz1,  'long_name', 'vertical layers (mid-layer depths)'))
      call check( nf90_put_att(fileid, varid_nz1,  'units', 'm'))
      call check( nf90_put_att(fileid, varid_iter, 'long_name', 'time month'))
      
      call check( nf90_put_att(fileid, varid_lon,  'long_name', 'longitude'))
      call check( nf90_put_att(fileid, varid_lat,  'long_name', 'latitude'))
      call check( nf90_put_att(fileid, varid_lon,  'units', 'degE [-180;180]'))
      call check( nf90_put_att(fileid, varid_lat,  'units', 'degN [-90;90]'))
      
      write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
      call check( nf90_put_att(fileid, varid_iter, 'long_name', 'time'))
      call check( nf90_put_att(fileid, varid_iter, 'standard_name', 'time'))
      call check( nf90_put_att(fileid, varid_iter, 'units', trim(att_text)))
      call check( nf90_put_att(fileid, varid_iter, 'axis', 'T'))
      call check( nf90_put_att(fileid, varid_iter, 'stored_direction', 'increasing'))
      
      ! fill dimension variables
      call check (nf90_enddef(fileid))
      call check (nf90_put_var(fileid, varid_nz1,  mesh_fesom% Z   (1:nlmax)))
      call check (nf90_put_var(fileid, varid_nod2, [(n,n=1,mesh_fesom% nod2D)]))
      
      call check (nf90_put_var(fileid, varid_lon,  REAL(180./pi * lon, 4)))
      call check (nf90_put_var(fileid, varid_lat,  REAL(180./pi * lat, 4)))
      
      firstdayofmonth(1)=0
      DO n=2,12
         firstdayofmonth(n) = firstdayofmonth(n-1)+num_day_in_month(fleapyear,n-1)*SecondsPerDay
      ENDDO
      call check( nf90_put_var(fileid, varid_iter, firstdayofmonth))
      
      ! define diagnostic variables
      call check (nf90_redef(fileid))

      ! set dimensions
      ndims = 3
      dimIDs(1) = dimID_nod2
      dimIDs(2) = dimID_nz1
      dimIDs(3) = dimID_iter
      
      call check( NF90_DEF_VAR(fileid, 't_u_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_v_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_w_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_u_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_v_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_w_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_u_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_v_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_w_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_sink_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_u_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_v_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_w_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol s$^{-1}$ m$^{-2}$'))
      
      call check( NF90_DEF_VAR(fileid, 't_export', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_hor_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_ver_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffV_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))

      call check( NF90_DEF_VAR(fileid, 's_diffH_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_bio_alk', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_hor_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_ver_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffV_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffH_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_bio_dic', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_hor_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_ver_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffV_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffH_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_sink_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_bio_livingmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_hor_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_ver_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffV_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_diffH_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_export', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check( NF90_DEF_VAR(fileid, 's_bio_deadmatter', NF90_FLOAT, dimIDs(1:ndims), varid))
      call check( nf90_put_att(fileid, varid, 'long_name', ''))
      call check( nf90_put_att(fileid, varid, 'units',     'mmol  m$^{-3}$ s$^{-1}$'))
      
      call check(NF90_ENDDEF(fileid))
      call check(NF90_CLOSE(fileid))
      
      ENDIF

      deallocate(lat,lon)

END SUBROUTINE init_carbonfluxes_diags_out

! ***********************************************
! ***                                         ***
! ***   write_carbonfluxes_diags_out          ***
! ***                                         ***
! ***********************************************
! Writes carbon flux diagnostics to netCDF. 
SUBROUTINE write_carbonfluxes_diags_out(writepos)

USE g_clock, &
   ONLY: month, cyearnew
USE mod_assim_pdaf, &
   ONLY: nlmax, mesh_fesom, DAoutput_path
USE g_parsup, &
   ONLY: myDim_nod2D
USE mod_parallel_pdaf, &
   ONLY: writepe
USE g_comm_auto, &
   ONLY: gather_nod
USE netcdf

IMPLICIT NONE

! ARGUMENTS:
INTEGER, INTENT(in) :: writepos                    ! Write position

! LOCAL VARIABLES:
REAL, allocatable  :: data3_g(:,:)                 ! Temporary array for global 3D-fields
character(len=200) :: filename                     ! Full name of output file
integer            :: fileid                       ! nc-file ID for output file
character(len=100) :: varname

filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.cfx.'//cyearnew//'.nc'
! print screen information:
IF (writepe) THEN
   WRITE (*, '(a, 8x, a, i9, a, a)') 'FESOM-PDAF', 'Write carbon flux diagnostics to NetCDF at for month ', &
   month, ' to file ', TRIM(filename)
   ! open file
   call check( nf90_open(TRIM(filename), nf90_write, fileid))
ENDIF ! writepe

! gather and write global ocean fields
allocate(data3_g(nlmax,mesh_fesom% nod2D))

varname =      't_u_alk'            
CALL gather_nod(tM_u_alk            , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_v_alk'            
CALL gather_nod(tM_v_alk            , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_w_alk'            
CALL gather_nod(tM_w_alk            , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_u_dic'            
CALL gather_nod(tM_u_dic            , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_v_dic'            
CALL gather_nod(tM_v_dic            , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_w_dic'            
CALL gather_nod(tM_w_dic            , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_u_livingmatter'   
CALL gather_nod(tM_u_livingmatter   , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_v_livingmatter'   
CALL gather_nod(tM_v_livingmatter   , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_w_livingmatter'   
CALL gather_nod(tM_w_livingmatter   , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_sink_livingmatter'
CALL gather_nod(tM_sink_livingmatter, data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_u_deadmatter'     
CALL gather_nod(tM_u_deadmatter     , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_v_deadmatter'     
CALL gather_nod(tM_v_deadmatter     , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_w_deadmatter'     
CALL gather_nod(tM_w_deadmatter     , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      't_export'           
CALL gather_nod(tM_export           , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_hor_alk'          
CALL gather_nod(sM_hor_alk          , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_ver_alk'          
CALL gather_nod(sM_ver_alk          , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffV_alk'         
CALL gather_nod(sM_diffV_alk         , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffH_alk'         
CALL gather_nod(sM_diffH_alk         , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_bio_alk'         
CALL gather_nod(sM_bio_alk         , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_hor_dic'          
CALL gather_nod(sM_hor_dic          , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_ver_dic'          
CALL gather_nod(sM_ver_dic          , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffV_dic'         
CALL gather_nod(sM_diffV_dic         , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffH_dic'         
CALL gather_nod(sM_diffH_dic         , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_bio_dic'          
CALL gather_nod(sM_bio_dic          , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_hor_livingmatter' 
CALL gather_nod(sM_hor_livingmatter , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_ver_livingmatter' 
CALL gather_nod(sM_ver_livingmatter , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffV_livingmatter'
CALL gather_nod(sM_diffV_livingmatter, data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffH_livingmatter'
CALL gather_nod(sM_diffH_livingmatter, data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_sink_livingmatter'
CALL gather_nod(sM_sink_livingmatter, data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_bio_livingmatter' 
CALL gather_nod(sM_bio_livingmatter , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_hor_deadmatter'   
CALL gather_nod(sM_hor_deadmatter   , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_ver_deadmatter'   
CALL gather_nod(sM_ver_deadmatter   , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffV_deadmatter'  
CALL gather_nod(sM_diffV_deadmatter  , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_diffH_deadmatter'  
CALL gather_nod(sM_diffH_deadmatter  , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_export'           
CALL gather_nod(sM_export           , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

varname =      's_bio_deadmatter'   
CALL gather_nod(sM_bio_deadmatter   , data3_g) 
IF (writepe) call putvar(fileid,varname,data3_g,writepos)

deallocate(data3_g)

! close file:
IF (writepe) call check (nf90_close(fileid))
END SUBROUTINE write_carbonfluxes_diags_out

! ***********************************************
! ***                                         ***
! ***   putvar                                ***
! ***                                         ***
! ***********************************************
! Puts variable to netCDF. 
SUBROUTINE putvar(fileid,varname,data3_g,writepos)
   USE g_clock, &
      ONLY: month
   USE mod_assim_pdaf, &
      ONLY: nlmax, mesh_fesom
   USE netcdf

   implicit none
   ! arguments:
   integer, intent(in) :: fileid
   character(len=100), intent(in) :: varname
   real, intent(in) :: data3_g(nlmax,mesh_fesom%nod2D)
   integer, intent(in) :: writepos
   ! local variables:
   integer :: varid
   
   ! inquire variable ID
   call check( nf90_inq_varid(fileid, TRIM(varname), varid))
   ! put data
   call check( nf90_put_var(fileid, varid, REAL(TRANSPOSE(data3_g),4), &
                                   start=(/ 1, 1, writepos /), &
                                   count=(/ mesh_fesom% nod2D, nlmax, 1 /) ))
                                   ! dims: 1-nod2, 2-nz / nz1, 3-iter
END SUBROUTINE putvar

! ***********************************************
! ***                                         ***
! ***   debug                                 ***
! ***                                         ***
! ***********************************************
! check if fluxes add up to conserve mass

! vertical fluxes

SUBROUTINE debug_vert(varname,vardata)
   USE g_parsup, &
      ONLY: myDim_nod2D
   USE mod_parallel_pdaf, &
      ONLY: writepe
   USE g_comm_auto, &
      ONLY: gather_nod 
   USE mod_assim_pdaf, &
      ONLY: istep_asml,nlmax, mesh_fesom, DAoutput_path
      
   implicit none
   
   character(len=100), intent(in) :: varname
   real, intent(in)               :: vardata(nlmax,myDim_nod2D)
   integer :: n, fid
   
   if (writepe) then
   
     open(fid,file=TRIM(DAoutput_path)//TRIM(varname)//'.out', status='replace', action='write')
     
     DO n=1,myDim_nod2D
       write(fid, '(i6,1x,F10.2)') n, sum(vardata(:,n))
     ENDDO
     
     close(fid)
   
   endif

END SUBROUTINE debug_vert

! horizontal fluxes

SUBROUTINE debug_hor(varname,vardata)
   USE g_parsup, &
      ONLY: myDim_nod2D
   USE mod_parallel_pdaf, &
      ONLY: writepe
   USE g_comm_auto, &
      ONLY: gather_nod 
   USE mod_assim_pdaf, &
      ONLY: istep_asml,nlmax, mesh_fesom, DAoutput_path
      
   implicit none
   
   character(len=100), intent(in) :: varname
   real, intent(in)               :: vardata(nlmax,myDim_nod2D)
   integer :: n, fid
   real :: vardata_glob(nlmax,mesh_fesom%nod2D)
   
   CALL gather_nod(vardata,vardata_glob)
   
   if (writepe) then
   
     open(fid,file=TRIM(DAoutput_path)//TRIM(varname)//'.out', status='replace', action='write')
     
     DO n=1,nlmax
       write(fid, '(i6,1x,F10.2)') n, sum(vardata_glob(n,:))
     ENDDO
     
     close(fid)
   
   endif ! writepe

END SUBROUTINE debug_hor

END MODULE mod_carbon_fluxes_diags





