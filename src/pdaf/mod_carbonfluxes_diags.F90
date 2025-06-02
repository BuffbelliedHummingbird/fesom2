MODULE mod_carbon_fluxes_diags

   use mod_assim_pdaf, only: nlmax

   IMPLICIT NONE

   save
 
   logical            :: cfdiags_debug ! whether to write debug output for tracer mass conservation
   character(len=100) :: vname_cfdiags ! variable name for debugging output
   
!  _________________________________________________________________________   
!  *** Transports, i.e. directed fluxes through area [mmol C / m2 / sec] ***

!  alkalinity
   integer :: id_t_u_alk             = 1  ! advective transport
   integer :: id_t_v_alk             = 2
   integer :: id_t_w_alk             = 3

!  DIC
   integer :: id_t_u_dic             = 4  ! advective transport
   integer :: id_t_v_dic             = 5
   integer :: id_t_w_dic             = 6
   
!  Living biomass (Phy, Dia, Het, Zoo)
   integer :: id_t_u_livingmatter    = 7  ! advective transport
   integer :: id_t_v_livingmatter    = 8
   integer :: id_t_w_livingmatter    = 9
   integer :: id_t_sink_livingmatter = 10 ! sinking in the water column

!  Dead biomass (Detritus and DOC)
   integer :: id_t_u_deadmatter     = 11 ! advective transport of detritus and DOC
   integer :: id_t_v_deadmatter     = 12
   integer :: id_t_w_deadmatter     = 13
   integer :: id_t_sink_deadmatter  = 14 ! sinking of detritus in the water column

!  _____________________________________________________   
!  *** Sources minus sinks: mass          [mmol C / m3 / sec] ***
!  ***                      concentration [mmol C / sec]      ***

!  alkalinity
   integer :: id_s_horLO_alk         = 15  ! advection
   integer :: id_s_verLO_alk         = 16
   integer :: id_s_horAD_alk         = 55 
   integer :: id_s_verAD_alk         = 56
   integer :: id_s_diffV_alk         = 17  ! vertical diffusion
   integer :: id_s_diffH_alk         = 18  ! horizontal diffusion
   integer :: id_s_bio_alk           = 19  ! remineralization and calcification
   integer :: id_s_benthos_alk       = 20  ! remineralization and diffusion from benthos
   integer :: id_s_surf_alk          = 49  ! surface flux due to alkalinity restoring
   integer :: id_s_vol_alk           = 51  ! concentration change due to cell volume change

!  DIC
   integer :: id_s_horLO_dic         = 21  ! advection
   integer :: id_s_verLO_dic         = 22
   integer :: id_s_horAD_dic         = 57
   integer :: id_s_verAD_dic         = 58
   integer :: id_s_diffV_dic         = 23  ! diffusion
   integer :: id_s_diffH_dic         = 24
   integer :: id_s_bio_dic           = 25  ! source of DIC from deadmatter (remineralization of DOC; calcite dissolution)
   integer :: id_s_benthos_dic       = 26  ! remineralization and diffusion from benthos
   integer :: id_s_surf_dic          = 50  ! surface flux due to air-sea exchange
   integer :: id_s_vol_dic           = 52  ! concentration change due to cell volume change


!  Living biomass (Phy, Dia, Het, Zoo)
   integer :: id_s_horLO_livingmatter   = 27  ! advection
   integer :: id_s_verLO_livingmatter   = 28
   integer :: id_s_horAD_livingmatter   = 59
   integer :: id_s_verAD_livingmatter   = 60
   integer :: id_s_diffV_livingmatter   = 29  ! diffusion
   integer :: id_s_diffH_livingmatter   = 30
   integer :: id_s_sink_livingmatter    = 31  ! sinking in the water column
   integer :: id_s_bio_livingmatter     = 32  ! net source of alive biomass from DIC (photosynthesis - respiration; calcification)
   integer :: id_s_benthos_livingmatter = 33  ! sinking into benthos
   integer :: id_s_vol_livingmatter     = 53  ! concentration change due to cell volume change


!  Dead organic carbon (Detritus and DOC)
   integer :: id_s_horLO_deadmatter   = 34  ! advection of detritus and DOC
   integer :: id_s_verLO_deadmatter   = 35
   integer :: id_s_horAD_deadmatter   = 61
   integer :: id_s_verAD_deadmatter   = 62
   integer :: id_s_diffV_deadmatter   = 36  ! diffusion 
   integer :: id_s_diffH_deadmatter   = 37
   integer :: id_s_sink_deadmatter    = 38  ! sinking of detritus in the water column
   integer :: id_s_bio_deadmatter     = 39  ! net source of dead organic C from living biomass (aggregation, mortality, excretion, sloppy feeding, etc.)
   integer :: id_s_benthos_deadmatter = 40  ! sinking into benthos
   integer :: id_s_vol_deadmatter     = 54  ! concentration change due to cell volume change


!  ___________________________________________
!  *** Tracer mass [mmol C ]               ***
!             concentration [mmol C / m3]  ***

   integer :: id_m_alk               = 41
   integer :: id_m_dic               = 42
   integer :: id_m_livingmatter      = 43
   integer :: id_m_deadmatter        = 44
   
!  ____________________________________________________________________________
!  *** Source minus sinks through assimilation mass [mmol C ]               ***
!                                              concentration [mmol C / m3]  ***

   integer :: id_s_asml_alk               = 45
   integer :: id_s_asml_dic               = 46
   integer :: id_s_asml_livingmatter      = 47
   integer :: id_s_asml_deadmatter        = 48
   
   ! Number of fields
   integer :: cfnfields = 62
   
   integer :: cffieldsflux(14)   = [  1, 2, 3, 4, 5, 6, 7, 8, 9,10, &             ! Directed fluxes through area
                                     11,12,13,14 ]
   integer :: cffieldssms(36)    = [             15,16,17,18,19,20, &             ! Sources minus sinks in model
                                     21,22,23,24,25,26,27,28,29,30, &
                                     31,32,33,34,35,36,37,38,39,40, &
                                                             49,50, &
                                                 55,56,57,58,59,60, &
                                     61,62 ]
   integer :: cffieldstracer(4)  = [ 41,42,43,44 ]                                ! Tracers
   integer :: cffieldsasml(4)    = [             45,46,47,48 ]                    ! Sources minus sinks assimilation step
   integer :: cffieldsvol(4)     = [ 51,52,53,54 ]                                ! Concentration change cell volume

   type cffieldtype
     real, allocatable  :: instantmass  (:,:) !  model data is collected into instanteous arrays (mass)
     real, allocatable  :: instantconc  (:,:) !  """                                             (concentration)
     real, allocatable  :: loc          (:)   !  recom-specific terms need to be collected into domain-local arrays
     real, allocatable  :: timemeanmass (:,:) !  time mean (mass)
     real, allocatable  :: timemeanconc (:,:) !  """       (concentration)
     real, allocatable  :: ensmmass     (:,:) !  ensemble mean (mass)
     real, allocatable  :: ensmconc     (:,:) !                (concentration)
     character(len=200) :: unitsmass          !  units mass
     character(len=200) :: unitsconc          !  units concentration
     character(len=100) :: varname            !  name to save output
     real, allocatable  :: fmass        (:,:) !  forecast field (mass)
     real, allocatable  :: fconc        (:,:) !  """            (concentration)
     real, allocatable  :: amass        (:,:) !  analysis field (mass)
     real, allocatable  :: aconc        (:,:) !  """            (concentration)
   end type cffieldtype
   
   type(cffieldtype), allocatable :: cffields(:)
  
   ! output frequency
   character(len=1) :: cfoutfreq_unit = 'm'
   integer          :: cfoutfreq = 1
   
   logical          :: outconc = .false. ! whether to write output for concentration
   logical          :: outmass = .true.  ! """                         mass

   
CONTAINS

! ***********************************************
! ***                                         ***
! ***   init_carbonfluxes_diags_arrays        ***
! ***                                         ***
! ***********************************************
SUBROUTINE init_carbonfluxes_diags_arrays()
USE g_PARSUP, ONLY: myDim_nod2D, eDim_nod2D            ! model grid dimensions
USE mod_assim_pdaf, ONLY: factor_mass, factor_conc     ! used in prepoststep
USE mod_parallel_pdaf, ONLY: writepe

implicit none
integer :: i, ids ! counters


   ! To write debugging output:
   cfdiags_debug = .False.
   
   allocate(cffields(cfnfields))
   
   ! __________________
   ! ___ give names ___
   
   cffields( id_t_u_alk               )%varname = 't_u_alk'
   cffields( id_t_v_alk               )%varname = 't_v_alk'
   cffields( id_t_w_alk               )%varname = 't_w_alk'
   cffields( id_t_u_dic               )%varname = 't_u_dic'
   cffields( id_t_v_dic               )%varname = 't_v_dic'
   cffields( id_t_w_dic               )%varname = 't_w_dic'
   cffields( id_t_u_livingmatter      )%varname = 't_u_livingmatter'
   cffields( id_t_v_livingmatter      )%varname = 't_v_livingmatter'
   cffields( id_t_w_livingmatter      )%varname = 't_w_livingmatter'
   cffields( id_t_sink_livingmatter   )%varname = 't_sink_livingmatter'
   cffields( id_t_u_deadmatter        )%varname = 't_u_deadmatter'
   cffields( id_t_v_deadmatter        )%varname = 't_v_deadmatter'
   cffields( id_t_w_deadmatter        )%varname = 't_w_deadmatter'
   cffields( id_t_sink_deadmatter     )%varname = 't_sink_deadmatter'
   cffields( id_s_horLO_alk           )%varname = 's_horLO_alk'
   cffields( id_s_verLO_alk           )%varname = 's_verLO_alk'
   cffields( id_s_horAD_alk           )%varname = 's_horAD_alk'
   cffields( id_s_verAD_alk           )%varname = 's_verAD_alk'
   cffields( id_s_diffV_alk           )%varname = 's_diffV_alk'
   cffields( id_s_diffH_alk           )%varname = 's_diffH_alk'
   cffields( id_s_bio_alk             )%varname = 's_bio_alk'
   cffields( id_s_horLO_dic           )%varname = 's_horLO_dic'
   cffields( id_s_verLO_dic           )%varname = 's_verLO_dic'
   cffields( id_s_horAD_dic           )%varname = 's_horAD_dic'
   cffields( id_s_verAD_dic           )%varname = 's_verAD_dic'
   cffields( id_s_diffV_dic           )%varname = 's_diffV_dic'
   cffields( id_s_diffH_dic           )%varname = 's_diffH_dic'
   cffields( id_s_bio_dic             )%varname = 's_bio_dic'
   cffields( id_s_horLO_livingmatter  )%varname = 's_horLO_livingmatter'
   cffields( id_s_verLO_livingmatter  )%varname = 's_verLO_livingmatter'
   cffields( id_s_horAD_livingmatter  )%varname = 's_horAD_livingmatter'
   cffields( id_s_verAD_livingmatter  )%varname = 's_verAD_livingmatter'
   cffields( id_s_diffV_livingmatter  )%varname = 's_diffV_livingmatter'
   cffields( id_s_diffH_livingmatter  )%varname = 's_diffH_livingmatter'
   cffields( id_s_sink_livingmatter   )%varname = 's_sink_livingmatter'
   cffields( id_s_bio_livingmatter    )%varname = 's_bio_livingmatter'
   cffields( id_s_horLO_deadmatter    )%varname = 's_horLO_deadmatter'
   cffields( id_s_verLO_deadmatter    )%varname = 's_verLO_deadmatter'
   cffields( id_s_horAD_deadmatter    )%varname = 's_horAD_deadmatter'
   cffields( id_s_verAD_deadmatter    )%varname = 's_verAD_deadmatter'
   cffields( id_s_diffV_deadmatter    )%varname = 's_diffV_deadmatter'
   cffields( id_s_diffH_deadmatter    )%varname = 's_diffH_deadmatter'
   cffields( id_s_sink_deadmatter     )%varname = 's_sink_deadmatter'
   cffields( id_s_bio_deadmatter      )%varname = 's_bio_deadmatter'
   cffields( id_s_benthos_alk         )%varname = 's_benthos_alk'
   cffields( id_s_benthos_dic         )%varname = 's_benthos_dic'
   cffields( id_s_benthos_livingmatter)%varname = 's_benthos_livingmatter'
   cffields( id_s_benthos_deadmatter  )%varname = 's_benthos_deadmatter'
   cffields( id_m_alk                 )%varname = 'm_alk'
   cffields( id_m_dic                 )%varname = 'm_dic'
   cffields( id_m_livingmatter        )%varname = 'm_livingmatter'
   cffields( id_m_deadmatter          )%varname = 'm_deadmatter'
   cffields( id_s_asml_alk            )%varname = 's_asml_alk'
   cffields( id_s_asml_dic            )%varname = 's_asml_dic'
   cffields( id_s_asml_livingmatter   )%varname = 's_asml_livingmatter'
   cffields( id_s_asml_deadmatter     )%varname = 's_asml_deadmatter'
   cffields( id_s_surf_alk            )%varname = 's_surf_alk'
   cffields( id_s_surf_dic            )%varname = 's_surf_dic'
   cffields( id_s_vol_alk             )%varname = 's_vol_alk'
   cffields( id_s_vol_dic             )%varname = 's_vol_dic'
   cffields( id_s_vol_livingmatter    )%varname = 's_vol_livingmatter'
   cffields( id_s_vol_deadmatter      )%varname = 's_vol_deadmatter'
   
   ! _________________________________
   ! ___ allocate and intialize fields
   
   ! directed fluxes through area
   DO i=1, size(cffieldsflux)
      ids=cffieldsflux(i)
      
      allocate(cffields(ids)%instantconc  (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%timemeanconc (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmconc     (nlmax,myDim_nod2D), source=0.0)
      
      cffields(ids)%unitsconc = 'mmol m$^{-2}$ s$^{-1}$'
   ENDDO ! (directed fluxes)
   
   ! sources minus sinks in model
   DO i=1, size(cffieldssms)
      ids=cffieldssms(i)
      
      ! default case: not horizontal diffusion or advection
      IF (ALL(ids /= [id_s_diffH_alk, id_s_diffH_deadmatter  , &
                      id_s_diffH_dic, id_s_diffH_livingmatter, &
                      id_s_horLO_alk, id_s_horLO_deadmatter  , &
                      id_s_horLO_dic, id_s_horLO_livingmatter, &
                      id_s_horAD_alk, id_s_horAD_deadmatter  , &
                      id_s_horAD_dic, id_s_horAD_livingmatter])) THEN
      
         allocate(cffields(ids)%instantconc  (nlmax,myDim_nod2D), source=0.0)
         allocate(cffields(ids)%instantmass  (nlmax,myDim_nod2D), source=0.0)
         
      ! cross-edge horizontal diffusion and advection requires halo edges 
      ELSE
                        
         allocate(cffields(ids)%instantconc (nlmax,myDim_nod2D+eDim_nod2D), source=0.0)
         allocate(cffields(ids)%instantmass (nlmax,myDim_nod2D+eDim_nod2D), source=0.0)
      ENDIF
      
      allocate(cffields(ids)%timemeanconc (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmconc     (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%timemeanmass (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmmass     (nlmax,myDim_nod2D), source=0.0)
      
      cffields(ids)%unitsconc = 'mmol m$^{-3}$ s$^{-1}$'
      cffields(ids)%unitsmass = 'mmol s$^{-1}$'
   ENDDO ! (sources minus sinks in model)

   ! tracer fields
   DO i=1, size(cffieldstracer)
      ids=cffieldstracer(i)
      
      allocate(cffields(ids)%instantconc  (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%instantmass  (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmconc     (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmmass     (nlmax,myDim_nod2D), source=0.0)
      
      cffields(ids)%unitsconc = 'mmol m$^{-3}$'
      cffields(ids)%unitsmass = 'mmol'
   ENDDO ! (tracer fields)
   
   ! sources minus sinks during assimilation step
   DO i=1, size(cffieldsasml)
      ids=cffieldsasml(i)
      
      allocate(cffields(ids)%fconc        (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%fmass        (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%aconc        (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%amass        (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmconc     (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmmass     (nlmax,myDim_nod2D), source=0.0)
      
      cffields(ids)%unitsconc = 'mmol m$^{-3}$ s$^{-1}$'
      cffields(ids)%unitsmass = 'mmol s$^{-1}$'

   ENDDO ! (assimilation fields)
   
   allocate(factor_mass (nlmax,myDim_nod2D), source=0.0) ! used in prepoststep
   allocate(factor_conc (nlmax,myDim_nod2D), source=0.0) ! used in prepoststep
   
   ! change of concentration due to cell volume change
   DO i=1, size(cffieldsvol)
      ids=cffieldsvol(i)
      
      allocate(cffields(ids)%instantconc  (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%timemeanconc (nlmax,myDim_nod2D), source=0.0)
      allocate(cffields(ids)%ensmconc     (nlmax,myDim_nod2D), source=0.0)
      
      cffields(ids)%unitsconc = 'mmol m$^{-3}$'
   ENDDO ! (cell volume change)

END SUBROUTINE init_carbonfluxes_diags_arrays


! ***********************************************
! ***                                         ***
! ***   carbonfluxes_diags_output_timemean     ***
! ***                                         ***
! ***********************************************
! - compute time mean and ensemble means of model fields and write output
! - called from model (fvom_main)

SUBROUTINE carbonfluxes_diags_output_timemean(mstep)

      USE g_clock, &
        ONLY: month, num_day_in_month, fleapyear, cyearnew, &
        daynew, timenew
      USE g_events
      USE g_config, &
        ONLY: step_per_day
      USE g_parsup, &
        ONLY: myDim_nod2D
      USE mod_assim_pdaf, &
        ONLY: nlmax, dim_ens, mesh_fesom
      USE mod_parallel_pdaf, &
        ONLY: mype_world, abort_parallel, task_id, mype_submodel, &
        COMM_COUPLE, filterpe, writepe, mype_model
      USE g_comm_auto, &
        ONLY: gather_nod
      USE o_arrays, &
        ONLY: tr_arr, hnode_new


      IMPLICIT NONE
      include 'mpif.h'
      
      ! arguments
      INTEGER, INTENT(in) :: mstep
      ! local variables
      LOGICAL :: now_to_write
      INTEGER :: month_iter, whichmonth,fleap
      REAL    :: weights
      INTEGER :: mpierror
      character(len=200) :: filename
      REAL, allocatable ::  f_mass(:,:)
      
      integer :: i, ids ! counters
      
      ! debugging:
      REAL, allocatable  :: data3_g(:,:)                 ! Temporary array for global 3D-fields
      character(len=3)   :: cday                         ! Day as character string
      character(len=5)   :: ctime                        ! Time as character string
      integer            :: fileid                       ! nc-file ID for output file
                  
      ! *** Debugging *** (write instantanteously collected fields for 1 ensemble member)
      IF (cfdiags_debug .and. filterpe) THEN
         allocate(data3_g(nlmax,mesh_fesom%nod2D))
         write(cday ,'(i3)') daynew
         write(ctime,'(i5)') int(timenew)
         
         DO i=1, cfnfields
            ! only model fields
            if (ALL(i /= cffieldsasml)) then
            ! concentration
            if (writepe) open(unit=1, file=cyearnew//cday//ctime//trim(cffields(i)%varname)//'_conc.txt', status='unknown')
            CALL gather_nod(cffields(i)%instantconc, data3_g)
            if (writepe) write(1,*) data3_g
            if (writepe) close(1)
            ! mass
            if (writepe) open(unit=1, file=cyearnew//cday//ctime//trim(cffields(i)%varname)//'_mass.txt', status='unknown')
            CALL gather_nod(cffields(i)%instantmass, data3_g)
            if (writepe) write(1,*) data3_g
            if (writepe) close(1)
            endif ! only model fields
         ENDDO ! i=1, cfnfields
      
         deallocate(data3_g)
      ENDIF ! (cfdiags_debug .and. filterpe)
      

      ! *** set output frequency ***
      now_to_write = .FALSE.
      IF (cfoutfreq_unit=='m') then
         call monthly_event(now_to_write,cfoutfreq)
      ELSEIF (cfoutfreq_unit=='d') then
         call daily_event(now_to_write,cfoutfreq)
      ELSE
         if (writepe) WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'CFDIAGS: Invalid output frequency'
      ENDIF
      ! stepwise debugging of output freq.:
      IF (cfdiags_debug .and. writepe) WRITE (*, '(/a, 1x, a, 1x, i5, 1x, a, 1x, l)') 'FESOM-PDAF', 'Call carbonflux diagnostics at step', mstep, 'nowtowrite', now_to_write
      
      ! ***
      ! *** add instantaneous data to timemean
      ! ***
      
      ! directed fluxes through area
      DO i=1, size(cffieldsflux)
         ids = cffieldsflux(i)
         cffields(ids)%timemeanconc = cffields(ids)%timemeanconc + cffields(ids)%instantconc(:,:myDim_nod2D)
      ENDDO
      ! sources minus sinks in model
      DO i=1, size(cffieldssms)
         ids=cffieldssms(i)
         cffields(ids)%timemeanconc = cffields(ids)%timemeanconc + cffields(ids)%instantconc(:,:myDim_nod2D)
         cffields(ids)%timemeanmass = cffields(ids)%timemeanmass + cffields(ids)%instantmass(:,:myDim_nod2D)
      ENDDO
      ! cell volume change
      DO i=1, size(cffieldsvol)
         ids = cffieldsvol(i)
         cffields(ids)%timemeanconc = cffields(ids)%timemeanconc + cffields(ids)%instantconc(:,:myDim_nod2D)
      ENDDO
      
      IF (now_to_write) THEN
      ! ***
      ! *** compute time mean, compute ensemble mean, write output and reset time mean to zero
      ! ***
      
      ! apply weighting to mean
      IF (cfoutfreq_unit=='m') then
         weights = 1.0/REAL(num_day_in_month(fleapyear,month)*step_per_day)/REAL(dim_ens)
      ELSEIF (cfoutfreq_unit=='d') then
         weights = 1.0/REAL(step_per_day)/REAL(dim_ens)
      ENDIF
      
      ! reduce ensemble to ensemble mean
      ! directed fluxes through area
      DO i=1, size(cffieldsflux)
         ids = cffieldsflux(i)
         cffields(ids)%timemeanconc = cffields(ids)%timemeanconc * weights
         !                - send -                   - receive -            - size -             - type -              - sum -  -receiver-   -ensemble-    -check-
         CALL MPI_REDUCE( cffields(ids)%timemeanconc, cffields(ids)%ensmconc, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
      ENDDO
      ! sources minus sinks in model
      DO i=1, size(cffieldssms)
         ids=cffieldssms(i)
         cffields(ids)%timemeanconc = cffields(ids)%timemeanconc * weights
         cffields(ids)%timemeanmass = cffields(ids)%timemeanmass * weights
         !                - send -                   - receive -            - size -             - type -              - sum -  -receiver-   -ensemble-    -check-
         CALL MPI_REDUCE( cffields(ids)%timemeanconc, cffields(ids)%ensmconc, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( cffields(ids)%timemeanmass, cffields(ids)%ensmmass, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
      ENDDO
      ! cell volume change
      DO i=1, size(cffieldsvol)
         ids = cffieldsvol(i)
         cffields(ids)%timemeanconc = cffields(ids)%timemeanconc * weights
         !                - send -                   - receive -            - size -             - type -              - sum -  -receiver-   -ensemble-    -check-
         CALL MPI_REDUCE( cffields(ids)%timemeanconc, cffields(ids)%ensmconc, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
      ENDDO
      
      ! tracer fields
      ! collect snapshots
      cffields(id_m_dic)         %instantconc =   tr_arr(:nlmax,:myDim_nod2D, 4)    ! DIC
      cffields(id_m_alk)         %instantconc =   tr_arr(:nlmax,:myDim_nod2D, 5)    ! Alk
      cffields(id_m_livingmatter)%instantconc = ( tr_arr(:nlmax,:myDim_nod2D, 7) &  ! PhyC
                                                + tr_arr(:nlmax,:myDim_nod2D,12) &  ! HetC
                                                + tr_arr(:nlmax,:myDim_nod2D,22) &  ! PhyCalc
                                                + tr_arr(:nlmax,:myDim_nod2D,16) &  ! DiaC
                                                + tr_arr(:nlmax,:myDim_nod2D,26))   ! Zoo2C
      cffields(id_m_deadmatter)  %instantconc = ( tr_arr(:nlmax,:myDim_nod2D,28) &  ! Det2C
                                                + tr_arr(:nlmax,:myDim_nod2D,30) &  ! Det2Calc
                                                + tr_arr(:nlmax,:myDim_nod2D,10) &  ! DetC
                                                + tr_arr(:nlmax,:myDim_nod2D,23) &  ! DetCalc
                                                + tr_arr(:nlmax,:myDim_nod2D,14))   ! DOC
      ! concentration to mass
      allocate(f_mass(nlmax,myDim_nod2D))
      f_mass = mesh_fesom%areasvol(:nlmax,:myDim_nod2D) * hnode_new(:nlmax,:myDim_nod2D)
      DO i=1, size(cffieldstracer)
         ids=cffieldstracer(i)
         cffields(ids)%instantmass = cffields(ids)%instantconc * f_mass
      ENDDO
      ! ensemble mean
      DO i=1, size(cffieldstracer)
         ids=cffieldstracer(i)
         cffields(ids)%instantconc = cffields(ids)%instantconc / REAL(dim_ens)
         cffields(ids)%instantmass = cffields(ids)%instantmass / REAL(dim_ens)
         !                - send -                   - receive -            - size -             - type -              - sum -  -receiver-   -ensemble-    -check-
         CALL MPI_REDUCE( cffields(ids)%instantconc, cffields(ids)%ensmconc, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
         CALL MPI_REDUCE( cffields(ids)%instantmass, cffields(ids)%ensmmass, nlmax * myDim_nod2D, MPI_DOUBLE_PRECISION, MPI_SUM, 0          , COMM_COUPLE, mpierror)
      ENDDO

      ! write ensemble mean output
      IF (filterpe) call write_carbonfluxes_diags_out()
      
      ! reset timemean data to zero
      ! directed fluxes through area
      DO i=1, size(cffieldsflux)
         ids = cffieldsflux(i)
         cffields(ids)%timemeanconc = 0.0
      ENDDO
      ! sources minus sinks in model
      DO i=1, size(cffieldssms)
         ids=cffieldssms(i)
         cffields(ids)%timemeanconc = 0.0
         cffields(ids)%timemeanmass = 0.0
      ENDDO
      ! cell volumne
      DO i=1, size(cffieldsvol)
         ids = cffieldsvol(i)
         cffields(ids)%timemeanconc = 0.0
      ENDDO
      
      ENDIF ! now_to_write
  
END SUBROUTINE carbonfluxes_diags_output_timemean


! ****************************************************
! ***                                              ***
! ***   carbonfluxes_diags_output_timemean_asml     ***
! ***                                              ***
! ****************************************************
! - compute time mean of SMS during assimilation step and write output
! - called from PDAF (prepoststep)

SUBROUTINE carbonfluxes_diags_output_timemean_asml()

      USE g_clock, &
        ONLY: month, num_day_in_month, fleapyear, cyearnew, &
        daynew, timenew
      USE g_events
      USE g_config, &
        ONLY: step_per_day
      USE g_parsup, &
        ONLY: myDim_nod2D
      USE mod_assim_pdaf, &
        ONLY: nlmax, dim_ens, mesh_fesom, assim_time
      USE mod_parallel_pdaf, &
        ONLY: mype_world, abort_parallel, task_id, mype_submodel, &
        COMM_COUPLE, filterpe, writepe
      USE g_comm_auto, &
        ONLY: gather_nod


      IMPLICIT NONE
      include 'mpif.h'
      
      ! local variables
      LOGICAL :: now_to_write
      INTEGER :: month_iter, whichmonth,fleap
      REAL    :: weights
      INTEGER :: mpierror
      character(len=200) :: filename
      
      integer :: i, ids ! counters
      
      ! debugging:
      REAL, allocatable  :: data3_g(:,:)                 ! Temporary array for global 3D-fields
      character(len=3)   :: cday                         ! Day as character string
      character(len=5)   :: ctime                        ! Time as character string
      integer            :: fileid                       ! nc-file ID for output file
      
      ! *** set output frequency ***
      now_to_write = .FALSE.
      IF (cfoutfreq_unit=='m') then
         call monthly_event_assimstep(now_to_write,cfoutfreq)
      ELSEIF (cfoutfreq_unit=='d') then
         if (mod(daynew,cfoutfreq)==0 .and. timenew==assim_time) then
            now_to_write=.true.
         else
            now_to_write=.false.
         endif
      ENDIF
      
      ! ***
      ! *** add instantaneous data to time mean, which is already ensemble mean
      
      DO i=1, size(cffieldsasml)
         ids = cffieldsasml(i)
         cffields(ids)%ensmconc = cffields(ids)%ensmconc + (cffields(ids)%aconc - cffields(ids)%fconc)
         cffields(ids)%ensmmass = cffields(ids)%ensmmass + (cffields(ids)%amass - cffields(ids)%fmass)
      ENDDO
      
      IF (now_to_write) THEN
      ! ***
      ! *** compute time mean, write output and reset time mean to zero
      ! ***
      
      ! apply weighting to mean; one assimilation step per day
      IF (cfoutfreq_unit=='m') then
         weights = 1.0/REAL(num_day_in_month(fleapyear,month))
      ELSEIF (cfoutfreq_unit=='d') then
         weights = 1.0
      ENDIF
      DO i=1, size(cffieldsasml)
         ids = cffieldsasml(i)
         cffields(ids)%ensmconc = cffields(ids)%ensmconc * weights
         cffields(ids)%ensmmass = cffields(ids)%ensmmass * weights
      ENDDO 

      ! write time mean output
      IF (filterpe) call write_carbonfluxes_diags_out_asml()
      
      ! reset time mean to zero
      DO i=1, size(cffieldsasml)
         ids = cffieldsasml(i)
         cffields(ids)%ensmconc = 0.0
         cffields(ids)%ensmmass = 0.0
      ENDDO  
      ENDIF ! now_to_write
  
END SUBROUTINE carbonfluxes_diags_output_timemean_asml


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
         ONLY: cyearnew, num_day_in_month, fleapyear, &
         yearold, yearnew, month, daynew, timenew
      USE g_parsup, &
         ONLY: myDim_nod2D
      USE g_comm_auto, &
         ONLY: gather_nod
      USE mod_parallel_pdaf, &
         ONLY: writepe
      USE netcdf

      IMPLICIT NONE
      
      character(len=200) :: filename
      character(len=200) :: units
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

      LOGICAL :: file_exists
      
      integer :: i, ids ! counters

      filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.cfx.'//cyearnew//'.nc'

      INQUIRE(file=filename, exist=file_exists)

      IF (.not. file_exists) THEN
      
      ! gather GEO coordinates (from all PEs)
      allocate(lon(mesh_fesom%nod2D),lat(mesh_fesom%nod2D))
      call gather_nod(mesh_fesom%geo_coord_nod2D(1, 1:myDim_nod2D), lon)
      call gather_nod(mesh_fesom%geo_coord_nod2D(2, 1:myDim_nod2D), lat)
      
      IF (writepe) THEN
      ! initialize file
      WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize carbon flux diags NetCDF file: '//TRIM(filename)
      ! open file
      call check(NF90_CREATE(trim(filename),NF90_NETCDF4,fileid))
      
      ! define dimensions
      call check( NF90_DEF_DIM(fileid, 'nod2', mesh_fesom%nod2d, dimID_nod2))
      call check( NF90_DEF_DIM(fileid, 'nz1',  nlmax,            dimID_nz1))
      call check( NF90_DEF_DIM(fileid, 'time', NF90_UNLIMITED,   dimID_iter))
      
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
      
      write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearnew, '-', 1, '-', 1, ' 0:0:0'
      call check( nf90_put_att(fileid, varid_iter, 'long_name', 'time'))
      call check( nf90_put_att(fileid, varid_iter, 'standard_name', 'time'))
      call check( nf90_put_att(fileid, varid_iter, 'units', trim(att_text)))
      call check( nf90_put_att(fileid, varid_iter, 'axis', 'T'))
      call check( nf90_put_att(fileid, varid_iter, 'stored_direction', 'increasing'))
      
      ! fill dimension variables
      call check (nf90_enddef(fileid))
      call check (nf90_put_var(fileid, varid_nz1,  mesh_fesom%Z   (1:nlmax)))
      call check (nf90_put_var(fileid, varid_nod2, [(n,n=1,mesh_fesom%nod2D)]))
      
      call check (nf90_put_var(fileid, varid_lon,  REAL(180./pi * lon, 4)))
      call check (nf90_put_var(fileid, varid_lat,  REAL(180./pi * lat, 4)))
      
      ! define diagnostic variables
      call check (nf90_redef(fileid))

      ! set dimensions
      ndims = 3
      dimIDs(1) = dimID_nod2
      dimIDs(2) = dimID_nz1
      dimIDs(3) = dimID_iter
      
      DO ids=1,cfnfields
         ! define concentrations
         IF (outconc) then
            call check( NF90_DEF_VAR(fileid, trim(cffields(ids)%varname)//'_c', NF90_FLOAT, dimIDs(1:ndims), varid))
            call check( nf90_put_att(fileid, varid, 'units',trim(cffields(ids)%unitsconc)))
         endif
         ! define mass
         if (outmass) then
            if (ANY(ids == cffieldsasml  )  .or. &
                ANY(ids == cffieldssms   )  .or. &
                ANY(ids == cffieldstracer)) then
            call check( NF90_DEF_VAR(fileid, trim(cffields(ids)%varname)//'_m', NF90_FLOAT, dimIDs(1:ndims), varid))
            call check( nf90_put_att(fileid, varid, 'units',trim(cffields(ids)%unitsmass)))
            endif
         endif
      ENDDO
      
      call check(NF90_ENDDEF(fileid))
      call check(NF90_CLOSE(fileid))
      
      ENDIF ! writepe

      deallocate(lat,lon)

      ENDIF ! .not. file_exists

END SUBROUTINE init_carbonfluxes_diags_out

! ***********************************************
! ***                                         ***
! ***   write_carbonfluxes_diags_out          ***
! ***                                         ***
! ***********************************************
! Writes carbon flux diagnostics to netCDF during model step
SUBROUTINE write_carbonfluxes_diags_out()

USE g_clock, &
   ONLY: month, cyearnew, daynew, num_day_in_month, fleapyear
USE recom_config, ONLY: SecondsPerDay
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

! LOCAL VARIABLES:
REAL, allocatable  :: data3_g(:,:)                 ! Temporary array for global 3D-fields
character(len=200) :: filename                     ! Full name of output file
integer            :: fileid                       ! nc-file ID for output file
INTEGER            :: writepos                     ! Write position
INTEGER            :: writetime                    ! Model time in seconds since beginning of year
INTEGER            :: varid_time
character(len=100) :: varname
integer            :: i, ids                       ! counters


! print screen information
filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.cfx.'//cyearnew//'.nc'
IF (writepe) THEN
   WRITE (*, '(a, 8x, a, i9, a, a)') 'FESOM-PDAF', 'Write carbon flux diagnostics to NetCDF during model step at day ', &
   daynew, ' to file ', TRIM(filename)
   ! open file
   call check( nf90_open(TRIM(filename), nf90_write, fileid))
ENDIF ! writepe

! specify time coordinate
IF (cfoutfreq_unit=='m') THEN
   writepos  = month/cfoutfreq
   writetime = (daynew+1-num_day_in_month(fleapyear,month))*SecondsPerDay
ELSEIF (cfoutfreq_unit=='d') THEN
   writepos  = daynew/cfoutfreq
   writetime = (daynew-1)*SecondsPerDay
ELSE
   if (writepe) WRITE (*, '(a, 8x, a)') 'FESOM-PDAF', 'CFDIAGS Invalid output frequency.'
ENDIF

if (writepe) call check( nf90_inq_varid(fileid, "time"  , varid_time))
if (writepe) call check( nf90_put_var  (fileid, varid_time, writetime, start=(/ writepos /)))

! gather and write global ocean fields
allocate(data3_g(nlmax,mesh_fesom%nod2D))

! directed fluxes through area
DO i=1, size(cffieldsflux)
 ids = cffieldsflux(i)
 ! conc
 if (outconc) then
    varname = trim(cffields(ids)%varname)//'_c'
    CALL gather_nod(cffields(ids)%ensmconc, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
ENDDO
! sources minus sinks in model
DO i=1, size(cffieldssms)
 ids=cffieldssms(i)
 ! conc
 if (outconc) then
    varname = trim(cffields(ids)%varname)//'_c'
    CALL gather_nod(cffields(ids)%ensmconc, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
 ! mass
 if (outmass) then
    varname = trim(cffields(ids)%varname)//'_m'
    CALL gather_nod(cffields(ids)%ensmmass, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
ENDDO
! tracer fields
DO i=1, size(cffieldstracer)
 ids=cffieldstracer(i)
 ! conc
 if (outconc) then
    varname = trim(cffields(ids)%varname)//'_c'
    CALL gather_nod(cffields(ids)%ensmconc, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
 ! mass
 if (outmass) then
    varname = trim(cffields(ids)%varname)//'_m'
    CALL gather_nod(cffields(ids)%ensmmass, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
ENDDO
! volume change
DO i=1, size(cffieldsvol)
 ids = cffieldsvol(i)
 ! conc
 if (outconc) then
    varname = trim(cffields(ids)%varname)//'_c'
    CALL gather_nod(cffields(ids)%ensmconc, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
ENDDO

deallocate(data3_g)

! close file:
IF (writepe) call check (nf90_close(fileid))
END SUBROUTINE write_carbonfluxes_diags_out


! ***********************************************
! ***                                         ***
! ***   write_carbonfluxes_diags_out_asml     ***
! ***                                         ***
! ***********************************************
! Writes carbon flux diagnostics to netCDF during assimilation step
SUBROUTINE write_carbonfluxes_diags_out_asml()

USE g_clock, &
   ONLY: month, cyearnew, daynew, num_day_in_month, fleapyear
USE recom_config, ONLY: SecondsPerDay
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

! LOCAL VARIABLES:
REAL, allocatable  :: data3_g(:,:)                 ! Temporary array for global 3D-fields
character(len=200) :: filename                     ! Full name of output file
integer            :: fileid                       ! nc-file ID for output file
INTEGER            :: writepos                     ! Write position
character(len=100) :: varname
integer            :: i, ids                       ! counters


! print screen information
filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.cfx.'//cyearnew//'.nc'
IF (writepe) THEN
   WRITE (*, '(a, 8x, a, i9, a, a)') 'FESOM-PDAF', 'Write carbon flux diagnostics to NetCDF during assimilation step at day ', &
   daynew, ' to file ', TRIM(filename)
   ! open file
   call check( nf90_open(TRIM(filename), nf90_write, fileid))
ENDIF ! writepe

! specify time coordinate
IF (cfoutfreq_unit=='m') THEN
   writepos  = month/cfoutfreq
ELSEIF (cfoutfreq_unit=='d') THEN
   writepos  = daynew/cfoutfreq
ENDIF

! sources minus sinks during assimilation step
allocate(data3_g(nlmax,mesh_fesom%nod2D))

DO i=1, size(cffieldsasml)
 ids=cffieldsasml(i)
 ! conc
 if (outconc) then
    varname = trim(cffields(ids)%varname)//'_c'
    CALL gather_nod(cffields(ids)%ensmconc, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
 ! mass
 if (outmass) then
    varname = trim(cffields(ids)%varname)//'_m'
    CALL gather_nod(cffields(ids)%ensmmass, data3_g) 
    IF (writepe) call putvar(fileid,varname,data3_g,writepos)
 endif
ENDDO

deallocate(data3_g)

! close file:
IF (writepe) call check (nf90_close(fileid))
END SUBROUTINE write_carbonfluxes_diags_out_asml


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
                                   count=(/ mesh_fesom%nod2D, nlmax, 1 /) ))
                                   ! dims: 1-nod2, 2-nz / nz1, 3-iter
END SUBROUTINE putvar

! ***********************************************
! ***                                         ***
! ***   cfdiags_computetransport              ***
! ***                                         ***
! ***********************************************
SUBROUTINE cfdiags_computetransport(tr_arr,Unode,wvel)

USE mod_assim_pdaf, &
   ONLY: nlmax, mesh_fesom
USE g_parsup, &
   ONLY: myDim_nod2D, eDim_nod2D
USE o_param, &
   ONLY: num_tracers

IMPLICIT NONE

! ARGUMENTS:
REAL, intent(in)  :: tr_arr(mesh_fesom%nl-1, myDim_nod2D+eDim_nod2D, num_tracers) ! tracers
REAL, intent(in)  :: Unode(2, mesh_fesom%nl-1, myDim_nod2D+eDim_nod2D)            ! horizontal velocities on nodes
REAL, intent(in)  :: wvel(mesh_fesom%nl, myDim_nod2D+eDim_nod2D)                  ! vertical velocity on levels

! LOCAL VARIABLES:
real, allocatable    :: wlayers(:,:)   ! vertical velocity on layers
integer, allocatable :: tracerlist(:)
integer              :: trcounter

! vertical velocities at layers
allocate(wlayers(nlmax,myDim_nod2D))
wlayers = 0.5*(wvel(1:nlmax,:myDim_nod2D) + wvel(2:nlmax+1,:myDim_nod2D)) ! positive upwards

! alkalinity (5)
cffields(id_t_u_alk)%instantconc = tr_arr(:nlmax,:myDim_nod2D,5) * Unode(1,:nlmax,:myDim_nod2D)
cffields(id_t_v_alk)%instantconc = tr_arr(:nlmax,:myDim_nod2D,5) * Unode(2,:nlmax,:myDim_nod2D)
cffields(id_t_w_alk)%instantconc = tr_arr(:nlmax,:myDim_nod2D,5) * wlayers(:nlmax,:myDim_nod2D)

! DIC (4)
cffields(id_t_u_dic)%instantconc = tr_arr(:nlmax,:myDim_nod2D,4) * Unode(1,:nlmax,:myDim_nod2D)
cffields(id_t_v_dic)%instantconc = tr_arr(:nlmax,:myDim_nod2D,4) * Unode(2,:nlmax,:myDim_nod2D)
cffields(id_t_w_dic)%instantconc = tr_arr(:nlmax,:myDim_nod2D,4) * wlayers(:nlmax,:myDim_nod2D)

! Living biomass
allocate(tracerlist(5))
tracerlist(1) =  7 ! PhyC
tracerlist(2) = 12 ! HetC
tracerlist(3) = 22 ! PhyCalc
tracerlist(4) = 16 ! DiaC
tracerlist(5) = 26 ! Zoo2C

DO trcounter=1,5
   cffields(id_t_u_livingmatter)%instantconc = cffields(id_t_u_livingmatter)%instantconc + tr_arr(:nlmax,:myDim_nod2D,tracerlist(trcounter)) * Unode(1,:nlmax,:myDim_nod2D)
   cffields(id_t_v_livingmatter)%instantconc = cffields(id_t_v_livingmatter)%instantconc + tr_arr(:nlmax,:myDim_nod2D,tracerlist(trcounter)) * Unode(2,:nlmax,:myDim_nod2D)
   cffields(id_t_w_livingmatter)%instantconc = cffields(id_t_w_livingmatter)%instantconc + tr_arr(:nlmax,:myDim_nod2D,tracerlist(trcounter)) * wlayers(:nlmax,:myDim_nod2D)
ENDDO
deallocate(tracerlist)

! Dead biomass
allocate(tracerlist(5))
tracerlist(1) = 28 ! DetZ2C
tracerlist(2) = 30 ! DetZ2Calc
tracerlist(3) = 10 ! DetC
tracerlist(4) = 23 ! DetCalc
tracerlist(5) = 14 ! DOC

DO trcounter=1,5
   cffields(id_t_u_deadmatter)%instantconc = cffields(id_t_u_deadmatter)%instantconc + tr_arr(:nlmax,:myDim_nod2D,tracerlist(trcounter)) * Unode(1,:nlmax,:myDim_nod2D)
   cffields(id_t_v_deadmatter)%instantconc = cffields(id_t_v_deadmatter)%instantconc + tr_arr(:nlmax,:myDim_nod2D,tracerlist(trcounter)) * Unode(2,:nlmax,:myDim_nod2D)
   cffields(id_t_w_deadmatter)%instantconc = cffields(id_t_w_deadmatter)%instantconc + tr_arr(:nlmax,:myDim_nod2D,tracerlist(trcounter)) * wlayers(:nlmax,:myDim_nod2D)
ENDDO
deallocate(tracerlist)

deallocate(wlayers)
END SUBROUTINE cfdiags_computetransport


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


