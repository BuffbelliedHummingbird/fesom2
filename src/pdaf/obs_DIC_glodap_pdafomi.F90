!$Id: obs_TSprof_EN4_pdafomi.F90 2543 2021-05-13 08:17:31Z lnerger $
!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
!! Module type here: temperature and salinity subsurface profiles data from EN4.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs)
!! and for the observation operator (obs_op).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! In addition, there are two optional routine, which are required if filters 
!! with localization are used:
!! * init_dim_obs_l_TYPE \n
!!           Only required if domain-localized filters (e.g. LESTKF, LETKF) are used:
!!           Count number of local observations of module-type according to
!!           their coordinates (distance from local analysis domain). Initialize
!!           module-internal distances and index arrays.
!! * localize_covar_TYPE \n
!!           Only required if the localized EnKF is used:
!!           Apply covariance localization in the LEnKF.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE obs_DIC_glodap_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter, writepe ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l         ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_o_DIC_glodap   !< Whether to assimilate temperature profiles

  ! Further variables specific for the EN4 profile observations
  CHARACTER(len=110) :: path_obs_DIC_glodap  = ''      !< Path to profile observations
  CHARACTER(len=110) :: file_DIC_glodap

  REAL    :: rms_obs_DIC_glodap      ! Observation error
  REAL    :: bias_obs_DIC_glodap     ! Observation bias

  REAL    :: lradius_DIC_glodap      ! Localization radius
  REAL    :: sradius_DIC_glodap      ! Support radius for localization function

  REAL    :: DIC_glodap_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)

  REAL, ALLOCATABLE :: mean_DIC_p (:)             ! ensemble mean for observation exclusion
  REAL, ALLOCATABLE :: loc_radius_DIC_glodap(:)   ! localization radius array


! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in INIT_DIM_OBS
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!                                           ! (0) Cartesian, (1) Cartesian periodic
!                                           ! (2) simplified geographic, (3) geographic haversine function
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! Indices of observed field in state vector (process-local)
!           
!           Optional variables - they can be set in INIT_DIM_OBS
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!
!           Variables with predefined values - they can be changed in INIT_DIM_OBS
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!      INTEGER :: use_global_obs=1          ! Whether to use (1) global full obs. 
!                                           ! or (0) obs. restricted to those relevant for a process domain
!
!           The following variables are set in the routine PDAFomi_gather_obs
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!      INTEGER :: off_obs_g                 ! Offset of this observation in overall global obs. vector
!      INTEGER :: obsid                     ! Index of observation over all assimilated observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!      INTEGER :: locweight                 ! Specify localization function
!      REAL :: lradius                      ! localization radius
!      REAL :: sradius                      ! support radius for localization function
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF (default: 0=Gaussian)
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: 1=use global full observations)
!!
!! Further variables are set when the routine PDAFomi_gather_obs is called.
!!
  SUBROUTINE init_dim_obs_DIC_glodap(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs, PDAFomi_set_debug_flag
    USE mod_assim_pdaf, &
         ONLY: offset, use_global_obs, id, &
               mesh_fesom, &
               local_range, srange
    USE mod_parallel_pdaf, &
         ONLY: MPI_SUM, MPIerr, COMM_filter, MPI_INTEGER
    USE g_parsup, &
         ONLY: myDim_nod2D
    USE g_clock, &
         ONLY: month, day_in_month, yearold, timenew
    USE o_param, &
         ONLY: pi
         
    USE netcdf

    IMPLICIT NONE

!~     INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in)    :: step      !< Current time step
    INTEGER, INTENT(inout) :: dim_obs   !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, s, k, j          ! Counters
    INTEGER :: dim_obs_p                      ! Number of process local observations
    CHARACTER(len=2) :: mype_string           ! String for process rank
    CHARACTER(len=4) :: year_string           ! String for yearly observation data path
    CHARACTER(len=2) :: mon_string            ! String for daily observation files
    CHARACTER(len=2) :: day_string            ! String for daily observation files
    
    REAL, ALLOCATABLE :: obs_p(:)             ! PE-local observed observation values
    REAL, ALLOCATABLE :: ocoord_p(:,:)        ! PE-local coordinates of observations
    REAL, ALLOCATABLE :: lon_p(:),lat_p(:)
    REAL, ALLOCATABLE :: ivariance_obs_p(:)   ! PE-local array of observation errors
    INTEGER, ALLOCATABLE :: nod1_p(:), &
                            nod2_p(:), &
                            nod3_p(:)         ! Array of observation pe-local indeces on FESOM grid
    INTEGER, ALLOCATABLE :: nl_p(:)           ! Array of observation layer indeces
    
    INTEGER :: ncstat                         ! Status for NetCDF functions
    INTEGER :: ncid, dimid                    ! NetCDF IDs
    INTEGER :: id_obs, id_nod1, id_nod2, &
               id_nod3, id_nl, id_lon, &
               id_lat

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    ! Store whether to assimilate this observation type
    IF (assim_o_DIC_glodap) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Initialize flag for type of full observations
    thisobs%use_global_obs = use_global_obs

    ! set localization radius, everywhere the same
    lradius_DIC_glodap = local_range
    sradius_DIC_glodap = srange


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Initialize complete file name
    WRITE(mype_string,'(i2.2)') mype_filter
    WRITE(year_string,'(i4.4)') yearold
    WRITE(mon_string, '(i2.2)') month
    WRITE(day_string, '(i2.2)') day_in_month
    
    file_DIC_glodap = TRIM(year_string//'-'//mon_string//'-'//day_string//'.'//mype_string//'.nc')
    
    ! Debugging message:
    IF (mype_filter == 0) THEN
       WRITE (*,'(a,5x,a,i2,a,i2,a,i4,a,f5.2,a,a)') &
            'FESOM-PDAF', 'Assimilate GLODAP DIC observations - OBS_DIC_glodap at ', &
            day_in_month, '.', month, '.', yearold, ' ', timenew/3600.0,&
            ' h; read from file: ', file_DIC_glodap
    END IF
    
    ! Open the pe-local NetCDF file for that day
    ncstat = nf90_open(TRIM(path_obs_DIC_glodap)//year_string//'/dist72/'//file_DIC_glodap, nf90_nowrite, ncid)
    if (ncstat /= nf90_noerr) then
     print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error opening NetCDF file'
    end if
    
    ! Get the number of observations
    ! 1. Get the dimension ID
    ncstat = nf90_inq_dimid(ncid, "index", dimid)
    if (ncstat /= nf90_noerr) then
       print *, "FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting dimension ID"
    end if
    ncstat = nf90_inquire_dimension(ncid,dimid,len=dim_obs_p)
    if (ncstat /= nf90_noerr) then
       print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting number of observations from NetCDF file'
    end if
    
    IF (dim_obs_p <= 0) THEN
    
      allocate(obs_p(1))
      allocate(ivariance_obs_p(1))
      allocate(ocoord_p(2, 1))
      allocate(thisobs%id_obs_p(3,1))
      
      obs_p=0.0
      ivariance_obs_p=0.0
      ocoord_p=0.0
      thisobs%id_obs_p=0
      
      allocate(nod1_p(0),nod2_p(0),nod3_p(0),nl_p(0),lon_p(0),lat_p(0))
      
      print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - No obs'
      
    ELSE ! (i.e. dim_obs_p > 0)
    
      ! Allocate memory before reading from netCDF
      allocate(obs_p(dim_obs_p))
      allocate(nod1_p(dim_obs_p),nod2_p(dim_obs_p),nod3_p(dim_obs_p))
      allocate(nl_p(dim_obs_p))
      allocate(lon_p(dim_obs_p),lat_p(dim_obs_p))
      allocate(ivariance_obs_p(dim_obs_p))
      allocate(ocoord_p(2,dim_obs_p))
      
      print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Have obs: ', dim_obs_p
      
      obs_p = 2300.0
      nod1_p = 1167
      nod2_p = 183
      nod3_p = 1198
      nl_p = 5
      ocoord_p(1,:)= 84.260/180*pi
      ocoord_p(2,:)=-56.059/190*pi
      
!~       ! Reading observations
!~       ncstat = nf90_inq_varid(ncid,'DIC', id_obs)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_obs from netCDF'
!~       ncstat = nf90_get_var(ncid, id_obs, obs_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading obs from netCDF'
      
!~       ! Reading nodes
!~       ncstat = nf90_inq_varid(ncid,'NOD1_P', id_nod1)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_nod1 from netCDF'
!~       ncstat = nf90_get_var(ncid, id_nod1, nod1_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading nod1 from netCDF'
      
!~       ncstat = nf90_inq_varid(ncid,'NOD2_P', id_nod2)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_nod2 from netCDF'
!~       ncstat = nf90_get_var(ncid, id_nod2, nod2_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading nod2 from netCDF'
      
!~       ncstat = nf90_inq_varid(ncid,'NOD3_P', id_nod3)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_nod3 from netCDF'
!~       ncstat = nf90_get_var(ncid, id_nod3, nod3_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading nod3 from netCDF'
      
!~       ! Reading layer
!~       ncstat = nf90_inq_varid(ncid,'NZ1', id_nl)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_layer from netCDF'
!~       ncstat = nf90_get_var(ncid, id_nl, nl_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading layer from netCDF'
      
!~       ! Reading coordinates
!~       ncstat = nf90_inq_varid(ncid,'LON', id_lon)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_lon from netCDF'
!~       ncstat = nf90_get_var(ncid, id_lon, lon_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading lon from netCDF'
      
!~       ncstat = nf90_inq_varid(ncid,'LAT', id_lat)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error getting id_lat from netCDF'
!~       ncstat = nf90_get_var(ncid, id_lat, lat_p)
!~       if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error reading lat from netCDF'
      
!~       ocoord_p(1,:) = lon_p / 180.0 * PI
!~       ocoord_p(2,:) = lat_p / 180.0 * PI
      
      ! *** Initialize index vector of observed surface nodes ***
      ! This array has as many rows as required for the observation operator
      ! 1 if observations are at grid points; >1 if interpolation is required
      allocate(thisobs%id_obs_p(3,dim_obs_p))
      DO i = 1, dim_obs_p
        thisobs%id_obs_p(1,i) = (mesh_fesom% nl-1) * (nod1_p(i)-1) + nl_p(i) + offset(id%DIC)
        thisobs%id_obs_p(2,i) = (mesh_fesom% nl-1) * (nod2_p(i)-1) + nl_p(i) + offset(id%DIC)
        thisobs%id_obs_p(3,i) = (mesh_fesom% nl-1) * (nod3_p(i)-1) + nl_p(i) + offset(id%DIC)
      END DO
      
      ! *** Set constant observation error *** 
      IF (mype_filter == 18) &
          WRITE (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
          '--- Use global GLODAP DIC observation error of ', rms_obs_DIC_glodap, 'mg chl m-3'
      ! Set inverse observation error variance
      ivariance_obs_p = 1.0 / (rms_obs_DIC_glodap ** 2)
      WRITE (*,*) 'Pe-local inverse observation error variance: ', ivariance_obs_p

    ENDIF
    
    ! **************************************
    ! *** Gather full observation arrays ***
    ! **************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivariance_obs_p, ocoord_p, &
                            thisobs%ncoord, lradius_DIC_glodap, dim_obs)
                              
    IF (mype_filter == 0) &
        WRITE (*, '(a, 5x, a, i)') 'FESOM-PDAF', &
        '--- Full GLODAP DIC observations have been gathered; dim_obs is ', dim_obs

    ! *** Clean-up ***
    ncstat = nf90_close(ncid)
    if (ncstat /= nf90_noerr) then
       print *, 'FESOM-PDAF - obs_DIC_glodap_pdafomi - Error closing NetCDF file'
    end if
    
    deallocate(obs_p,ivariance_obs_p,ocoord_p)
    deallocate(nod1_p,nod2_p,nod3_p,nl_p,lon_p,lat_p) 
    ! this_obs%id_obs_p is deallocated in PDAFomi_obs_l.F90

  END SUBROUTINE init_dim_obs_DIC_glodap



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_DIC_glodap(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_gridavg, &
               PDAFomi_set_debug_flag

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state

! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

! For EN4 profile observations handled here, the observation
! operator has to average the values of 3 grid points.
! For this the observation operator OBS_OP_F_GRIDAVG is used.

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_obs_op_gridavg(thisobs, 3, state_p, ostate)
    END IF

  END SUBROUTINE obs_op_DIC_glodap



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_DIC_glodap(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l,&
                       PDAFomi_set_debug_flag

    ! Include localization radius and local coordinates
    USE mod_assim_pdaf, ONLY: coords_l, locweight, loctype

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

! *** OMI-Debug:
!~   IF (mype_filter==64 .AND. domain_p==776) THEN
!~     CALL PDAFomi_set_debug_flag(domain_p)
!~   ELSE
!~     CALL PDAFomi_set_debug_flag(0)
!~   ENDIF

! **********************************************
! *** Initialize local observation dimension ***
! **********************************************
    IF (thisobs%doassim == 1) THEN
       
       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_DIC_glodap, sradius_DIC_glodap, dim_obs_l)
    END IF

  END SUBROUTINE init_dim_obs_l_DIC_glodap

END MODULE obs_DIC_glodap_pdafomi
