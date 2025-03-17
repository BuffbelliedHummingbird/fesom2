!> PDAF-OMI template observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
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
MODULE obs_o2_argo_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter, writepe ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l         ! Declaration of observation data types
  USE mod_assim_pdaf, &
       ONLY: n_sweeps, obs_PP     ! Variables for coupled data assimilation
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_o_o2_argo   !< Whether to assimilate

  ! Further variables specific for the profile observations
  CHARACTER(len=110) :: path_obs_o2_argo  = ''      !< Path to profile observations
  CHARACTER(len=110) :: file_o2_argo

  REAL    :: rms_obs_o2_argo      ! Observation error
  REAL    :: bias_obs_o2_argo     ! Observation bias

  REAL    :: lradius_o2_argo      ! Localization radius
  REAL    :: sradius_o2_argo      ! Support radius for localization function

  REAL    :: o2_argo_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)

  REAL, ALLOCATABLE :: mean_O2_p(:)            ! ensemble mean for observation exclusion
  REAL, ALLOCATABLE :: loc_radius_o2_argo(:)   ! localization radius array
  
  REAL, ALLOCATABLE :: ivariance_obs_g(:)      ! global-earth inverse observation variances


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
  
  type(obs_PP) :: thisobs_PP
  type(obs_PP) :: thisobs_PP_f
  
  LOGICAL :: isPP = .false.                 ! T: for postprocessing of simulation output
                                            ! F: for assimilation

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
  SUBROUTINE init_dim_obs_o2_argo(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs, PDAFomi_set_debug_flag
    USE mod_assim_pdaf, &
         ONLY: offset, use_global_obs, id, &
               mesh_fesom, nlmax, &
               local_range, srange
    USE mod_parallel_pdaf, &
         ONLY: MPI_SUM, MPIerr, COMM_filter, MPI_INTEGER
    USE g_parsup, &
         ONLY: myDim_nod2D
    USE g_clock, &
         ONLY: month, day_in_month, yearnew, timenew
    USE o_param, &
         ONLY: pi
         
    USE netcdf

    IMPLICIT NONE

!~     INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in)    :: step      !< Current time step
    INTEGER, INTENT(inout) :: dim_obs   !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, s, k, j, e, e_reps, e_found ! Counters
    CHARACTER(len=2) :: mype_string           ! String for process rank
    CHARACTER(len=4) :: year_string           ! String for yearly observation data path
    CHARACTER(len=2) :: mon_string            ! String for daily observation files
    CHARACTER(len=2) :: day_string            ! String for daily observation files
    
    CHARACTER(len=300) :: filename
    LOGICAL            :: FileExists
    
    ! observation data containing repetitions in case of repeatedly sampled elements
    INTEGER :: dim_obs_p_reps                         ! Number of process local observations
    REAL, ALLOCATABLE :: obs_p_reps(:)                ! PE-local observation values
    REAL, ALLOCATABLE :: lon_p_reps(:),lat_p_reps(:)  ! PE-local observed coords
    INTEGER, ALLOCATABLE :: nod1_p_reps(:), &
                            nod2_p_reps(:), &
                            nod3_p_reps(:)            ! PE-local observed node/element indeces
    INTEGER, ALLOCATABLE :: elem_p_reps(:)
    INTEGER, ALLOCATABLE :: nod1_g_reps(:), &
                            nod2_g_reps(:), &
                            nod3_g_reps(:)            ! Global observed node/element indeces
    INTEGER, ALLOCATABLE :: elem_g_reps(:)
    INTEGER, ALLOCATABLE :: nl_p_reps(:)
    REAL, ALLOCATABLE    :: depth_p_reps(:)
    REAL, ALLOCATABLE    :: std_obs_p_reps(:)         ! PE-local array of observation errors
    LOGICAL :: is_unique
    
    ! unique observation data: for each element, sorted array holds sum of all samples
    REAL, ALLOCATABLE :: obs_p_sorted(:)
    REAL, ALLOCATABLE :: x_p(:), y_p(:), z_p(:)
    INTEGER, ALLOCATABLE :: nod1_p_sorted(:), &
                            nod2_p_sorted(:), &
                            nod3_p_sorted(:)
    INTEGER, ALLOCATABLE :: elem_p_sorted(:)
    INTEGER, ALLOCATABLE :: nod1_g_sorted(:), &
                            nod2_g_sorted(:), &
                            nod3_g_sorted(:)
    INTEGER, ALLOCATABLE :: elem_g_sorted(:)
    INTEGER, ALLOCATABLE :: nl_p_sorted(:)
    REAL, ALLOCATABLE    :: depth_p_sorted(:)
    REAL, ALLOCATABLE    :: std_obs_p_sorted(:)
    INTEGER, ALLOCATABLE :: numrep_p(:)       ! number of observations at each element
    
    ! observation data after averaging over samples for each element, trimmed arrays:
    INTEGER :: dim_obs_p                      ! Number of process local observations
    
    REAL, ALLOCATABLE :: obs_p(:)             ! PE-local observation values
    REAL, ALLOCATABLE :: std_obs_p(:)         ! PE-local observation errors
    REAL, ALLOCATABLE :: ocoord_p(:,:)        ! PE-local coordinates of observations
    REAL, ALLOCATABLE :: lon_p(:),lat_p(:)
    REAL, ALLOCATABLE :: ivariance_obs_p(:)   ! PE-local array of inverse observation errors
    INTEGER, ALLOCATABLE :: nod1_p(:), &
                            nod2_p(:), &
                            nod3_p(:)         ! Array of observation pe-local indeces on FESOM grid
    INTEGER, ALLOCATABLE :: elem_p(:)
    INTEGER, ALLOCATABLE :: nl_p(:)           ! Array of observation layer indeces
    
    INTEGER :: ncstat                         ! Status for NetCDF functions
    INTEGER :: ncid, dimid                    ! NetCDF IDs
    INTEGER :: id_obs, id_nod1, id_nod2, &
               id_nod3, id_nl, id_lon, &
               id_lat, id_depth, id_elem, &
               id_err
               
    INTEGER :: cnt_ex_dry_p, cnt_ex_dry       ! Number of excluded observations due to model topography ("dry nodes")
    INTEGER :: nzmin                          ! Number of wet vertical model layers at observation location considering model topography
    
    REAL :: deg2rad
    REAL :: rad2deg
    
    deg2rad = PI/180.0
    rad2deg = 180.0/PI

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! set localization radius, everywhere the same
    lradius_o2_argo = 1.0e6 ! 1000 km
    sradius_o2_argo = 1.0e6 ! 1000 km
    
    ! set limit factor for omitted observation due to high innovation
    ! omit observation if innovation larger than this factor times
    ! observation error (only active for >0)
    thisobs%inno_omit =  o2_argo_exclude_diff / rms_obs_O2_argo

! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Daily files: nitialize complete file name
    WRITE(mype_string,'(i2.2)') mype_filter
    WRITE(year_string,'(i4.4)') yearnew
    WRITE(mon_string, '(i2.2)') month
    WRITE(day_string, '(i2.2)') day_in_month
    
    file_o2_argo = TRIM(year_string//'-'//mon_string//'-'//day_string//'.'//mype_string//'.nc')
    
    ! Message:
    IF (step < 0) then
        ! ---
        ! Postprocessing
        ! ---
        isPP = .true.
        if (writepe) WRITE (*,'(a,5x,a,1x,i2,a1,i2,a1,i4,1x,f5.2,a,1x,a)') &
            'FESOM-PDAF', 'Postprocessing of ARGO O2 observations at', &
            day_in_month, '.', month, '.', yearnew, timenew/3600.0,&
            'h; read from file:', file_o2_argo
        thisobs%use_global_obs = 1
        thisobs%doassim = 1
    ELSE
       ! ---
       ! Assimilation
       ! ---
       if (writepe) WRITE (*,'(a,5x,a,1x,i2,a1,i2,a1,i4,1x,f5.2,a,1x,a)') &
            'FESOM-PDAF', 'Assimilate ARGO O2 observations at', &
            day_in_month, '.', month, '.', yearnew, timenew/3600.0,&
            'h; read from file:', file_o2_argo
       ! Store whether to use global observations
       thisobs%use_global_obs = use_global_obs
       ! Store whether to assimilate this observation type
       IF (assim_o_o2_argo) thisobs%doassim = 1
    END IF
    
    ! complete filename
    filename=TRIM(path_obs_o2_argo)//year_string//'/dist72/'//TRIM(file_o2_argo)
    INQUIRE(file=TRIM(filename), exist=FileExists)
    
    IF (FileExists) THEN
       ! Open the pe-local NetCDF file for that day
       ncstat = nf90_open(TRIM(filename), nf90_nowrite, ncid)
       if (ncstat /= nf90_noerr) then
        print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error opening NetCDF file'
       end if
       
       ! Get the number of observations dim_obs_p_reps
       ncstat = nf90_inq_dimid(ncid, "INDEX", dimid)
       if (ncstat /= nf90_noerr) then
          print *, "FESOM-PDAF - obs_o2_argo_pdafomi - Error getting dimension ID"
       end if
       ncstat = nf90_inquire_dimension(ncid,dimid,len=dim_obs_p_reps)
       if (ncstat /= nf90_noerr) then
          print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting number of observations from NetCDF file'
       end if
    ELSE
       dim_obs_p_reps=0
    ENDIF ! FileExists
    
    cnt_ex_dry_p = 0
    cnt_ex_dry   = 0
    
    IF (dim_obs_p_reps <= 0) THEN
    
      allocate(obs_p(1))
      allocate(ivariance_obs_p(1))
      allocate(ocoord_p(2, 1))
      allocate(thisobs%id_obs_p(3,1))
      allocate(std_obs_p(1))
      
      obs_p=0.0
      ivariance_obs_p=0.0
      ocoord_p=0.0
      thisobs%id_obs_p=0
      
      allocate(lon_p(0),lat_p(0))
      
      if (isPP) then
        allocate(thisobs_PP%nod1_g(0),thisobs_PP%nod2_g(0),thisobs_PP%nod3_g(0))
        allocate(thisobs_PP%elem_g(0))
        allocate(thisobs_PP%depth(0),thisobs_PP%nz(0))
        allocate(thisobs_PP%isExclObs(0))
        allocate(thisobs_PP%lon(0),thisobs_PP%lat(0))
        allocate(thisobs_PP%numrep(0))
      endif
            
    ELSE ! (i.e. dim_obs_p > 0)
    
      ! Allocate memory before reading from netCDF
      allocate(obs_p_reps(dim_obs_p_reps))
      allocate(nod1_p_reps(dim_obs_p_reps),nod2_p_reps(dim_obs_p_reps),nod3_p_reps(dim_obs_p_reps))
      allocate(nl_p_reps(dim_obs_p_reps))
      allocate(lon_p_reps(dim_obs_p_reps),lat_p_reps(dim_obs_p_reps))
      allocate(elem_p_reps(dim_obs_p_reps))
      allocate(std_obs_p_reps(dim_obs_p_reps))
      
      if (isPP) then
         allocate(nod1_g_reps(dim_obs_p_reps),nod2_g_reps(dim_obs_p_reps),nod3_g_reps(dim_obs_p_reps))
         allocate(elem_g_reps(dim_obs_p_reps))
         allocate(depth_p_reps(dim_obs_p_reps))
      endif
            
      ! Reading observations
      ncstat = nf90_inq_varid(ncid,'OBS', id_obs)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_obs from netCDF'
      ncstat = nf90_get_var(ncid, id_obs, obs_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading obs from netCDF'
      
      ! Reading observation error
      ncstat = nf90_inq_varid(ncid,'ERR', id_err)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_err from netCDF'
      ncstat = nf90_get_var(ncid, id_err, std_obs_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading err from netCDF'
   
      ! Reading nodes
      ncstat = nf90_inq_varid(ncid,'NOD1_P', id_nod1)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_nod1 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod1, nod1_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading nod1 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD2_P', id_nod2)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_nod2 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod2, nod2_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading nod2 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD3_P', id_nod3)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_nod3 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod3, nod3_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading nod3 from netCDF'
      
      ncstat = nf90_inq_varid(ncid,'ELEM_P', id_elem)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_elem from netCDF'
      ncstat = nf90_get_var(ncid, id_elem, elem_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading elem from netCDF'
   
      ! Reading layer
      ncstat = nf90_inq_varid(ncid,'NZ1', id_nl)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_layer from netCDF'
      ncstat = nf90_get_var(ncid, id_nl, nl_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading layer from netCDF'
   
      ! Reading coordinates
      ncstat = nf90_inq_varid(ncid,'LON', id_lon)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_lon from netCDF'
      ncstat = nf90_get_var(ncid, id_lon, lon_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading lon from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'LAT', id_lat)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_lat from netCDF'
      ncstat = nf90_get_var(ncid, id_lat, lat_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading lat from netCDF'
      
      if (isPP) then
      
      ! Reading nodes on globe
      ncstat = nf90_inq_varid(ncid,'NOD1_G', id_nod1)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_nod1 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod1, nod1_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading nod1 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD2_G', id_nod2)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_nod2 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod2, nod2_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading nod2 from netCDF'
   
      ncstat = nf90_inq_varid(ncid,'NOD3_G', id_nod3)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_nod3 from netCDF'
      ncstat = nf90_get_var(ncid, id_nod3, nod3_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading nod3 from netCDF'
      
      ! Reading mesh element
      ncstat = nf90_inq_varid(ncid,'ELEM_G', id_elem)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_elem from netCDF'
      ncstat = nf90_get_var(ncid, id_elem, elem_g_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading elem from netCDF'
      
      ! Reading depth
      ncstat = nf90_inq_varid(ncid,'DEPTH', id_depth)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error getting id_depth from netCDF'
      ncstat = nf90_get_var(ncid, id_depth, depth_p_reps)
      if (ncstat /= nf90_noerr) print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error reading depth from netCDF'
      endif ! isPP
      
      ! **********************
      ! * Finding duplicates *
      ! **********************
      
      dim_obs_p = 0 ! Number of observations without duplicates
      e_found = 0
      
      ALLOCATE(obs_p_sorted(dim_obs_p_reps))
      ALLOCATE(elem_p_sorted(dim_obs_p_reps))
      ALLOCATE(nl_p_sorted(dim_obs_p_reps))
      ALLOCATE(nod1_p_sorted(dim_obs_p_reps),nod2_p_sorted(dim_obs_p_reps),nod3_p_sorted(dim_obs_p_reps))
      ALLOCATE(x_p(dim_obs_p_reps),y_p(dim_obs_p_reps),z_p(dim_obs_p_reps))
      ALLOCATE(std_obs_p_sorted(dim_obs_p_reps))
      ALLOCATE(numrep_p(dim_obs_p_reps))
      
      obs_p_sorted(:)     = 0
      std_obs_p_sorted(:) = 0
      elem_p_sorted(:) = 0
      nl_p_sorted(:)   = 0
      nod1_p_sorted(:) = 0
      nod2_p_sorted(:) = 0
      nod3_p_sorted(:) = 0
      x_p(:) = 0
      y_p(:) = 0
      z_p(:) = 0
      numrep_p(:) = 0
      
      if (isPP) then
         ALLOCATE(elem_g_sorted(dim_obs_p_reps))
         ALLOCATE(nod1_g_sorted(dim_obs_p_reps),nod2_g_sorted(dim_obs_p_reps),nod3_g_sorted(dim_obs_p_reps))
         ALLOCATE(depth_p_sorted(dim_obs_p_reps))
         elem_g_sorted(:)=0.0
         nod1_g_sorted(:) = 0
         nod2_g_sorted(:) = 0
         nod3_g_sorted(:) = 0
         depth_p_sorted(:) = 0
      endif ! isPP
      
      ! go through repreated observations one-by-one:
      do_elem_p_reps: DO e_reps=1, dim_obs_p_reps
         
         IF (obs_p_reps(e_reps)>=10.0) THEN
         ! O2 observations apparently spurious/small: exclusion
         ! else, continue here:
         
            is_unique = .true.
            
            ! searching for previous occurence of elem(e_reps) and nl(e_reps) within sorted arrays:
            do_elem_p: DO e=1,dim_obs_p
               
               IF ((elem_p_reps(e_reps)==elem_p_sorted(e)) .and. &
                   (nl_p_reps  (e_reps)==nl_p_sorted(e))) THEN
                  is_unique = .false.
                  e_found   = e
                  EXIT ! --> found previous occurence at e; stop searching.
               ENDIF
            ENDDO do_elem_p
            
            IF (.not. is_unique) THEN
            ! --> already have observations of elem(e_reps) at e_found:
                
                ! number of observations at e:
                numrep_p(e_found) = numrep_p(e_found) + 1
                
                ! adding to average of observations at e:
                obs_p_sorted(e_found)     = obs_p_sorted(e_found)     + obs_p_reps(e_reps)
                std_obs_p_sorted(e_found) = std_obs_p_sorted(e_found) + std_obs_p_reps(e_reps)
                
                ! adding to mean coordinate,
                ! transforming to cartesian coordinates:
                x_p(e_found) = x_p(e_found) &
                            + COS(deg2rad*lat_p_reps(e_reps)) * COS(deg2rad*lon_p_reps(e_reps))
                y_p(e_found) = y_p(e_found) &
                            + COS(deg2rad*lat_p_reps(e_reps)) * SIN(deg2rad*lon_p_reps(e_reps))
                z_p(e_found) = z_p(e_found) &
                            + SIN(deg2rad*lat_p_reps(e_reps))
                            
                ! adding to mean depth:
                IF (isPP) depth_p_sorted(e_found) = depth_p_sorted(e_found) + depth_p_reps(e_reps)
               
            ELSE
            ! --> is the first or only observation at elem(e_reps)
               
               dim_obs_p = dim_obs_p + 1
               numrep_p(dim_obs_p) = 1
               
               elem_p_sorted   (dim_obs_p) = elem_p_reps   (e_reps)
               nl_p_sorted     (dim_obs_p) = nl_p_reps     (e_reps)
               nod1_p_sorted   (dim_obs_p) = nod1_p_reps   (e_reps)
               nod2_p_sorted   (dim_obs_p) = nod2_p_reps   (e_reps)
               nod3_p_sorted   (dim_obs_p) = nod3_p_reps   (e_reps)
               obs_p_sorted    (dim_obs_p) = obs_p_reps    (e_reps)
               std_obs_p_sorted(dim_obs_p) = std_obs_p_reps(e_reps)
               
               if (isPP) then
                  elem_g_sorted (dim_obs_p) = elem_g_reps(e_reps)
                  nod1_g_sorted (dim_obs_p) = nod1_g_reps(e_reps)
                  nod2_g_sorted (dim_obs_p) = nod2_g_reps(e_reps)
                  nod3_g_sorted (dim_obs_p) = nod3_g_reps(e_reps)
                  depth_p_sorted(dim_obs_p) = depth_p_reps(e_reps)
               endif ! isPP

               ! transforming to cartesian coordinates:
               x_p(dim_obs_p) = COS(deg2rad*lat_p_reps(e_reps)) * COS(deg2rad*lon_p_reps(e_reps))
               y_p(dim_obs_p) = COS(deg2rad*lat_p_reps(e_reps)) * SIN(deg2rad*lon_p_reps(e_reps))
               z_p(dim_obs_p) = SIN(deg2rad*lat_p_reps(e_reps))

            ENDIF ! is_unique
         ENDIF ! obs >= 10.0
      ENDDO do_elem_p_reps
      
      ! Averaging repetitive samples
      allocate(obs_p(dim_obs_p))
      allocate(std_obs_p(dim_obs_p))
      allocate(lon_p(dim_obs_p),lat_p(dim_obs_p))
      
      DO e=1,dim_obs_p
        
        ! observation mean
        obs_p(e)     = obs_p_sorted(e)    /REAL(numrep_p(e))
        std_obs_p(e) = std_obs_p_sorted(e)/REAL(numrep_p(e))
        
        ! mean location (in cartesian coordinates)
        x_p(e) = x_p(e)/REAL(numrep_p(e))
        y_p(e) = y_p(e)/REAL(numrep_p(e))
        z_p(e) = z_p(e)/REAL(numrep_p(e))
        
        ! mean depth
        if (isPP) then
           depth_p_sorted(e) = depth_p_sorted(e)/REAL(numrep_p(e))
        endif
        
        ! transforming to spherical coordinates
        lat_p(e) = ATAN2( z_p(e), SQRT(x_p(e)*x_p(e) + y_p(e)*y_p(e)) )
        lon_p(e) = ATAN2( y_p(e), x_p(e))
        
      ENDDO
            
      allocate(ivariance_obs_p(dim_obs_p))
      allocate(ocoord_p(2,dim_obs_p))
      
      if (isPP) then
          allocate(thisobs_PP%nod1_g(dim_obs_p),thisobs_PP%nod2_g(dim_obs_p),thisobs_PP%nod3_g(dim_obs_p))
          allocate(thisobs_PP%elem_g(dim_obs_p))
          allocate(thisobs_PP%isExclObs(dim_obs_p))
          allocate(thisobs_PP%depth(dim_obs_p),thisobs_PP%nz(dim_obs_p))
          allocate(thisobs_PP%lon(dim_obs_p),thisobs_PP%lat(dim_obs_p))
          allocate(thisobs_PP%numrep(dim_obs_p))
          
          thisobs_PP%nod1_g = nod1_g_sorted(:dim_obs_p)
          thisobs_PP%nod2_g = nod2_g_sorted(:dim_obs_p)
          thisobs_PP%nod3_g = nod3_g_sorted(:dim_obs_p)
          
          thisobs_PP%elem_g = elem_g_sorted(:dim_obs_p)
          
          thisobs_PP%lon = lon_p(:dim_obs_p) / pi * 180.0
          thisobs_PP%lat = lat_p(:dim_obs_p) / pi * 180.0
          
          thisobs_PP%depth = depth_p_sorted(:dim_obs_p)
          thisobs_PP%nz    = nl_p_sorted(:dim_obs_p)
          
          thisobs_PP%numrep = numrep_p(:dim_obs_p)
          
      endif
   
      ocoord_p(1,:) = lon_p
      ocoord_p(2,:) = lat_p
      
      ! Set inverse observation error variance
      DO i = 1, dim_obs_p
         ivariance_obs_p(i) = 1.0 / (std_obs_p(i) * std_obs_p(i))
      ENDDO
      if (isPP) then
         thisobs_PP%isExclObs = 0
      endif
      
      ! *** Initialize index vector of observed surface nodes ***
      ! This array has as many rows as required for the observation operator
      ! 1 if observations are at grid points; >1 if interpolation is required
      allocate(thisobs%id_obs_p(3,dim_obs_p))
      DO i = 1, dim_obs_p
        thisobs%id_obs_p(1,i) = (nlmax) * (nod1_p_sorted(i)-1) + nl_p_sorted(i) + offset(id%O2)
        thisobs%id_obs_p(2,i) = (nlmax) * (nod2_p_sorted(i)-1) + nl_p_sorted(i) + offset(id%O2)
        thisobs%id_obs_p(3,i) = (nlmax) * (nod3_p_sorted(i)-1) + nl_p_sorted(i) + offset(id%O2)
        
      ! *** exclude observations at dry nodes ***
      ! number of layers at nodes considering bottom topography: mesh_fesom% nlevels_nod2D
        nzmin = MIN ( mesh_fesom% nlevels_nod2D( nod1_p_sorted(i) ), &
                      mesh_fesom% nlevels_nod2D( nod2_p_sorted(i) ), &
                      mesh_fesom% nlevels_nod2D( nod3_p_sorted(i) ))
                      
        IF (nl_p_sorted(i) >= nzmin) THEN
           ivariance_obs_p(i) = 1e-12
           cnt_ex_dry_p = cnt_ex_dry_p + 1
           if (isPP) thisobs_PP%isExclObs(i) = 1
        ENDIF
       
      END DO
      
    ENDIF ! IF (dim_obs_p = 0) ELSEIF (dim_obs_p > 0)
    
    CALL MPI_Allreduce(cnt_ex_dry_p, cnt_ex_dry, 1, MPI_INTEGER, MPI_SUM, &
            COMM_filter, MPIerr)

    IF (mype_filter == 0) &
            WRITE (*,'(a,5x,a,2x,i7)') 'REcoM-PDAF', &
            '--- O2 ARGO observations excluded due to topography: ', cnt_ex_dry
    
    ! **************************************
    ! *** Gather full observation arrays ***
    ! **************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivariance_obs_p, ocoord_p, &
                            thisobs%ncoord, lradius_o2_argo, dim_obs)
                              
    IF (mype_filter == 0) &
        WRITE (*, '(a, 5x, a, i)') 'FESOM-PDAF', &
        '--- Full ARGO O2 observations have been gathered; dim_obs is ', dim_obs
    
    ! Global inverse variance array (thisobs%ivar_obs_f)
    ! has been gathered, but, in case of coupled DA / "double-sweep",
    ! it will be reset during each sweep. Thus, save a copy:
    if (n_sweeps>1) then
       if (allocated(ivariance_obs_g)) deallocate(ivariance_obs_g)
       allocate(ivariance_obs_g(dim_obs))
       ivariance_obs_g = thisobs%ivar_obs_f
    end if


    ! *** Clean-up ***
    if (FileExists) THEN
       ncstat = nf90_close(ncid)
       if (ncstat /= nf90_noerr) then
          print *, 'FESOM-PDAF - obs_o2_argo_pdafomi - Error closing NetCDF file'
       end if
    ENDIF
    
    deallocate(ivariance_obs_p,ocoord_p,obs_p)
    deallocate(lat_p,lon_p)
    deallocate(std_obs_p)
    IF (dim_obs_p_reps>0) THEN
      deallocate(nod1_p_sorted,nod2_p_sorted,nod3_p_sorted,obs_p_sorted,elem_p_sorted,nl_p_sorted)
      deallocate(std_obs_p_sorted)
      deallocate(nod1_p_reps,nod2_p_reps,nod3_p_reps,lon_p_reps,lat_p_reps,obs_p_reps,elem_p_reps,nl_p_reps)
      deallocate(std_obs_p_reps)
      deallocate(x_p,y_p,z_p,numrep_p)
    ENDIF
    ! this_obs%id_obs_p is deallocated in PDAFomi_obs_l.F90

  END SUBROUTINE init_dim_obs_o2_argo



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
  SUBROUTINE obs_op_o2_argo(dim_p, dim_obs, state_p, ostate)

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

! For profile observations handled here, the observation
! operator has to average the values of 3 grid points.
! For this the observation operator OBS_OP_F_GRIDAVG is used.

    IF (thisobs%doassim == 1) THEN
       CALL PDAFomi_obs_op_gridavg(thisobs, 3, state_p, ostate)
    END IF

  END SUBROUTINE obs_op_o2_argo



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
  SUBROUTINE init_dim_obs_l_o2_argo(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l,&
                       PDAFomi_set_debug_flag

    ! Include localization radius and local coordinates
    USE mod_assim_pdaf, ONLY: coords_l, locweight, loctype
    
    ! Number of domains per sweep:
    USE g_parsup, ONLY: myDim_nod2D

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain: containing repititive sweeps
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector

! *** OMI-Debug:
!~   IF (mype_filter==0 .AND. domain_p==5) THEN
!~     CALL PDAFomi_set_debug_flag(domain_p)
!~   ELSE
!~     CALL PDAFomi_set_debug_flag(0)
!~   ENDIF


    IF (thisobs%doassim == 1) THEN
    
       ! ************************************************************
       ! *** Adapt observation error for coupled DA (double loop) ***
       ! ************************************************************
    
       if (n_sweeps>1) then
       
          ! Physics observations sweep.
          ! Set inverse observation error to small value
          if (domain_p==1) then

             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                   '--- PHY sweep: set ivar_obs_f for ARGO O2 to 1.0e-12'
             thisobs%ivar_obs_f = 1.0e-12
             
          ! BGC observations sweep.
          elseif (domain_p==myDim_nod2D+1) then
             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                  '--- BIO sweep: set ivar_obs_f for ARGO O2 to normal'
             thisobs%ivar_obs_f(:) = ivariance_obs_g
          end if
       end if ! n_sweeps

       ! **********************************************
       ! *** Initialize local observation dimension ***
       ! **********************************************
       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_o2_argo, sradius_o2_argo, dim_obs_l)
    
     END IF ! thisobs%doassim
  END SUBROUTINE init_dim_obs_l_o2_argo

END MODULE obs_o2_argo_pdafomi
