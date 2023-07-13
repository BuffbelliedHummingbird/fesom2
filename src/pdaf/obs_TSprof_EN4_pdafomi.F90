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
MODULE obs_TSprof_EN4_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter          ! Rank of filter process
  USE PDAFomi, &
       ONLY: obs_f, obs_l         ! Declaration of observation data types
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_o_en4_t                      !< Whether to assimilate temperature profiles
  LOGICAL :: assim_o_en4_s                      !< Whether to assimilate salinity profiles

  ! Further variables specific for the EN4 profile observations
  CHARACTER(len=110) :: path_obs_prof  = ''      !< Path to profile observations
  CHARACTER(len=110) :: file_prof_prefix  = ''   !< file name prefix for profile observations 
  CHARACTER(len=110) :: file_prof_suffix  = '.nc'!< file name suffix for profile observations 
  CHARACTER(len=110) :: file_syntobs_prof = 'syntobs_en4.nc' !< File name for synthetic observations

  REAL    :: rms_obs_T         ! RMS error for temperature profiles
  REAL    :: rms_obs_S         ! RMS error for salinity profiles
  REAL    :: bias_obs_prof     ! profile observation bias

  REAL    :: lradius_prof      ! Localization radius in the ocean for profiles
  REAL    :: sradius_prof      ! Support radius for localization function

  REAL    :: prof_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)

  REAL, ALLOCATABLE :: mean_temp_p (:)      ! mean temperature for observation exclusion
  REAL, ALLOCATABLE :: loc_radius_prof(:)   ! Localization radius array for profiles


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
  SUBROUTINE init_dim_obs_prof(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs, PDAFomi_set_debug_flag
    USE mod_assim_pdaf, &
         ONLY: offset, twin_experiment, use_global_obs, id, &
               delt_obs_ocn, mesh_fesom, &
               local_range, srange
    USE mod_parallel_pdaf, &
         ONLY: MPI_SUM, MPIerr, COMM_filter, MPI_INTEGER
    USE g_parsup, &
         ONLY: myDim_nod2D
    USE g_clock, &
         ONLY: month, day_in_month, yearold, timenew
    USE g_rotate_grid, &
         ONLY: r2g

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in)    :: step      !< Current time step
    INTEGER, INTENT(inout) :: dim_obs   !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, iter_file, s, k, j          ! Counters
    INTEGER :: dim_obs_p                      ! Number of process local observations
    CHARACTER(len=4) :: mype_string           ! String for process rank
    CHARACTER(len=4) :: year_string           ! String for yearly profile data path
    INTEGER :: stat(100)                      ! Status for NetCDF functions
    INTEGER :: status                         ! Status flag for PDAF gather operation
    INTEGER, ALLOCATABLE :: n2d_temp(:,:), &  ! Grid point indices read from file
                            n2d_sal(:,:), &
                            nl1_temp(:), &
                            nl1_sal(:)
    REAL, ALLOCATABLE :: obs_p(:)             ! PE-local observed T/S profile values
    REAL, ALLOCATABLE :: ocoord_n2d_p(:,:)    ! PE-local coordinates of observed profiles
    REAL, ALLOCATABLE :: ivariance_obs_p(:)   ! PE-local coordinates of observed profiles
    REAL, ALLOCATABLE :: obs_depth_p(:)       ! PE-local observation depths read from file
    CHARACTER(len=100) :: prof_file = ''      ! Complete name of profile observation file without path
    INTEGER :: ncid_in, id_temp, id_sal                ! IDs 
    INTEGER :: id_time, id_nprof_temp, id_nprof_sal    ! IDs
    INTEGER :: id_nflagged_temp, id_nflagged_sal       ! IDs
    INTEGER :: id_nobs_temp, id_nobs_sal               ! IDs
    INTEGER :: id_depth_temp, id_depth_sal             ! IDs
    INTEGER :: id_coord_temp, id_coord_sal             ! IDs
    INTEGER :: id_offset_temp, id_offset_sal           ! IDs
    INTEGER :: id_node_temp, id_node_sal               ! IDs
    INTEGER :: id_layer_temp, id_layer_sal             ! IDs
    INTEGER :: startv, countv, startv2(2), countv2(2)  ! Vectors for file reading
    INTEGER :: day                          ! Day read from file
    INTEGER :: offset_temp, offset_sal      ! Offset of profiles for the chosen day
    INTEGER :: cnt_prof_temp, cnt_prof_sal  ! Number of assimilated profiles for temperature and salinity
    INTEGER :: cnt_pro_temp_g, cnt_pro_sal_g! Total number of profiles at one day
    INTEGER :: cnt_temp, cnt_sal            ! Number of assimilated observations for temperature and salinity
    INTEGER :: cnt_flagged_temp_p, cnt_flagged_sal_p  ! PE-local number of flagged profiles
    INTEGER :: cnt_flagged_temp, cnt_flagged_sal      ! Global number of flagged profiles
    INTEGER :: cnt_ex_T_diff_p              ! PE-local count of excluded points due to temperature difference
    INTEGER :: cnt_ex_T_diff                ! Global count of excluded points due to temperature difference
    REAL, ALLOCATABLE :: obs_g(:)           ! Global full observation vector (used in case of limited obs.)
    REAL, ALLOCATABLE :: ivar_obs_g(:)      ! Global full inverse variances (used in case of limited obs.)
    REAL, ALLOCATABLE :: ocoord_g(:,:)      ! Global full observation coordinates (used in case of limited obs.)
    
    REAL    :: test1
    REAL    :: test2
    INTEGER :: my_debug_id_nod2

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    ! Store whether to assimilate this observation type
    IF (assim_o_en4_t .OR. assim_o_en4_s) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Initialize flag for type of full observations
    thisobs%use_global_obs = use_global_obs

    ! set localization radius
    lradius_prof = local_range
    sradius_prof = srange
    
    IF (.NOT.ALLOCATED(loc_radius_prof)) ALLOCATE(loc_radius_prof(mydim_nod2d))
    loc_radius_prof(:) = lradius_prof


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Initialize complete file name
    WRITE(mype_string,'(i4.4)') mype_filter
    WRITE(year_string,'(i4.4)') yearold

    prof_file=TRIM(file_prof_prefix)//TRIM(mype_string)//TRIM(file_prof_suffix)

    ! Position to read from file (which day)
    iter_file = step / delt_obs_ocn
    
    ! Debugging message:
    IF (mype_filter == 0) THEN
       WRITE (*,'(a,5x,a,i2,a,i2,a,i4,a,f5.2,a,i)') &
            'FESOM-PDAF', 'Assimilate EN4 profile observations - OBS_TSPROF_EN4 at ', &
            day_in_month, '.', month, '.', yearold, ' ', timenew/3600.0,&
            ' h; read at day: ', iter_file
    END IF

    ! Initialize no. of temperature and salinity observations
    cnt_temp = 0
    cnt_sal = 0
    cnt_prof_temp = 0
    cnt_prof_sal = 0
    cnt_flagged_temp_p = 0
    cnt_flagged_sal_p = 0

    ! Read profile observation
    s = 1
    stat(s) = NF_OPEN(TRIM(path_obs_prof)//'/'//year_string//'/'//prof_file, NF_NOWRITE,ncid_in)
    s = s + 1
    
    IF (mype_filter == 0) &
         WRITE (*,'(a,5x,a)') 'FESOM-PDAF', '--- Read PE-local profile observations from file: ', &
                               TRIM(path_obs_prof)//'/'//year_string//'/'//prof_file

    ! Get IDs for variables
    stat(s) = NF_INQ_VARID(ncid_in, 'day', id_time)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nprof_temp', id_nprof_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nprof_sal', id_nprof_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nflagged_temp', id_nflagged_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nflagged_sal', id_nflagged_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nobs_temp', id_nobs_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nobs_sal', id_nobs_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'offset_temp', id_offset_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'offset_sal', id_offset_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'depth_temp', id_depth_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'depth_sal', id_depth_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'temperature', id_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'salinity', id_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'n2d_temp', id_node_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'n2d_sal', id_node_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nl1_temp', id_layer_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'nl1_sal', id_layer_sal)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'coordinate_temp', id_coord_temp)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'coordinate_sal', id_coord_sal)
  
    ! Read day and number of observations
    startv = iter_file ! (day in the year)
    countv = 1

    s = s + 1
    stat(s) = NF_GET_VARA_INT(ncid_in, id_time, startv, countv, day)

    IF (mype_filter == 0) &
         WRITE (*,'(a,5x,a,i3)') 'FESOM-PDAF', '--- Read PE-local profile observations on day:', day

    IF (assim_o_en4_t) THEN
       s = s + 1
       stat(s) = NF_GET_VARA_INT(ncid_in, id_nprof_temp, startv, countv, cnt_prof_temp)
       s = s + 1
       stat(s) = NF_GET_VARA_INT(ncid_in, id_nobs_temp, startv, countv, cnt_temp)
       s = s + 1
       stat(s) = NF_GET_VARA_INT(ncid_in, id_nflagged_temp, startv, countv, cnt_flagged_temp_p)
    END IF

    IF (assim_o_en4_s) THEN
       s = s + 1
       stat(s) = NF_GET_VARA_INT(ncid_in, id_nprof_sal, startv, countv, cnt_prof_sal)
       s = s + 1
       stat(s) = NF_GET_VARA_INT(ncid_in, id_nobs_sal, startv, countv, cnt_sal)
       s = s + 1
       stat(s) = NF_GET_VARA_INT(ncid_in, id_nflagged_sal, startv, countv, cnt_flagged_sal_p)
    END IF
  
    DO i = 1, s
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in opening and reading profile variable, no.', i
    END DO


! ***************************************************
! *** Count available observations and initialize ***
! *** PE-local index and coordinate arrays.       ***
! ***************************************************

    ! *** Count PE-local number of observations ***
    dim_obs_p = cnt_temp + cnt_sal

    haveobs: IF (dim_obs_p > 0) THEN

       ! Allocate PE-local observation arrays
       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(obs_depth_p(dim_obs_p))
       ALLOCATE(ocoord_n2d_p(2,dim_obs_p))
     
       ! *** Initialize index vector of observed nodes ***
       ! This array has a many rows as required for the observation operator
       ! 1 if observations are at grid points; >1 if interpolation is required
       ALLOCATE(thisobs%id_obs_p(3,dim_obs_p))

       ! *** Temperature ***
          
       havetemp: IF (cnt_temp > 0) THEN

          ALLOCATE(n2d_temp(3,cnt_temp))
          ALLOCATE(nl1_temp(  cnt_temp))

          startv = iter_file ! (day in the year)
          countv = 1
        
          ! Read offset, depth, temperature and node index
          ! offset_temp:      Daily offset of values in nc-file
          ! n2d_temp:         Index refering to model node
          ! nl1_temp:         Index refering to model layer
          ! obs_depth_p:      Depth of the model layer
          ! obs_p:            Quality-checked temperature averaged on model depth
          ! o_coord_n2d_p:    Observation coordinate
          
          s = 1
          stat(s) = NF_GET_VARA_INT(ncid_in, id_offset_temp, startv, countv, offset_temp)

          startv = offset_temp
          countv = cnt_temp

          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(ncid_in, id_depth_temp, startv, countv, obs_depth_p(1:cnt_temp))
          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(ncid_in, id_temp, startv, countv, obs_p(1:cnt_temp))
          s = s + 1
          stat(s) = NF_GET_VARA_INT(ncid_in, id_layer_temp, startv, countv, nl1_temp(1:cnt_temp))

          startv2(1) = 1
          countv2(1) = 3       
          startv2(2) = offset_temp
          countv2(2) = cnt_temp

          s = s + 1
          stat(s) = NF_GET_VARA_INT(ncid_in, id_node_temp, startv2, countv2, n2d_temp(:,1:cnt_temp))
          
          DO i = 1, cnt_temp
             DO k = 1, 3
                thisobs%id_obs_p(k,i) = (mesh_fesom% nl-1) * (n2d_temp(k,i)-1) + nl1_temp(i) + offset(id%temp)
             END DO
          END DO
        
          countv2(1) = 2

          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(ncid_in, id_coord_temp, startv2, countv2, ocoord_n2d_p(:,1:cnt_temp))

          DO i = 1, s
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in reading profile variable --- temperature, no.', i
          END DO
        
       END IF havetemp
     
        
       ! *** Salinity ***

       havesal: IF (cnt_sal > 0) THEN

          ALLOCATE(n2d_sal(3,cnt_sal))
          ALLOCATE(nl1_sal(  cnt_sal))

          ! Read offset, depth, salinity and node index

          startv = iter_file
          countv = 1

          s = 1
          stat(s) = NF_GET_VARA_INT(ncid_in, id_offset_sal, startv, countv, offset_sal)

          startv = offset_sal
          countv = cnt_sal

          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(ncid_in, id_depth_sal, startv, countv, obs_depth_p(cnt_temp+1:cnt_temp+cnt_sal))
          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(ncid_in, id_sal, startv, countv, obs_p(cnt_temp+1 : cnt_temp+cnt_sal))
          s = s + 1
          stat(s) = NF_GET_VARA_INT(ncid_in, id_layer_sal, startv, countv, nl1_sal(1:cnt_sal))

          startv2(1) = 1
          countv2(1) = 3       
          startv2(2) = offset_sal
          countv2(2) = cnt_sal

          s = s + 1
          stat(s) = NF_GET_VARA_INT(ncid_in, id_node_sal, startv2, countv2, n2d_sal(:,1:cnt_sal))

          ! *** plus offset here instead of in obs_op_f_pdaf
          DO i = cnt_temp + 1, cnt_temp + cnt_sal
             DO k = 1, 3
                thisobs%id_obs_p(k,i) = (mesh_fesom% nl-1) * (n2d_sal(k,i-cnt_temp)-1) + nl1_sal(i-cnt_temp) + offset(id%salt)
             END DO
          END DO
 
!~           ! Debugging:
!~           IF (mype_filter==19) THEN
!~           open(3, file = 'id_obs_p.dat')
!~           write(3,*) thisobs%id_obs_p
!~           close(3)
!~           END IF

          countv2(1) = 2

          s = s + 1
          stat(s) = NF_GET_VARA_DOUBLE(ncid_in, id_coord_sal, startv2, countv2, ocoord_n2d_p(:,cnt_temp+1:cnt_temp+cnt_sal))

          DO i = 1, s
             IF (stat(i) /= NF_NOERR) &
                  WRITE(*, *) 'NetCDF error in reading profile variable --- salinity, no.', i
          END DO

       END IF havesal
 

! ***************************************
! *** Define local observation errors ***
! ***************************************

       ! *** Use different observation errors for different variables:
       ! *** rms_obs_T: Temperature  --> profile
       ! *** rms_obs_S: Salinity     --> profile

       ! Define observation error
       ! Currently use constant value for all grid points

!~        IF (assim_o_en4_t .AND. mype_filter == 0) &
!~             WRITE (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
!~             '--- global observation error for profile T: ', rms_obs_T, ' degC'
!~        IF (assim_o_en4_s .AND. mype_filter == 0) &
!~             WRITE (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
!~             '--- global observation error for profile S: ', rms_obs_S, ' psu'

!~        ALLOCATE(ivariance_obs_p(dim_obs_p))

!~        ! For temperature
!~        DO i = 1, cnt_temp
!~           ivariance_obs_p(i) = 1.0 / rms_obs_T**2
!~        END DO

!~        ! For salinity
!~        DO i = 1+cnt_temp, cnt_temp+cnt_sal
!~           ivariance_obs_p(i) = 1.0 / rms_obs_S**2
!~        END DO

          IF (assim_o_en4_t .AND. mype_filter == 0) &
                WRITE (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
                '--- global observation minimum error for profile T: ', rms_obs_T, ' degC. ! NOTE, it is scaled by depth'
           IF (assim_o_en4_s .AND. mype_filter == 0) &
                WRITE (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
                '--- global observation minimum error for profile S: ', rms_obs_S, ' psu.  ! NOTE, it is scaled by depth'
    
           ALLOCATE(ivariance_obs_p(dim_obs_p))
    
           ! For temperature
           DO i = 1, cnt_temp
              !ivariance_obs_p(i) = 1.0 / rms_obs_T**2
              ! Xie and Zhu, Ocean Modelling, 2010, vertical decay
              ivariance_obs_p(i) = 1.0 / (rms_obs_T+0.45*exp(0.002*obs_depth_p(i)))**2
           END DO
    
           ! For salinity
           ! https://scielo.conicyt.cl/scielo.php?script=sci_arttext&pid=S0717-65382004000300025
           DO i = 1+cnt_temp, cnt_temp+cnt_sal
              !ivariance_obs_p(i) = 1.0 / rms_obs_S**2
              ivariance_obs_p(i) = 1.0 / (rms_obs_S+0.03*exp(0.008*obs_depth_p(i)))**2
           END DO
        
       ! *** Exclude observations if difference from ensemble mean is beyond limit PROF_EXCLUDE_DIFF ***

       cnt_ex_T_diff_p = 0
       IF (prof_exclude_diff > 0.0) THEN
          DO i = 1, cnt_temp
             IF( (((mean_temp_p(thisobs%id_obs_p(1,i)) + mean_temp_p(thisobs%id_obs_p(2,i)) &
                  + mean_temp_p(thisobs%id_obs_p(3,i))) / 3.0 ) - obs_p(i)) > prof_exclude_diff) THEN
                ivariance_obs_p(i) = 1.0 / 1.0E12
                cnt_ex_T_diff_p = cnt_ex_T_diff_p+1
             END IF
          END DO
       END IF
     
    ELSE haveobs  ! IF (dim_obs_p > 0)
       ! No valid observations
     
       ALLOCATE(ocoord_n2d_p(2,1), obs_depth_p(1), thisobs%id_obs_p(3,1), obs_p(1), &
            ivariance_obs_p(1))

       ocoord_n2d_p = 0.0 
       obs_depth_p = 0.0
       thisobs%id_obs_p = 0
       obs_p = 0.0 
       ivariance_obs_p = 0.0
       cnt_ex_T_diff_p = 0

    END IF haveobs ! IF (dim_obs_p > 0)
  
    stat(1) = NF_CLOSE(ncid_in)


! ****************************
! *** De-bias observations ***
! ****************************

    IF (bias_obs_prof/=0.0 .AND. mype_filter == 0) &
         WRITE (*, '(a, 5x, a, f12.3)') &
         'FESOM-PDAF', '--- Use global observation bias of ', bias_obs_prof

    obs_p = obs_p - bias_obs_prof


! *************************************************
! *** Gather global information on observations ***
! *************************************************

    ! Get total number of profiles for temperature and salinity
    CALL MPI_Allreduce(cnt_prof_temp, cnt_pro_temp_g, 1, MPI_INTEGER, MPI_SUM, &
         COMM_filter, MPIerr)
    CALL MPI_Allreduce(cnt_prof_sal, cnt_pro_sal_g, 1, MPI_INTEGER, MPI_SUM, &
         COMM_filter, MPIerr)

    ! Get global number of flagged profiles
    CALL MPI_Allreduce(cnt_flagged_temp_p, cnt_flagged_temp, 1, MPI_INTEGER, MPI_SUM, &
         COMM_filter, MPIerr)
    CALL MPI_Allreduce(cnt_flagged_sal_p, cnt_flagged_sal, 1, MPI_INTEGER, MPI_SUM, &
         COMM_filter, MPIerr)

    ! Write number of profiles
    IF (mype_filter == 0) THEN
       IF (assim_o_en4_t .AND. (.NOT.assim_o_en4_s)) THEN
          WRITE (*,'(a,5x,a,i3,x,a,i5)') 'FESOM-PDAF', '--- Day', iter_file, &
               'number of T profiles', cnt_pro_temp_g
       ELSEIF ((.NOT.assim_o_en4_t) .AND. assim_o_en4_s) THEN
          WRITE (*,'(a,5x,a,i3,x,a,i5)') 'FESOM-PDAF', '--- Day', iter_file, &
               'number of S profiles', cnt_pro_sal_g
       ELSEIF (assim_o_en4_t .AND. assim_o_en4_s) THEN
          WRITE (*,'(a,5x,a,i3,x,a,i5,x,a,i5)') 'FESOM-PDAF', '--- Day', iter_file, &
               'number of profiles: T', cnt_pro_temp_g,'S',cnt_pro_sal_g
       ENDIF
       WRITE (*,'(a,5x,a,2i7)') 'FESOM-PDAF', &
            '--- Observations flagged at external nodes (T, S)', cnt_flagged_temp, cnt_flagged_sal
    END IF

    ! Get total number of observations if difference from ensemble mean is beyond limit PROFT_EXCLUDE_DIFF
    CALL MPI_Allreduce(cnt_ex_T_diff_p, cnt_ex_T_diff, 1, MPI_INTEGER, MPI_SUM, &
         COMM_filter, MPIerr)

    IF (mype_filter == 0) &
         WRITE (*,'(a,5x,a,f6.2,a,i7)') 'FESOM-PDAF', &
         '--- Observations excluded due to difference >',prof_exclude_diff,'degC:', cnt_ex_T_diff


! ***************************************************************
! *** Gather full observation arrays (obs, ivariance, occord) ***
! ***************************************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivariance_obs_p, ocoord_n2d_p, &
         thisobs%ncoord, lradius_prof, dim_obs)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

    ! This only works with use_full_obs=.true.
    IF (twin_experiment) THEN
       CALL read_syn_obs(file_syntobs_prof, dim_obs, thisobs%obs_f, 0, 1-mype_filter)
    END IF


! ********************
! *** Finishing up ***
! ********************

    ! Clean up arrays
    DEALLOCATE(obs_depth_p)
    DEALLOCATE(obs_p, ocoord_n2d_p, ivariance_obs_p)

  END SUBROUTINE init_dim_obs_prof



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
  SUBROUTINE obs_op_prof(dim_p, dim_obs, state_p, ostate)

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

  END SUBROUTINE obs_op_prof



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
  SUBROUTINE init_dim_obs_l_prof(domain_p, step, dim_obs, dim_obs_l)

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
    
       lradius_prof = loc_radius_prof(domain_p)
       
       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_prof, sradius_prof, dim_obs_l)
    END IF

  END SUBROUTINE init_dim_obs_l_prof

END MODULE obs_TSprof_EN4_pdafomi
