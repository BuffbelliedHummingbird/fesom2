!> PDAF-OMI observation module 
!!
!! This module handles operations for one data type (called 'module-type' below).
!! 
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
MODULE obs_chl_cci_pdafomi

  USE mod_parallel_pdaf, &
       ONLY: mype_filter, writepe
  USE PDAFomi, &
       ONLY: obs_f, obs_l, & ! Declaration of observation data types
       PDAFomi_set_debug_flag
  USE mod_assim_pdaf, &
       ONLY: n_sweeps             ! Variables for coupled data assimilation

  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_o_chl_cci      ! Whether to assimilate data

  ! Further variables specific for the OC CCO chlorophyll observations
  CHARACTER(len=100) :: path_obs_chl_cci  = ''      ! Path to observations
  CHARACTER(len=110) :: file_chl_cci_prefix  = ''   ! file name prefix for observations 
  CHARACTER(len=110) :: file_chl_cci_suffix  = '.nc'! file name suffix for observations 

  REAL    :: rms_obs_chl_cci      ! Observation error standard deviation
  REAL    :: bias_obs_chl_cci     ! Observation bias

  REAL    :: lradius_chl_cci      ! Localization radius in the ocean
  REAL    :: sradius_chl_cci      ! Support radius for localization function

  LOGICAL :: chl_cci_exclude_ice  ! Whether to exclude observations at grid points with ice
  REAL    :: chl_cci_exclude_diff ! Limit difference beyond which observations are excluded (0.0 to deactivate)
  LOGICAL :: chl_cci_fixed_rmse   ! Whether to use a fixed RMS error or the error provided with the data

  REAL, ALLOCATABLE :: mean_chl_cci_p (:)    ! Mean value for observation exclusion
  REAL, ALLOCATABLE :: loc_radius_chl_cci(:) ! Localization radius array
  REAL, ALLOCATABLE :: ivariance_obs_g(:)    ! Global-earth inverse observation variances
  
  LOGICAL :: chl_logarithmic = .true.        ! Whether to apply log-transformation
  
  CHARACTER(len=300) :: path_bias
  CHARACTER(len=100) :: file_bias_prefix
  
  REAL, PARAMETER :: base10toE = 2.302585092994 ! = 1.0/(LOG10(EXP(1.0)))

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
  SUBROUTINE init_dim_obs_chl_cci(step, dim_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs
    USE mod_assim_pdaf, &
         ONLY: offset, twin_experiment, use_global_obs, id, mesh_fesom, &
         local_range, srange, nlmax
    USE mod_assim_pdaf, &
         ONLY: delt_obs_ocn
    USE mod_parallel_pdaf, &
         ONLY: MPI_SUM, MPIerr, COMM_filter, MPI_INTEGER
    USE g_parsup, &
         ONLY: mydim_nod2d, myList_nod2d
    USE g_rotate_grid, &
         ONLY: r2g
    USE g_clock, &
         ONLY: month, day_in_month, yearnew, timenew, daynew
    USE obs_sst_pdafomi, &
         ONLY: mean_ice_p

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs    !< Dimension of full observation vector

! *** Local variables ***
    INTEGER :: i, iter_file, i_obs       ! Counters
    INTEGER :: dim_obs_p                 ! number of PE-local observations
    INTEGER :: fileid                    ! ID for NetCDF file
    INTEGER :: id_state, id_std          ! ID for state
    INTEGER :: stat(100)                 ! Status for NetCDF functions
    INTEGER :: startv(2),countv(2)       ! Vectors for reading fields
    CHARACTER(len=5)   :: mype_string    ! String for process rank
    CHARACTER(len=100) :: obs_file = ''     ! Complete name of observation file without path
    REAL(4), ALLOCATABLE :: all_obs_p(:)    ! PE-local complete observation field read from file
    REAL(4), ALLOCATABLE :: all_std_p(:)    ! PE-local complete observation error field read from file
    REAL, ALLOCATABLE :: obs_error_p(:)     ! PE-local observation error
    REAL, ALLOCATABLE :: obs_p(:)           ! PE-local observed field
    REAL, ALLOCATABLE :: ivariance_obs_p(:) ! PE-local inverse observation error variance
    REAL, ALLOCATABLE :: ocoord_n2d_p(:,:)  ! PE-local coordinates of observations
    INTEGER :: status                       ! Status flag for PDAF gather operation
    INTEGER :: cnt_ex_ice_p, cnt_ex_diff_p  ! PE-local counts of excluded points due to ice and difference
    INTEGER :: cnt_ex_ice, cnt_ex_diff      ! Global counts of excluded points due to ice and difference
    REAL, ALLOCATABLE :: obs_g(:)           ! Global full observation vector (used in case of limited obs.)
    REAL, ALLOCATABLE :: ivar_obs_g(:)      ! Global full inverse variances (used in case of limited obs.)
    REAL, ALLOCATABLE :: ocoord_g(:,:)      ! Global full observation coordinates (used in case of limited obs.)
    INTEGER, ALLOCATABLE :: obs_include_index(:)   ! Index of observed (not excluded!) surface nodes on process domain
    INTEGER :: dim_obs_f                      ! Global full observation number
    
    CHARACTER(len=300) :: bias_file = ''     ! Complete name of observation file without path
    INTEGER :: fileidbias                    ! ID for NetCDF file
    INTEGER :: id_bias                       ! ID for state
    REAL(4), ALLOCATABLE :: all_bias_p(:)    ! PE-local complete observation field read from file

! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    ! Store whether to assimilate this observation type
    IF (assim_o_chl_cci) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 2   ! 2=Geographic

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2

    ! Initialize flag for type of full observations
    if (mype_filter==0) write(*,*) 'obs_chl_cci_pdafomi: use_global_obs', use_global_obs
    thisobs%use_global_obs = use_global_obs

    ! set localization radius
    lradius_chl_cci = local_range
    sradius_chl_cci = srange
    if (allocated(loc_radius_chl_cci)) deallocate(loc_radius_chl_cci)
    allocate(loc_radius_chl_cci(mydim_nod2d))
    loc_radius_chl_cci(:) = lradius_chl_cci


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Initialize complete file name
    WRITE(mype_string,'(i5.5)') mype_filter

    obs_file =TRIM(file_chl_cci_prefix)//TRIM(mype_string)//TRIM(file_chl_cci_suffix)
    bias_file=TRIM(path_bias)//TRIM(file_bias_prefix)//TRIM(mype_string)//'.nc'

    ! Allocate array
    ALLOCATE(all_obs_p (myDim_nod2D))
    ALLOCATE(all_std_p (myDim_nod2D))
    ALLOCATE(all_bias_p(myDim_nod2D))

    ! Position to read from file
    iter_file = daynew
    
    ! Debugging message:
    IF (mype_filter==0) THEN
       WRITE (*,'(a,5x,a,i2,a,i2,a,i4,a,f5.2,a,i)') &
            'FESOM-PDAF', 'Assimilate OC-CCI observations - OBS_CHL_CCI_PDAFOMI at ', &
            day_in_month, '.', month, '.', yearnew, ' ', timenew/3600.0,&
            ' h; read at ', iter_file
    END IF

    ! Read observation
         
    stat(1) = NF_OPEN(TRIM(path_obs_chl_cci)//TRIM(obs_file), NF_NOWRITE, fileid)        
 
    ! *** Read state estimate ***

    stat(2) = NF_INQ_VARID(fileid, 'obs', id_state)

    startv(2) = iter_file
    countv(2) = 1
    startv(1) = 1
    countv(1) = myDim_nod2D 

    stat(3) = NF_GET_VARA_REAL (fileid, id_state, startv, countv, all_obs_p)
    
    ! *** Read standard deviation ***
    ! Root-mean-square-difference of log10-transformed chlorophyll-a concentration in seawater.
   
    stat(4) = NF_INQ_VARID(fileid, 'std', id_std)
   
    startv(2) = iter_file
    countv(2) = 1
    startv(1) = 1
    countv(1) = myDim_nod2D 
   
    stat(5) = NF_GET_VARA_REAL (fileid, id_std, startv, countv, all_std_p)

    ! *** close file  ***
    stat(6) = NF_CLOSE(fileid)
 
    ! check status flag
    DO i=1,6
       IF (stat(i).NE.NF_NOERR) WRITE(*,*) &
            'NetCDF error in reading full OC_CCI, no.',i, &
            ' file ',obs_file
       stat(i)=0
    END DO
    
! ***********************************
! *** Read PE-local bias FREE-OBS ***
! ***********************************

    ! Debugging message:
    IF ((mype_filter==0) .and. (iter_file==1)) THEN
       WRITE (*,'(a,5x,a,i2,a,i2,a,i4,a,f5.2,a,i)') &
            'FESOM-PDAF', 'De-Bias OC-CCI observations - OBS_CHL_CCI_PDAFOMI at ', &
            day_in_month, '.', month, '.', yearnew, ' ', timenew/3600.0,&
            ' h; read at ', iter_file
    END IF
    
    IF (bias_obs_chl_cci/=0.0) THEN
       
       ! Read bias
         
        stat(1) = NF_OPEN(TRIM(bias_file), NF_NOWRITE, fileidbias)        
	    
        ! *** Read state estimate ***
	    
        stat(2) = NF_INQ_VARID(fileidbias, 'bias', id_bias)
	    
        startv(2) = iter_file
        countv(2) = 1
        startv(1) = 1
        countv(1) = myDim_nod2D 
	    
        stat(3) = NF_GET_VARA_REAL (fileidbias, id_bias, startv, countv, all_bias_p)
	    
        ! *** close file  ***
        stat(4) = NF_CLOSE(fileidbias)
	    
        ! check status flag
        DO i=1,4
           IF (stat(i).NE.NF_NOERR) WRITE(*,*) &
                'NetCDF error in reading OC_CCI Bias, no.',i, &
                ' file ', TRIM(bias_file)
           stat(i)=0
        END DO
       
    ELSE
       all_bias_p(:)=0
    ENDIF


! ****************************
! *** Exclude observations ***
! ****************************

    ! *** Exclude observations if mean_ice is not zero ***

    exclude_ice: IF (chl_cci_exclude_ice) THEN

       cnt_ex_ice_p = 0
       DO i = 1, myDim_nod2D
          IF ((mean_ice_p(i) > 0.0) &
             .AND. &
              (abs(all_obs_p(i))<999.0)) &
          THEN
             all_obs_p(i) = 1.0e6
             cnt_ex_ice_p = cnt_ex_ice_p + 1
          END IF
       END DO
       
       ! *** Sum PE-local excluded nodes to global number of nodes excluded ***

       CALL MPI_Allreduce(cnt_ex_ice_p, cnt_ex_ice, 1, MPI_INTEGER, MPI_SUM, &
            COMM_filter, MPIerr)

       IF (mype_filter == 0) &
            WRITE (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
            '--- Observations excluded because of ice', cnt_ex_ice

       ! *** Set localization radius to zero for grid points with ice ***

       IF (mype_filter == 0) &
            WRITE (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
            '--- Set localization radius to zero for points with ice'
       DO i = 1, myDim_nod2D
          IF (mean_ice_p (i) > 0.0) THEN
             loc_radius_chl_cci(i) = 0.0
          END IF
       END DO

    END IF exclude_ice


    ! *** Exclude observations if difference from ensemble mean is beyond limit CHL_EXCLUDE_DIFF ***

    exclude_diff: IF (chl_cci_exclude_diff > 0.0) THEN

       cnt_ex_diff_p = 0
       DO i = 1, myDim_nod2D
          IF (ABS(mean_chl_cci_p(i) - all_obs_p(i)) > chl_cci_exclude_diff .AND. abs(all_obs_p(i))<=999.0) THEN
             all_obs_p(i) = 1.0e6
             cnt_ex_diff_p = cnt_ex_diff_p+1
          END IF
       END DO

       CALL MPI_Allreduce(cnt_ex_diff_p, cnt_ex_diff, 1, MPI_INTEGER, MPI_SUM, &
            COMM_filter, MPIerr)

       IF (mype_filter == 0) &
            WRITE (*,'(a,5x,a,f6.2,a,i7)') 'FESOM-PDAF', &
            '--- Chlorophyll CCI Observations excluded due to difference >',chl_cci_exclude_diff,'mg chl m-3:', cnt_ex_diff

    END IF exclude_diff
  

! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count PE-local number of observations ***
    dim_obs_p = 0
    DO i = 1, myDim_nod2d
       IF (ABS(all_obs_p(i)) < 999.0) dim_obs_p=dim_obs_p+1
    ENDDO

    haveobs: IF(dim_obs_p>0) THEN
    
		! *** Initialize index vector of observed surface nodes ***
		! This array has a many rows as required for the observation operator
		! 1 if observations are at grid points; >1 if interpolation is required
		! 1 if one model field is required; >1 if value is to be derived from multiple model fields
		
		ALLOCATE(thisobs%id_obs_p(2, dim_obs_p))
		ALLOCATE(obs_include_index(dim_obs_p))

		i_obs=0
		DO i = 1, myDim_nod2d
		   IF (ABS(all_obs_p(i)) < 999.0) THEN
			  i_obs = i_obs + 1
			  
			  ! index for state vector
			  ! row 1: DiaChl
			  thisobs%id_obs_p(1, i_obs) = &
			  (i-1) * (nlmax) + 1 + offset(id% DiaChl)
			  
			  ! row 2: PhyChl
			  thisobs%id_obs_p(2, i_obs) = &
			  (i-1) * (nlmax) + 1 + offset(id% PhyChl)
			  
			  ! index for all_obs_p and surface nod2d vector, respectively. 
			  obs_include_index(i_obs) = i
		   ENDIF
		ENDDO

		! *** Initialize PE-local vectors of observations ***
		ALLOCATE(obs_p(dim_obs_p))
		ALLOCATE(obs_error_p(dim_obs_p))
		
		DO i = 1, dim_obs_p
		   IF (chl_logarithmic) THEN
		                       !      ** log-transformation of chlorophyll **    ** de-bias observations **
		      obs_p(i)       = REAL(  LOG(all_obs_p(obs_include_index(i)))       + all_bias_p(obs_include_index(i)), 8)
		                       ! logarithmic error is provided
		      obs_error_p(i) = base10toE * REAL(all_std_p(obs_include_index(i)), 8)
		   ELSE
		                       ! chlorophyll is provided
		      obs_p(i)       = REAL(all_obs_p(obs_include_index(i)), 8)
		                      ! logarithmic error from file (all_std_p) serves as relative error here,
		                      ! which holds for relative errors << 1
		      obs_error_p(i) = REAL(all_std_p(obs_include_index(i)), 8) * REAL(all_obs_p(obs_include_index(i)), 8)

		   ENDIF
		ENDDO

		! *** Initialize coordinate arrays for PE-local observations
		ALLOCATE(ocoord_n2d_p(2, dim_obs_p))
		DO i = 1, dim_obs_p
		   ! Rotate to geographic coordinates and store
		   CALL r2g(ocoord_n2d_p(1, i), ocoord_n2d_p(2, i), &
				mesh_fesom%coord_nod2d(1, obs_include_index(i)), mesh_fesom%coord_nod2d(2, obs_include_index(i)))
		ENDDO


	! ****************************************************************
	! *** Define observation errors for process-local observations ***
	! ****************************************************************

		IF (chl_cci_fixed_rmse) THEN

		   ! *** Set constant observation error *** 
		   IF (mype_filter == 0) &
				WRITE (*, '(a, 5x, a, f12.3, a)') 'FESOM-PDAF', &
				'--- Use global OC-CCI observation error of ', rms_obs_chl_cci, 'mg chl m-3'

		   obs_error_p(:) = rms_obs_chl_cci
		ELSE
		
          ! *** Use variable error from file
          IF (mype_filter == 0) &
          		WRITE (*,'(a,5x,a,i7)') 'FESOM-PDAF', &
                '--- Use variable OC-CCI observation error from file'

		END IF

		! Set inverse observation error variance
		ALLOCATE(ivariance_obs_p(dim_obs_p))
		DO i = 1, dim_obs_p
		   ivariance_obs_p(i) = 1.0 / obs_error_p(i)**2
		END DO

		
	ELSE ! (i.e. dim_obs_p==0)
	
		ALLOCATE(obs_p(1))
		ALLOCATE(ivariance_obs_p(1))
		ALLOCATE(ocoord_n2d_p(2, 1))
		ALLOCATE(thisobs%id_obs_p(2,1))
		thisobs%id_obs_p = 0
		
		ALLOCATE(obs_include_index(1))
		ALLOCATE(obs_error_p(1))
		
	ENDIF haveobs

! *******************************************
! *** No global observations? - Fictional ***
! *******************************************

    IF (dim_obs_p==0) THEN
       obs_p=1.0
       ivariance_obs_p=1E-12
       ocoord_n2d_p(1,1)=1.57
       ocoord_n2d_p(2,1)=0.0
       thisobs%id_obs_p(1,1)=offset(id% PhyChl)+1
       thisobs%id_obs_p(2,1)=offset(id% DiaChl)+1
       dim_obs_p=1
    ENDIF

! **************************************
! *** Gather full observation arrays ***
! **************************************

    CALL PDAFomi_gather_obs(thisobs, dim_obs_p, obs_p, ivariance_obs_p, ocoord_n2d_p, &
         thisobs%ncoord, lradius_chl_cci, dim_obs)
         
    ! Global inverse variance array (thisobs%ivar_obs_f)
    ! has been gathered, but, in case of coupled DA / "double-sweep",
    ! it will be reset during each sweep. Thus, save a copy:
    if (n_sweeps>1) then
       if (allocated(ivariance_obs_g)) deallocate(ivariance_obs_g)
       allocate(ivariance_obs_g(dim_obs))
       ivariance_obs_g = thisobs%ivar_obs_f
    end if


! ********************
! *** Finishing up ***
! ********************

    ! Clean up arrays
    DEALLOCATE(all_obs_p, all_std_p, all_bias_p)
    DEALLOCATE(obs_error_p)
    DEALLOCATE(obs_p, ocoord_n2d_p, ivariance_obs_p)
    if (allocated(obs_include_index)) deallocate(obs_include_index)

  END SUBROUTINE init_dim_obs_chl_cci



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
  SUBROUTINE obs_op_chl_cci(dim_p, dim_obs, state_p, ostate)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obsstate

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs               !< Dimension of full observed state (all observed fields)
    REAL, INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL, INTENT(inout) :: ostate(dim_obs)       !< Full observed state
    
! *** Local ***
    REAL, ALLOCATABLE   :: ostate_p(:)           !< Pe-local observed state
    INTEGER :: i                                 !< Counter


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim == 1) THEN
       
       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF
       
       ! Initialize observed pe-local state vector by sum of rows 1 and 2 (DiaChl and PhyChl)
       DO i = 1, thisobs%dim_obs_p
          IF (chl_logarithmic) THEN
             ! log-transform
             ostate_p(i) = LOG(state_p(thisobs%id_obs_p(1,i)) + state_p(thisobs%id_obs_p(2,i)))
          ELSE
             ostate_p(i) = state_p(thisobs%id_obs_p(1,i)) + state_p(thisobs%id_obs_p(2,i))
          ENDIF
       END DO
       
       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate(thisobs, ostate_p, ostate)

       ! *** Clean up
       DEALLOCATE(ostate_p)
       
    END IF

  END SUBROUTINE obs_op_chl_cci


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
  SUBROUTINE init_dim_obs_l_chl_cci(domain_p, step, dim_obs, dim_obs_l)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l
    ! Include localization radius and local coordinates
    USE mod_assim_pdaf, ONLY: coords_l, locweight, loctype, &
                              mype_debug, node_debug
    ! Number of domains per sweep:
    USE g_parsup, ONLY: myDim_nod2D

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs      !< Full dimension of observation vector
    INTEGER, INTENT(inout) :: dim_obs_l  !< Local dimension of observation vector


    IF (thisobs%doassim == 1) THEN
       IF (loctype == 1) THEN
          ! *** Variable localization radius for fixed effective observation dimension ***
          CALL get_adaptive_lradius_pdaf(domain_p, lradius_chl_cci, loc_radius_chl_cci)
       END IF
       lradius_chl_cci = loc_radius_chl_cci(modulo(domain_p,myDim_nod2D))

!~        IF (mype_filter==mype_debug .AND. domain_p==node_debug) THEN
!~          CALL PDAFomi_set_debug_flag(domain_p)
!~        ELSE
!~          CALL PDAFomi_set_debug_flag(0)
!~        ENDIF
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: thisobs_l% dim_obs_l', thisobs_l% dim_obs_l
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: thisobs_l% id_obs_l', thisobs_l% id_obs_l
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: thisobs_l% distance_l', thisobs_l% distance_l
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: coords_l', coords_l
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: locweight', locweight
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: lradius_chl_cci', lradius_chl_cci
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: sradius_chl_cci', sradius_chl_cci
!~        if (mype_filter==mype_debug .and. domain_p==node_debug) write(*,*) 'Frauke: dim_obs_l', dim_obs_l
       
       ! ************************************************************
       ! *** Adapt observation error for coupled DA (double loop) ***
       ! ************************************************************
    
       if (n_sweeps>1) then
       
          ! Physics observations sweep.
          ! Set inverse observation error to small value
          if (domain_p==1) then

             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                   '--- PHY sweep: set ivar_obs_f for CHL to 1.0e-12'
             thisobs%ivar_obs_f = 1.0e-12
             
          ! BGC observations sweep.
          elseif (domain_p==myDim_nod2D+1) then
             if (mype_filter==0) &
                  write (*,'(a,4x,a)') 'FESOM-PDAF', &
                  '--- BIO sweep: set ivar_obs_f for CHL to original'
             thisobs%ivar_obs_f(:) = ivariance_obs_g
          end if
       end if ! n_sweeps

       ! **********************************************
       ! *** Initialize local observation dimension ***
       ! **********************************************

       CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
            locweight, lradius_chl_cci, sradius_chl_cci, dim_obs_l)


       
    END IF

  END SUBROUTINE init_dim_obs_l_chl_cci

END MODULE obs_chl_cci_pdafomi
