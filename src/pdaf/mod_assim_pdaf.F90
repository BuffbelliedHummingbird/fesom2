! PDAF in AWI-CM2 / Fesom 2.0

MODULE mod_assim_pdaf

! DESCRIPTION:
! This module provides variables needed for the 
! assimilation within the routines of the dummy model.
! For simplicity, all assimilation-related variables
! are stored here, even if they are only used in
! the main program for the filter initialization.
! Most variables can be specified as a command line 
! argument.
!
!
!
! REVISION HISTORY:
! 2013-02 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022-03 - Frauke       - Removed sea ice, added ocean velocities for FESOM2.1
!
! USES:
USE MOD_MESH
IMPLICIT NONE
SAVE
! *** Below are the generic variables used for configuring PDAF ***
! *** Their values are set in init_PDAF                         ***

! PUBLIC MEMBER FUNCTIONS:

! Settings for time stepping - available as namelist read-in
INTEGER :: step_null = 0 ! initial time step of assimilation
INTEGER :: days_since_DAstart = 1 ! days since start of assimilation (continuous counter through restarts)
! Settings for observations - available as command line options
INTEGER :: delt_obs_ocn     ! time step interval between assimilation steps - Ocean
INTEGER :: dim_obs          ! Number of observations
REAL    :: peak_obs_error   ! Peak value used to define the observation error
INTEGER :: proffiles_o      ! (0) don't generate profile observation files; 
                            ! (1) generate distributed profile files
                            ! (2) generate global profile file
INTEGER :: start_year_o, &  ! Which years to generate profile files
           end_year_o
! LOGICAL :: use_global_obs   ! Whether to use global full obs, of full obs limited to process domains
INTEGER :: use_global_obs
LOGICAL :: twin_experiment = .false.   ! Whether to perform a twin experiment with synthetic observations
INTEGER :: dim_obs_max      ! Expect max. number of observations for synthetic obs.

! General control of PDAF - available as command line options
INTEGER :: screen       ! Control verbosity of PDAF
                        ! (0) no outputs, (1) progess info, (2) add timings
                        ! (3) debugging output
INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
                        ! Number of EOFs to be used for SEEK
INTEGER :: filtertype   ! Select filter algorithm:
                        ! SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4), LETKF (5)
INTEGER :: subtype      ! Subtype of filter algorithm
                        !   SEIK:
                        !     (0) ensemble forecast; new formulation
                        !     (1) ensemble forecast; old formulation
                        !     (2) fixed error space basis
                        !     (3) fixed state covariance matrix
                        !     (4) SEIK with ensemble transformation
                        !   LSEIK:
                        !     (0) ensemble forecast;
                        !     (2) fixed error space basis
                        !     (3) fixed state covariance matrix
                        !     (4) LSEIK with ensemble transformation
                        !   ETKF:
                        !     (0) ETKF using T-matrix like SEIK
                        !     (1) ETKF following Hunt et al. (2007)
                        !       There are no fixed basis/covariance cases, as
                        !       these are equivalent to SEIK subtypes 2/3
                        !   LETKF:
                        !     (0) ETKF using T-matrix like SEIK
                        !     (1) LETKF following Hunt et al. (2007)
                        !       There are no fixed basis/covariance cases, as
                        !       these are equivalent to LSEIK subtypes 2/3
INTEGER :: incremental  ! Perform incremental updating in LSEIK
INTEGER :: dim_lag      ! Number of time instances for smoother
INTEGER :: DA_couple_type ! (0) for weakly-coupled, (1) for strongly-coupled assimilation
! Filter settings - available as command line options
! General
INTEGER :: type_forget  ! Type of forgetting factor
REAL    :: forget       ! Forgetting factor for filter analysis
INTEGER :: dim_bias     ! dimension of bias vector
REAL    :: varscale=1.0 ! Scaling factor for initial ensemble variance
! SEIK/ETKF/LSEIK/ETKFS
INTEGER :: type_trans    ! Type of ensemble transformation
                         ! SEIK/LSEIK:
                         ! (0) use deterministic omega
                         ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
                         ! (2) use product of (0) with random orthonormal matrix with
                         !     eigenvector (1,...,1)^T
                         ! ETKF/LETKF with subtype=4:
                         ! (0) use deterministic symmetric transformation
                         ! (2) use product of (0) with random orthonormal matrix with
                         !     eigenvector (1,...,1)^T
! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
                         !   (0) symmetric square root, (1) Cholesky decomposition
! Localization - LSEIK/LETKF/LESTKF
REAL    :: local_range   ! Range for local observation domain
INTEGER :: locweight     ! Type of localizing weighting of observations
                  !   (0) constant weight of 1
                  !   (1) exponentially decreasing with SRANGE
                  !   (2) use 5th-order polynomial
                  !   (3) regulated localization of R with mean error variance
                  !   (4) regulated localization of R with single-point error variance
REAL    :: srange        ! Support range for 5th order polynomial
                         !   or radius for 1/e for exponential weighting
! Specific for FESOM
INTEGER :: dim_state              ! Global size of model state
INTEGER :: dim_state_p            ! PE-local size of model state
INTEGER :: step_final             ! Final time step
INTEGER :: istep_asml             ! Time step at end of an forecast phase
LOGICAL :: flag_final=.false.     ! Whether the current is the final analysis step

! Declare Fortran type holding the indices of model fields in the state vector
TYPE field_ids
   INTEGER :: ssh
   INTEGER :: u 
   INTEGER :: v 
   INTEGER :: w 
   INTEGER :: temp 
   INTEGER :: salt
   INTEGER :: a_ice
   INTEGER :: MLD1
   INTEGER :: PhyChl
   INTEGER :: DiaChl
   INTEGER :: DIC
   INTEGER :: DOC
   INTEGER :: Alk
   INTEGER :: DIN
   INTEGER :: DON
   INTEGER :: O2
   INTEGER :: pCO2s
   INTEGER :: CO2f
   INTEGER :: PhyN
   INTEGER :: PhyC
   INTEGER :: DiaN
   INTEGER :: DiaC
   INTEGER :: DiaSi
   INTEGER :: PAR
   INTEGER :: NPPn
   INTEGER :: NPPd
   !   INTEGER :: TChl   ! Total chlorophyll = PhyChl + DiaChl
   !   INTEGER :: TDN    ! Total dissolved N = DIN + DON
   INTEGER :: Zo1C
   INTEGER :: Zo1N
   INTEGER :: Zo2C
   INTEGER :: Zo2N
   INTEGER :: DetC
   !   INTEGER :: TOC    ! Total organic carbon: PhyC + DiaC + DetC + DOC + HetC
   INTEGER :: PhyCalc
   INTEGER :: export
   INTEGER :: alphaCO2
   INTEGER :: PistonVel
END TYPE field_ids

! Type variable holding field IDs in state vector
TYPE(field_ids) :: id

INTEGER :: nfields          ! Number of fields in state vector
INTEGER :: phymin, phymax   ! First and last physics field in state vector
INTEGER :: bgcmin, bgcmax   ! First and last biogeochemistry field in state vector
! Specific for local filters
INTEGER, ALLOCATABLE :: id_lobs_in_fobs(:)     ! Indices of local observations in full obs. vector
INTEGER, ALLOCATABLE :: id_lstate_in_pstate(:) ! Indices of local state vector in PE-local global state vector
REAL, ALLOCATABLE    :: ivariance_obs_l(:)     ! Local inverse variance of observations
REAL, ALLOCATABLE :: distance(:)               ! Distances of local observations
! Variables for adaptive localization radius
REAL, ALLOCATABLE :: eff_dim_obs(:)            ! Effective observation dimension
REAL, ALLOCATABLE :: loc_radius(:)             ! Effective observation dimension
INTEGER :: loctype       ! Type of localization
                         !   (0) Fixed radius defined by local_range
                         !   (1) Variable radius for constant effective observation dimension
REAL :: loc_ratio        ! Choose local_range so the effective observation dim. is loc_ratio times dim_ens
INTEGER, ALLOCATABLE :: id_nod2D_ice(:)        ! IDs of nodes with ice
INTEGER :: depth_excl_no
INTEGER, ALLOCATABLE :: depth_excl(:)        ! nodes excluded in each pe 
! File output and input - available as as namelist read-in
LOGICAL :: read_inistate = .false.            ! Whether to read initial state from separate file
CHARACTER(len=120) :: DAoutput_path  = '.'      ! Path of DAoutput
CHARACTER(len=120) :: path_init = '.'         ! Path to initialization files
CHARACTER(len=120) :: file_init = 'covar_'    ! netcdf file holding distributed initial
                                              ! state and covariance matrix (added is _XX.nc)
CHARACTER(len=120) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                              ! state (added is _XX.nc)
CHARACTER(len=120) :: file_syntobs = 'syntobs.nc' ! File name for synthetic observations
CHARACTER(len=120) :: path_obs_rawprof  = ''      ! Path to profile observations
CHARACTER(len=120) :: file_rawprof_prefix  = ''   ! file name prefix for profile observations 
CHARACTER(len=120) :: file_rawprof_suffix  = '.nc'! file name suffix for profile observations 
LOGICAL :: ASIM_START_USE_CLIM_STATE = .true.
LOGICAL :: this_is_pdaf_restart = .false.

CHARACTER(len=120) :: path_atm_cov

! Other variables - NOT available as command line options / in the namelist:
REAL    :: time      ! model time
INTEGER, ALLOCATABLE :: offset(:)          ! PE-local offsets of fields in state vector
INTEGER, ALLOCATABLE :: dim_fields(:)      ! PE-local dimensions of fields in state vector
INTEGER, ALLOCATABLE :: offset_glob(:)     ! Global offsets of fields in state vector
INTEGER, ALLOCATABLE :: dim_fields_glob(:) ! Global dimensions of fields in state vector
REAL :: coords_l(2)                        ! Coordinates of local analysis domain
! REAL, PARAMETER :: pi=3.141592653589793
REAL, ALLOCATABLE :: state_fcst(:,:) ! State prior to assimilation, which is saved to use for correction
REAL, ALLOCATABLE :: var_p(:)        ! Estimated local model state variances
REAL, ALLOCATABLE :: std_p(:)        ! Estimated local model std variances
INTEGER :: num_day_in_month(0:1,12), endday_of_month_in_year(0:1,12), startday_of_month_in_year(0:1,12)
REAL, ALLOCATABLE :: monthly_state_f(:)       ! forecasted monthly state
REAL, ALLOCATABLE :: monthly_state_a(:)       ! analyzed monthly state
REAL, ALLOCATABLE :: monthly_state_m(:)       ! (analyzed) monthly time-mean state
REAL, ALLOCATABLE :: monthly_state_ens_f(:,:)
REAL, ALLOCATABLE :: monthly_state_ens_a(:,:)
LOGICAL :: write_monthly_mean =.false.   ! set to true if writing 3D fields monthly;
                                         ! otherwise, set to false to write daily 3D;
                                         ! 2D fields are always written daily.
INTEGER :: mon_snapshot_mem =0
REAL, ALLOCATABLE :: timemean(:)     ! Daily mean local state vector (analysis)

DATA num_day_in_month(0,:) /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
DATA num_day_in_month(1,:) /31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
DATA endday_of_month_in_year(0,:) /31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/
DATA endday_of_month_in_year(1,:) /31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366/
DATA startday_of_month_in_year(0,:) /1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/
DATA startday_of_month_in_year(1,:) /1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336/

type(t_mesh), pointer, save      :: mesh_fesom

! For debugging:
INTEGER :: debug_id_depth, & ! Location for debugging output
           debug_id_nod2           
INTEGER :: ens_member_debug
INTEGER :: mype_debug = 30
INTEGER :: node_debug = 888

END MODULE mod_assim_pdaf
