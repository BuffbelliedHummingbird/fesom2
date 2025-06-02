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

! Variables for time stepping:
INTEGER :: step_null = 0          ! at initialization of PDAF, specify which time step we have, to assimilate the time-corresponding observations:
                                  !    - if model (re)starts on 1st-Jan, step_null must be zero
                                  !    - if model (re)starts later during year, step_null must be adapted (done in slurm-job-script)
INTEGER :: days_since_DAstart = 1 ! days since start of assimilation; counted onwards at model restarts
INTEGER :: delt_obs_ocn           ! time step interval between assimilation steps
INTEGER :: istep_asml             ! assimilation time step at end of an forecast phase (FESOM's "mstep" + step_null);
                                  ! FESOM's "mstep" starts counting from zero even at each restart (even mid-of year);
                                  ! istep_asml starts counting from zero on 1st-Jan
INTEGER :: assim_time             ! assimlation step time of day (UTC) in seconds

! Settings for observations:
INTEGER :: dim_obs          ! Number of observations
REAL    :: peak_obs_error   ! Peak value used to define the observation error
INTEGER :: proffiles_o      ! (0) don't generate profile observation files; 
                            ! (1) generate distributed profile files
                            ! (2) generate global profile file
INTEGER :: start_year_o, &  ! Which years to generate profile files
           end_year_o
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
LOGICAL :: resetforget  ! Whether to reset forgetting factor after initial phase
INTEGER :: dim_bias     ! dimension of bias vector
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


! Declare Fortran type holding the indices of model fields in the state vector
TYPE field_ids
   INTEGER :: ssh        ! physics
   INTEGER :: u 
   INTEGER :: v 
   INTEGER :: w 
   INTEGER :: temp 
   INTEGER :: salt
   INTEGER :: a_ice
   INTEGER :: MLD1
   INTEGER :: MLD2
   INTEGER :: PhyChl     ! chlorophyll
   INTEGER :: DiaChl
   INTEGER :: DIC        ! dissolved tracers
   INTEGER :: DOC
   INTEGER :: Alk
   INTEGER :: DIN
   INTEGER :: DON
   INTEGER :: O2
   INTEGER :: pCO2s      ! surface carbon diagnostics
   INTEGER :: CO2f
   INTEGER :: alphaCO2
   INTEGER :: PistonVel
   INTEGER :: PhyN       ! small phyto
   INTEGER :: PhyC
   INTEGER :: PhyCalc
   INTEGER :: DiaN       ! diatoms
   INTEGER :: DiaC
   INTEGER :: DiaSi
   INTEGER :: Zo1C       ! zooplankton
   INTEGER :: Zo1N
   INTEGER :: Zo2C
   INTEGER :: Zo2N
   INTEGER :: DetC       ! detritus
   INTEGER :: DetCalc
   INTEGER :: DetSi
   INTEGER :: DetN
   INTEGER :: Det2C
   INTEGER :: Det2Calc
   INTEGER :: Det2Si
   INTEGER :: Det2N
   INTEGER :: PAR        ! diags
   INTEGER :: NPPn
   INTEGER :: NPPd
   INTEGER :: export
   INTEGER :: sigma
   !   INTEGER :: TChl   ! Total chlorophyll = PhyChl + DiaChl
   !   INTEGER :: TDN    ! Total dissolved N = DIN + DON
   !   INTEGER :: TOC    ! Total organic carbon: PhyC + DiaC + DetC + DOC + HetC
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
INTEGER, ALLOCATABLE :: depth_excl(:)          ! nodes excluded in each pe

! File output and input - available as as namelist read-in
LOGICAL :: read_inistate = .false.            ! Whether to read initial state from separate file
CHARACTER(len=150) :: DAoutput_path  = '.'    ! Path of DAoutput
CHARACTER(len=150) :: path_init = '.'         ! Path to initialization files
CHARACTER(len=150) :: file_init = 'covar_'    ! netcdf file holding distributed initial
                                              ! state and covariance matrix (added is _XX.nc)
CHARACTER(len=150) :: file_inistate = 'state_ini_' ! netcdf file holding distributed initial
                                                   ! state (added is _XX.nc)
CHARACTER(len=150) :: file_syntobs = 'syntobs.nc'  ! File name for synthetic observations
CHARACTER(len=150) :: path_obs_rawprof  = ''       ! Path to profile observations
CHARACTER(len=150) :: file_rawprof_prefix  = ''    ! file name prefix for profile observations 
CHARACTER(len=150) :: file_rawprof_suffix  = '.nc' ! file name suffix for profile observations 
LOGICAL :: ASIM_START_USE_CLIM_STATE = .true.

! Initial ensemble covariance
REAL    :: varscale=1.0           ! scaling factor for initial ensemble variance
! which fields to perturb
LOGICAL :: perturb_ssh   = .true.
LOGICAL :: perturb_u     = .true.
LOGICAL :: perturb_v     = .true.
LOGICAL :: perturb_temp  = .true.
LOGICAL :: perturb_salt  = .true.
LOGICAL :: perturb_DIC   = .true.
LOGICAL :: perturb_Alk   = .true.
LOGICAL :: perturb_DIN   = .true.
LOGICAL :: perturb_O2    = .true.

CHARACTER(len=150) :: path_atm_cov

! Restart information - set in slurm-job-script:
LOGICAL :: this_is_pdaf_restart = .false.            ! init_pdaf:        - at every start, initialize PDAF-netCDF-output
                                                     !                   - at restart, set forget from restart info
                                                     ! init_ens_pdaf:    - at restart, skip perturbation of initial fields
                                                     ! distribute_state: - at restart, skip distribution of initial fields
                                                     ! add_atmos_ens_st: - at restart, read perturbed atmospheric state
                                                     !                   - at restart, read forget and target value for RMSE
LOGICAL :: start_from_ENS_spinup = .false.           ! init_ens_pdaf:    - at start from perturbed ensemble, skip perturbation of initial fields
                                                     ! add_atmos_ens_st: - at start from perturbed ensemble, read perturbed atmospheric state
                                                     ! distribute_state: - at start from perturbed ensemble, skip distribution of initial fields


LOGICAL :: assimilateBGC = .false. ! whether to do a BGC assimilation step
LOGICAL :: assimilatePHY = .false. ! whether to do a physics assimilation step

! Other variables - NOT available as command line options / in the namelist:
REAL                 :: time               ! model time
INTEGER, ALLOCATABLE :: offset(:)          ! PE-local offsets of fields in state vector
INTEGER, ALLOCATABLE :: dim_fields(:)      ! PE-local dimensions of fields in state vector
INTEGER, ALLOCATABLE :: offset_glob(:)     ! Global offsets of fields in state vector
INTEGER, ALLOCATABLE :: dim_fields_glob(:) ! Global dimensions of fields in state vector
REAL                 :: coords_l(2)        ! Coordinates of local analysis domain
INTEGER, ALLOCATABLE :: dim_fields_l(:)    ! Field dimensions for local domain (i.e. field of vertical water column at 1 node)
INTEGER, ALLOCATABLE :: offset_l(:)        ! Field offsets for local domain

REAL, ALLOCATABLE :: state_fcst(:,:)       ! state prior to assimilation, saved to use for correction
REAL, ALLOCATABLE :: state_fcst_SSH_p(:,:) ! state prior to assimilation, saved to use for correction
REAL, ALLOCATABLE :: stdev_SSH_f_p(:)      ! forecast ensemble standard deviation at grid points for SSH field, saved to use for correction
REAL, ALLOCATABLE :: monthly_state_f(:)       ! forecasted monthly state
REAL, ALLOCATABLE :: monthly_state_a(:)       ! analyzed monthly state
REAL, ALLOCATABLE :: monthly_state_m(:)       ! monthly time-mean state
REAL, ALLOCATABLE :: monthly_state_ens_f(:,:)
REAL, ALLOCATABLE :: monthly_state_ens_a(:,:)
REAL, ALLOCATABLE :: timemean(:)     ! daily mean local state vector (mean of model forecast steps and analysis step)

REAL, ALLOCATABLE :: monthly_state_sf(:)       ! forecasted monthly standard deviation
REAL, ALLOCATABLE :: monthly_state_sa(:)       ! analyzed monthly standard deviation
REAL, ALLOCATABLE :: monthly_state_sm(:)       ! monthly time-mean standard deviation
REAL, ALLOCATABLE :: timemean_s(:)

! whether to compute monthly means:
LOGICAL :: compute_monthly_ff
LOGICAL :: compute_monthly_aa
LOGICAL :: compute_monthly_mm

LOGICAL :: compute_monthly_sf
LOGICAL :: compute_monthly_sa
LOGICAL :: compute_monthly_sm


! Julian-Gregorian date transformation of EN4 raw data 
INTEGER :: num_day_in_month(0:1,12), endday_of_month_in_year(0:1,12), startday_of_month_in_year(0:1,12)

! FESOM mesh:
type(t_mesh), pointer, save :: mesh_fesom
INTEGER, PARAMETER :: nlmax = 46            ! CORE2 mesh: deepest wet cells at mesh_fesom%nl-2
REAL, ALLOCATABLE :: topography3D(:,:)      ! topography: 1 for wet nodes and 0 for dry nodes (array shape as in model)
REAL, ALLOCATABLE :: topography_p(:)        ! """                                             (array shape as state_p)
REAL, ALLOCATABLE :: topography3D_g(:,:)    ! """                                             (array shape as in model globally)
REAL :: area_surf_glob(nlmax)               ! ocean area and standard volume to calculate area-/volume weighted means
REAL :: inv_area_surf_glob(nlmax)
REAL :: volo_full_glob, inv_volo_full_glob
REAL, ALLOCATABLE :: cellvol(:,:)           ! standard volume of cells, NOT considering time-varying ALE layerwidth


! For weak coupling:
integer :: n_sweeps                 !< Number of sweeps in local analysis loop
character(len=3) :: type_sweep(2)   !< Type of sweep in local analysis loop
integer :: isweep                   !< Index of sweep during the local analysis loop
character(len=6) :: cda_phy   ! Flag whether strongly-coupled DA is done
character(len=6) :: cda_bio   ! Flag whether strongly-coupled DA is done

! Initial state in case of restarts:
real, allocatable :: state_p_init(:)
real, allocatable :: ens_p_init(:,:)

! For carbon diagnostics:
real, allocatable :: factor_mass(:,:)
real, allocatable :: factor_conc(:,:)

! Type variable for postprocessing:
type obs_PP
   real, allocatable :: isExclObs (:)! whether to exclude observation due to model topography
   real, allocatable :: isInnoOmit(:)! whether to exclude observation due to Inno Omit
   real, allocatable :: nod1_g(:)    ! observation indeces on global FESOM grid (nodes)
   real, allocatable :: nod2_g(:)    ! """
   real, allocatable :: nod3_g(:)    ! """
   real, allocatable :: elem_g(:)    ! observation indeces on global FESOM grid (elements)
   real, allocatable :: lon(:)       ! observation coordinates
   real, allocatable :: lat(:)       ! """
   real, allocatable :: nz(:)        ! observation layer indeces
   real, allocatable :: depth(:)     ! observation depth
   real, allocatable :: numrep(:)    ! number of observations on one single element
   real, allocatable :: volelem(:)   ! volume of FESOM element
end type obs_PP


! For debugging:
INTEGER :: debug_id_depth, & ! Location for debugging output
           debug_id_nod2           
INTEGER :: ens_member_debug
INTEGER :: mype_debug = 18
INTEGER :: node_debug = 1167
! DIC debugging: Have_obs on mype_debug=18; node_debug=1167 on 2010-01-05.
! Alk debugging: Have_obs on mype_debug=65; node_debug=915 .OR. 885 .OR. 268 on 2010-01-08.

END MODULE mod_assim_pdaf
