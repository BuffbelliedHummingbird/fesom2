! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE init_pdaf()

! !DESCRIPTION:
! This routine collects the initialization of variables for PDAF.
! In addition, the initialization routine PDAF_init is called
! such that the internal initialization of PDAF is performed.
! This variant is for the online mode of PDAF.
!
! This routine is generic. However, it assumes a constant observation
! error (rms_obs). Further, with parallelization the local state
! dimension dim_state_p is used.
!
! !REVISION HISTORY:
! 2008-10 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022-02 - Frauke B     - Adapted for FESOM2.1
!
! !USES:
  USE mod_parallel_pdaf, &                                                ! Parallelization variables for assimilation
       ONLY: n_modeltasks, task_id, COMM_filter, COMM_couple, filterpe, &
       mype_world, COMM_model, abort_parallel, MPI_COMM_WORLD, MPIerr, &
       mype_model, mype_filter, npes_filter, writepe, mype_submodel, &
       MPI_INTEGER, MPI_SUM
  USE mod_assim_pdaf, &                                                   ! Variables for assimilation
       ONLY: dim_state, dim_state_p, dim_ens, dim_lag, &
       step_null, offset, dim_fields, screen, filtertype, subtype, &
       delt_obs_ocn, write_monthly_mean, &
       DA_couple_type, incremental, type_forget, &
       forget, locweight, local_range, srange, &
       type_trans, type_sqrt, eff_dim_obs, loc_radius, loctype, &
       twin_experiment, dim_obs_max, use_global_obs, mesh_fesom, nlmax, &
       offset_glob, dim_fields_glob, nfields, id, &
       phymin, phymax, bgcmin, bgcmax, &
       path_atm_cov, this_is_pdaf_restart, timemean, &
       ! Domain-local:
       offset_l, dim_fields_l, &
       ! Debugging:
       debug_id_depth, debug_id_nod2, ens_member_debug, &
       ! EN4 profile data processing:
       proffiles_o, path_obs_rawprof, file_rawprof_prefix, file_rawprof_suffix, &
       ! Weak coupling of FESOM-REcoM:
       n_sweeps, type_sweep, assimilateBGC, assimilatePHY, &
       cda_phy, cda_bio
  ! BGC parameter perturbation:
  USE mod_perturbation_pdaf, &
       ONLY: perturb_scale, perturb_parameters, perturb_lognormal, &
       perturb_scaleD
  USE recom_config, &
       ONLY: alfa, alfa_d, P_cm, P_cm_d, Chl2N_max, Chl2N_max_d, &
       deg_Chl, deg_Chl_d, graz_max, graz_max2, grazEff, grazEff2, &
       VDet, VDet_zoo2, Vdet_a, k_din, k_din_d, res_phy, res_phy_d, &
       rho_N, rho_C1, lossN, lossN_d, lossC, lossC_d, reminN, &
       reminC, calc_prod_ratio, res_het, res_zoo2, &
       biosynth, calc_diss_rate, calc_diss_rate2
  USE obs_sss_smos_pdafomi, &
       ONLY: assim_o_sss, rms_obs_sss, path_obs_sss, file_sss_prefix, file_sss_suffix, &
       sss_exclude_ice, sss_exclude_diff, bias_obs_sss, sss_fixed_rmse
  USE obs_sss_cci_pdafomi, &
       ONLY: assim_o_sss_cci, rms_obs_sss_cci, path_obs_sss_cci, file_sss_cci_prefix, file_sss_cci_suffix, &
       sss_cci_exclude_ice, sss_cci_exclude_diff, bias_obs_sss_cci, sss_cci_fixed_rmse
  USE obs_ssh_cmems_pdafomi, &
       ONLY: assim_o_ssh, rms_obs_ssh, path_obs_ssh, file_ssh_prefix, file_ssh_suffix, &
       ssh_exclude_ice, ssh_exclude_diff, bias_obs_ssh, ssh_fixed_rmse
  USE obs_sst_pdafomi, &
       ONLY: assim_o_sst, rms_obs_sst, path_obs_sst, file_sst_prefix, file_sst_suffix, &
       sst_exclude_ice, sst_exclude_diff, bias_obs_sst, sst_fixed_rmse
  USE obs_TSprof_EN4_pdafomi, &
       ONLY: assim_o_en4_s, assim_o_en4_t, &
       rms_obs_S, rms_obs_T, &
       path_obs_prof, file_prof_prefix, file_prof_suffix, &
       bias_obs_prof, prof_exclude_diff
  USE mod_atmos_ens_stochasticity, &
       ONLY: disturb_xwind, disturb_ywind, disturb_humi, &
       disturb_qlw, disturb_qsr, disturb_tair, &
       disturb_prec, disturb_snow, disturb_mslp, &
       init_atmos_ens_stochasticity, init_atmos_stochasticity_output,&
       atmos_stochasticity_ON

  USE mod_obs_f_pdaf, &
       ONLY: get_domain_limits_unstr
  USE g_PARSUP, &
       ONLY: myDim_nod2D, MPI_COMM_FESOM, myList_edge2D, myDim_edge2D, myDim_elem2D, mype
  USE MOD_MESH
  USE g_clock, &
       ONLY: timeold, daynew, cyearold, yearnew, yearold
  USE g_rotate_grid, &
       ONLY: r2g
  USE timer, only: timeit
  USE PDAFomi, &
       ONLY: PDAFomi_get_domain_limits_unstr
  USE mod_nc_out_routines, &
       ONLY: netCDF_init, netCDF_STD_init
  USE mod_nc_out_variables, &
       ONLY: write_ens, init_sfields, sfields
  USE mod_carbon_fluxes_diags, &
       ONLY: init_carbonfluxes_diags_out


  IMPLICIT NONE

! !CALLING SEQUENCE:
! Called by: main
! Calls: init_pdaf_parse
! Calls: init_pdaf_info
! Calls: PDAF_init
! Calls: PDAF_get_state
!EOP

! Local variables
  INTEGER :: i,b               ! Counter
  INTEGER :: filter_param_i(7) ! Integer parameter array for filter
  REAL    :: filter_param_r(2) ! Real parameter array for filter
  INTEGER :: status_pdaf       ! PDAF status flag
  INTEGER :: doexit, steps     ! Not used in this implementation
  REAL    :: timenow           ! Not used in this implementation
  INTEGER :: show_cpus_for_barrier
  INTEGER :: iseed(4)          ! Seed for random number generation to perturb BGC parameters
  character(len=6) :: cdaval   ! Flag whether strongly-coupled DA is done
  
  REAL    :: dummy

  ! External subroutines
  EXTERNAL :: init_ens_pdaf            ! Ensemble initialization
  EXTERNAL :: next_observation_pdaf, & ! Provide time step, model time, 
                                       ! and dimension of next observation
       distribute_state_pdaf, &        ! Routine to distribute a state vector to model fields
       prepoststep_pdaf                ! User supplied pre/poststep routine

  dummy = 2.0

! ***************************
! ***   Initialize PDAF   ***
! ***************************

  ! Get process-ID in task of model compartment
  CALL MPI_Comm_Rank(MPI_COMM_FESOM, mype_submodel, MPIerr)

  IF (mype_submodel==0) THEN
     WRITE (*,'(1x,a, i5)') 'FESOM-PDAF: INITIALIZE PDAF, task: ', task_id
  END IF


! **********************************************************
! ***                  CONTROL OF PDAF                   ***
! ***              used in call to PDAF_init             ***
! **********************************************************

! *** IO options ***
  screen     = 2    ! Write screen output (1) for output, (2) add timings

! *** Filter specific variables
  filtertype = 7    ! Type of filter
                    !   (6) ESTKF
                    !   (7) LESTKF
                    !   (11) GENOBS: Generate synthetic observations
  dim_ens = n_modeltasks ! Size of ensemble for all ensemble filters
                    ! Number of EOFs to be used for SEEK
  dim_lag = 0       ! Size of lag in smoother
  subtype = 0       ! subtype of filter: 
                    !   ESTKF:
                    !     (0) Standard form of ESTKF
                    !   LESTKF:
                    !     (0) Standard form of LESTKF
  type_trans = 0    ! Type of ensemble transformation
                    !   SEIK/LSEIK and ESTKF/LESTKF:
                    !     (0) use deterministic omega
                    !     (1) use random orthonormal omega orthogonal to (1,...,1)^T
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
                    !   ETKF/LETKF:
                    !     (0) use deterministic symmetric transformation
                    !     (2) use product of (0) with random orthonormal matrix with
                    !         eigenvector (1,...,1)^T
  type_forget = 0   ! Type of forgetting factor in SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
                    !   (0) fixed
                    !   (1) global adaptive
                    !   (2) local adaptive for LSEIK/LETKF/LESTKF
  forget  = 1.0     ! Forgetting factor
  type_sqrt = 0     ! Type of transform matrix square-root
                    !   (0) symmetric square root, (1) Cholesky decomposition
  incremental = 0   ! (1) to perform incremental updating (only in SEIK/LSEIK!)
 

! **********************************************************
! ***     CONTROL OF USER ROUTINES FOR ASSIMILATION      ***
! **********************************************************

! *** Forecast length (time interval between analysis steps) ***
  delt_obs_ocn = 32  ! Number of time steps between analysis/assimilation steps

! *** Set weakly- or strongly-coupled DA
  DA_couple_type = 0 ! (0) for weakly- (1) for strongly-coupled DA

! *** Whether to generate profile observation files
  proffiles_o = 0  ! (0) don't generate them; 
                   ! (1) generate distributed profile files
                   ! (2) generate global profile file
                   
! *** Set assimilation variables
  assim_o_sst     = .false.
  assim_o_sss     = .false.
  assim_o_sss_cci = .false.
  assim_o_ssh     = .false.

  write_ens = .false.
  
  step_null = 0

! *** specifications for observations ***
  ! This error is the standard deviation for the Gaussian distribution 
  rms_obs_sst = 0.8 ! error for satellite SST observations
  rms_obs_sss = 0.5 ! error for satellite SSS observations (SMOS)
  rms_obs_sss_cci = 0.5 ! error for satellite SSS observations (CCI)
  rms_obs_ssh = 0.05 ! error for satellite SSH observations
  rms_obs_T = 0.8    ! error for temperature profile observations
  rms_obs_S = 0.5    ! error for salinity profile observations
  bias_obs_ssh = 0.0    ! observation bias  
  bias_obs_prof = 0.0   ! observation bias  
  sst_exclude_ice = .true.  ! Exclude SST observations at point with sea ice and T>0
  sst_exclude_diff = 0.0     ! Exclude SST observations if difference from ensemble mean is >sst_exclude_diff
  prof_exclude_diff = 0.0    ! Exclude profile T observations if difference from ensemble mean is >prof_exclude_diff
  use_global_obs = 1 ! Use global full obs. (1) or full obs. limited to process domains (0)
  twin_experiment = .false.  ! Whether to run a twin experiment assimilating synthetic observations
  dim_obs_max = 80000        ! Expected maximum number of observations for synthetic obs.

! *** Localization settings
  locweight = 0     ! Type of localizating weighting
                    !   (0) constant weight of 1
                    !   (1) exponentially decreasing with SRANGE
                    !   (2) use 5th-order polynomial
                    !   (3) regulated localization of R with mean error variance
                    !   (4) regulated localization of R with single-point error variance
  local_range = 0  ! Range in grid points for observation domain in local filters
  srange = local_range  ! Support range for 5th-order polynomial
                    ! or range for 1/e for exponential weighting

! *** File names - available as namelist read-in
  path_obs_sst = ''        ! Path to SST observation files
  file_sst_prefix = ''     ! Prefix of file holding SST observations
  file_sst_suffix = '.nc'  ! Suffix of file SST observations
  
  path_obs_sss = ''        
  file_sss_suffix = '.nc'  
  file_sss_prefix = ''

  path_obs_sss_cci = ''
  file_sss_cci_suffix = '.nc'
  file_sss_cci_prefix = ''  
  
  path_obs_ssh = ''
  file_ssh_prefix = ''     
  file_ssh_suffix = '.nc'  

  path_obs_prof = ''       ! Path to file holding profile observations
  file_prof_prefix = ''    ! Prefix of file holding profile observations
  file_prof_suffix = '.nc' ! Suffix of file holding profile observations
  
  path_obs_rawprof = ''       ! Path to file holding raw profile observations
  file_rawprof_prefix = ''    ! Prefix of file holding rawprofile observations
  file_rawprof_suffix = '.nc' ! Suffix of file holding raw profile observations

! *** Configuration for atmospheric stochasticity:
disturb_xwind=.true.
disturb_ywind=.true.
disturb_humi=.true.
disturb_qlw=.true.
disturb_qsr=.true.
disturb_tair=.true.
disturb_prec=.true.
disturb_snow=.true.
disturb_mslp=.true.

! *** Read PDAF configuration from namelist ***
  CALL read_config_pdaf()
  
! **************************************************************************
! *** Configuration for FESOM-REcom coupling: Define local analysis loop ***
! **************************************************************************

    cda_phy = 'weak'
    cda_bio = 'weak'

    if ((assimilateBGC) .and. (assimilatePHY)) then
       ! Observations of both physics and BGC are assimilated
       n_sweeps = 2
       type_sweep(1) = 'phy'
       type_sweep(2) = 'bio'
    else
       ! Less than two observation categories
       n_sweeps = 1
       if (assimilatePHY) then
          ! Only observations of physics are assimilated
          type_sweep(1) = 'phy'
       elseif (assimilateBGC) then
          ! Only observations of BGC are assimilated
          type_sweep(1) = 'bio'
       else
          ! No observation active (free run); set sweep to physics
          type_sweep(1) = 'phy'
       end if
    end if

    if (mype_world == 0) then
       write (*,'(a,2x,a)') 'FESOM-PDAF', '*** Setup for coupled DA FESOM-REcoM ***'
       write (*, '(a,4x,a,i5)') 'FESOM-PDAF', 'Number of local analysis sweeps', n_sweeps
       write (*, '(a,4x,a)') 'FESOM-PDAF','Type of sweeps:'
       do i = 1, n_sweeps
          if (trim(type_sweep(i))=='phy') then
             cdaval = cda_phy
          else
             cdaval = cda_bio
          end if
          write (*, '(a,8x,a,i3,3x,a,a,3x,a,a)') &
               'FESOM-PDAF', 'sweep', i, ' observation type: ', trim(type_sweep(i)), 'CDA: ', trim(cdaval)
       end do
    end if

! ***********************************************************************************************
! ***   For weakly-coupled assimilation of FESOM and atmosphere re-define filter communicator ***
! ***********************************************************************************************

  IF (DA_couple_type == 0) THEN

     !//TODO: here is the switcher for strongly coupled DA
     ! Set filter communicator to the communicator of FESOM
     
     COMM_filter = MPI_COMM_FESOM

     IF (filterpe) THEN
        CALL MPI_Comm_Size(COMM_filter, npes_filter, MPIerr)
        CALL MPI_Comm_Rank(COMM_filter, mype_filter, MPIerr)

        IF (mype_filter==0) THEN
           WRITE (*,'(a)') 'FESOM-PDAF: Initialize weakly-coupled data assimilation'
        ENDIF
     ENDIF
  ELSE
     IF (filterpe) THEN
        IF (mype_filter==0) THEN
           WRITE (*,'(a)') 'FESOM-PDAF: Initialize strongly-coupled data assimilation'
        END IF
     END IF
  END IF

! ***************************
! *** Define state vector ***
! ***************************
              
  id% ssh    =  1 ! sea surface height
  id% u      =  2 ! zonal velocity
  id% v      =  3 ! meridional velocity
  id% w      =  4 ! vertical velocity
  id% temp   =  5 ! temperature
  id% salt   =  6 ! salinity
  id% a_ice  =  7 ! sea-ice concentration
  id% MLD1   =  8 ! boundary layer depth (criterion after Large et al., 1997)
  id% MLD2   =  9 ! mixed layer depth (density treshold)
  
  id% PhyChl = 10 ! chlorophyll-a small phytoplankton
  id% DiaChl = 11 ! chlorophyll-a diatoms
  
  id% DIC    = 12 ! dissolved tracers
  id% DOC    = 13
  id% Alk    = 14
  id% DIN    = 15
  id% DON    = 16
  id% O2     = 17
  
  id% pCO2s     = 18 ! surface carbon diags
  id% CO2f      = 19
  id% alphaCO2  = 20
  id% PistonVel = 21
    
  id% PhyN   = 22 ! small phyto
  id% PhyC   = 23
  id% PhyCalc= 24
  
  id% DiaN   = 25 ! diatoms
  id% DiaC   = 26
  id% DiaSi  = 27
  
  id% Zo1N   = 28 ! zooplankton
  id% Zo2C   = 29
  id% Zo2N   = 30
  id% Zo1C   = 31
  
  id% DetC      = 32 ! detritus
  id% DetCalc   = 33
  id% DetSi     = 34
  id% DetN      = 35
  id% Det2C     = 36
  id% Det2Calc  = 37
  id% Det2Si    = 38
  id% Det2N     = 39

  id% PAR    = 40 ! diags
  id% NPPn   = 41
  id% NPPd   = 42
  id% export = 43
  
  nfields = 43

  phymin = 1
  phymax = 9
  
  bgcmin = 10
  bgcmax = nfields
  
  CALL init_sfields()

! *** Specify offset of fields in pe-local state vector ***

!    . . . A . . . . . B . . . . . C . . . . . D
!         . .         / .         / .         . .
!        .   .   2   /   .   3   /   .   5   .   .
! .     .     .     /     .     /     .     .     .
!  .   .   1   .   /   3   .   /   4   .   .       
!   . .         . /         . /         . .        
!    A . . . . . B . . . . . C . . . . . D . . . . 
!
!  A:    Internal nodes of left PE
!  B:    Internal nodes of left PE, simultanesously external nodes of right PE
!  C:    External nodes of left PE, simultanesously internal nodes of right PE
!  D:    Internal nodes of right PE
!  1:    Internal element of left PE, simultanesously wide-halo element of right PE (shares node B with right PE)
!  2:    Internal element of left PE, simultanesously small-halo element of right PE (shares edge BB with right PE)
!  3:    Internal (shared) elements of both PEs
!  4:    Small-halo element of left PE (shares edge CC with left PE), simultanesously internal element of right PE
!  5:    Wide-halo element of left PE (shares node C with left PE), simultanesously internal element of right PE
!
!  myDim_nod2D:      Number of internal nodes (A+B)
!  eDim_nod2D:       Number of external nodes (C)
!
!  myDim_elem2D:     Number of internal elements (1+2+3)
!  eDim_elem2D:      Number of small-halo elements (4)
!  xDim_elem2D:      Number of wide-halo elements (5)
!
!  mesh_fesom%nl:    Maximum number of fesom levels (1 is air-sea interface)
!  mesh_fesom%nl-1:  Maximum number of fesom layers (1 is surface layer, 0-5m)

  nlmax = mesh_fesom%nl - 2 ! CORE2 mesh: deepest wet cells at mesh_fesom%nl-2

  ALLOCATE(dim_fields(nfields))
  dim_fields(id% ssh   )   = myDim_nod2D           ! 1 SSH
  dim_fields(id% u     )   = myDim_nod2D*(nlmax)   ! 2 u (interpolated on nodes)
  dim_fields(id% v     )   = myDim_nod2D*(nlmax)   ! 3 v (interpolated on nodes)
  dim_fields(id% w     )   = myDim_nod2D*(nlmax)   ! 4 w
  dim_fields(id% temp  )   = myDim_nod2D*(nlmax)   ! 5 temp
  dim_fields(id% salt  )   = myDim_nod2D*(nlmax)   ! 6 salt
  dim_fields(id% a_ice )   = myDim_nod2D           ! 7 a_ice
  dim_fields(id% MLD1  )   = myDim_nod2D
  dim_fields(id% MLD2  )   = myDim_nod2D
  
  ! dim_fields biogeochemistry:
  do b=bgcmin,bgcmax
    ! 3D fields:
    if     (sfields(b)% ndims == 2) then
      dim_fields(b) = myDim_nod2D*(nlmax)
    ! surface fields:
    elseif (sfields(b)% ndims == 1) then
      dim_fields(b) = myDim_nod2D
    endif
  enddo

  ALLOCATE(offset(nfields))
  offset(id% ssh   )   = 0                                                ! 1 SSH
  offset(id% u     )   = offset(id% u     -1) + dim_fields(id% u     -1)  ! 2 u
  offset(id% v     )   = offset(id% v     -1) + dim_fields(id% v     -1)  ! 3 v
  offset(id% w     )   = offset(id% w     -1) + dim_fields(id% w     -1)  ! 4 w
  offset(id% temp  )   = offset(id% temp  -1) + dim_fields(id% temp  -1)  ! 5 temp
  offset(id% salt  )   = offset(id% salt  -1) + dim_fields(id% salt  -1)  ! 6 salt
  offset(id% a_ice )   = offset(id% a_ice -1) + dim_fields(id% a_ice -1)  ! 7 a_ice
  offset(id% MLD1  )   = offset(id% MLD1  -1) + dim_fields(id% MLD1  -1)
  offset(id% MLD2  )   = offset(id% MLD2  -1) + dim_fields(id% MLD2  -1)
  
  ! offset biogeochemistry
  do b=bgcmin,bgcmax
	offset(b)   = offset(b-1) + dim_fields(b-1)
  enddo
  
  dim_state_p = sum(dim_fields)
  
! *** Domain/node-local dimensions ***
    ALLOCATE(dim_fields_l(nfields)) ! Field dimensions for node-local domain; their value is set in init_dim_l_pdaf.
    ALLOCATE(offset_l(nfields))     ! Field offsets for node-local domain; their value is set in init_dim_l_pdaf.
    
  if (mype_model==0 .and. task_id==1) then
  do b=1,nfields
      WRITE (*,'(1X,A30,1X,A5,1X,I7,1X,I7)') 'PE0-dim and offset of field ', &
                                              trim(sfields(b)%variable), &
                                              dim_fields(b), &
                                              offset(b)
  enddo
  endif

! *******************************************************
! *** Specify offset of fields in global state vector ***
! *******************************************************

	ALLOCATE(dim_fields_glob(nfields))
	dim_fields_glob(id% ssh   ) = mesh_fesom%nod2D              ! SSH
	dim_fields_glob(id% u     ) = mesh_fesom%nod2D * (nlmax)    ! u (interpolated on nodes)
	dim_fields_glob(id% v     ) = mesh_fesom%nod2D * (nlmax)    ! v (interpolated on nodes)
	dim_fields_glob(id% w     ) = mesh_fesom%nod2D * (nlmax)    ! w
	dim_fields_glob(id% temp  ) = mesh_fesom%nod2D * (nlmax)    ! temp
	dim_fields_glob(id% salt  ) = mesh_fesom%nod2D * (nlmax)    ! salt
	dim_fields_glob(id% a_ice ) = mesh_fesom%nod2D              ! a_ice
	dim_fields_glob(id% MLD1  ) = mesh_fesom%nod2D
	dim_fields_glob(id% MLD2  ) = mesh_fesom%nod2D
	
	! dim_fields biogeochemistry:
    do b=bgcmin,bgcmax
      ! 3D fields:
      if     (sfields(b)% ndims == 2) then
        dim_fields_glob(b) = mesh_fesom%nod2D * (nlmax)
      ! surface fields:
      elseif (sfields(b)% ndims == 1) then
        dim_fields_glob(b) = mesh_fesom%nod2D
      endif
    enddo
	
	ALLOCATE(offset_glob(nfields))
	offset_glob(id% ssh   ) = 0                                                                 ! SSH
	offset_glob(id% u     ) = offset_glob(id% u     -1) + dim_fields_glob(id% u     -1)         ! u
	offset_glob(id% v     ) = offset_glob(id% v     -1) + dim_fields_glob(id% v     -1)         ! v
	offset_glob(id% w     ) = offset_glob(id% w     -1) + dim_fields_glob(id% w     -1)         ! w
	offset_glob(id% temp  ) = offset_glob(id% temp  -1) + dim_fields_glob(id% temp  -1)         ! temp
	offset_glob(id% salt  ) = offset_glob(id% salt  -1) + dim_fields_glob(id% salt  -1)         ! salt
	offset_glob(id% a_ice ) = offset_glob(id% a_ice -1) + dim_fields_glob(id% a_ice -1)         ! a_ice
	offset_glob(id% MLD1  ) = offset_glob(id% MLD1  -1) + dim_fields_glob(id% MLD1  -1)
	offset_glob(id% MLD2  ) = offset_glob(id% MLD2  -1) + dim_fields_glob(id% MLD2  -1)
	
	! offset_glob biogeochemistry
    do b=bgcmin,bgcmax
	  offset_glob(b)   = offset_glob(b-1) + dim_fields_glob(b-1)
    enddo
	
	dim_state = sum(dim_fields_glob)

  
  if (mype_world==0) THEN
      WRITE(*,*) 'init_pdaf - myDim_elem2D:      ', myDim_elem2D
	  WRITE(*,*) 'init_pdaf - mesh_fesom%elem2D: ', mesh_fesom%elem2D
	  WRITE(*,*) 'init_pdaf - mesh_fesom%nl:     ', mesh_fesom%nl
	  WRITE(*,*) 'init_pdaf - nlmax:             ', nlmax
	  WRITE(*,*) 'init_pdaf - dim_state:         ', dim_state
	  WRITE(*,*) 'init_pdaf - dim_state_p:       ', dim_state_p
  endif

! *** Initial Screen output ***

  IF (mype_model==0 .AND. task_id==1) CALL init_pdaf_info()

! *** Check ensemble size
  IF (dim_ens /= n_modeltasks) THEN
     WRITE (*,*) 'ERROR: Ensemble size (',dim_ens, &
          ') needs to be identical to the number of model tasks (',n_modeltasks,')'
     CALL abort_parallel()
  END IF
  
  ALLOCATE(timemean(dim_state_p))
  timemean = 0.0
  
! *************************************************
! *** Perturb model parameters for ensemble run ***
! *************************************************
! Selected parameters to disturb chlorophyll and biomass:
! alfa       = 0.14   ! Initial slope of P I curve small phytoplankton
! alfa_d     = 0.19   ! Initial slope of P I curve Diatoms
! P_cm        = 3.0   ! Small phytoplankton maximum rate of phtotosynthesis
! P_cm_d      = 3.5   ! Diatom maximum rate of phtotosynthesis
! Chl2N_max   = 3.78  ! Small phytoplankton maximum Chlorophyll a to nitrogen ratio
! Chl2N_max_d = 4.2   ! Diatom maximum Chlorophyll a to nitrogen ratio
! deg_Chl     = 0.1   ! degradation rate constant
! deg_Chl_d   = 0.1   ! degradation rate constant
! graz_max, graz_max2 ! maximum grazing rates
! grazEff, grazEff2   ! grazing effeciency of ZooPlankton
! Selected parameters to disturb dissolved tracers, most importantly DIC:
! VDet, VDet_zoo2, Vdet_a   ! Sinking speed detritus
! k_din, k_din_d            ! Nitrate uptake
! res_phy, res_phy_d        ! Respiration
! res_het, res_zoo2
! biosynth, (biosynthSi=0)
! calc_diss_rate, calc_diss_rate2 ! Dissolution during sinking
! rho_N, rho_C1             ! Remineralization of dissolved organic material
!                             (DON -> DIN; DOC --> DIC)
! lossN, lossN_d            ! Loss terms phytoplankton (PhyN --> DON)
! lossC, lossC_d
! reminN, reminC            ! Remineralization of detritus
! calc_prod_ratio           ! How much of small phytoplankton are calcifiers


IF (perturb_parameters) THEN

        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'FESOM-PDAF', 'Perturbing BGC parameters'

        ! iseed is the seed of the random number generator
        ! be aware that elements must be between 0 and 4095,
        ! and ISEED(4) must be odd.

        ! alfa
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'alfa  = ', alfa, ' (unperturbed)'
        iseed(1)=1
        iseed(2)=3
        iseed(3)=31+task_id ! different seed for each ensemble member
        iseed(4)=37
        CALL perturb_lognormal(alfa, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'alfa  = ', alfa, ' on task ', task_id
        
        ! alfa_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'alfa_d  = ', alfa_d, ' (unperturbed)'
        iseed(1)=2
        iseed(2)=5
        iseed(3)=29+task_id
        iseed(4)=41
        CALL perturb_lognormal(alfa_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'alfa_d  = ', alfa_d, ' on task ', task_id
       
        ! P_cm
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'P_cm  = ', P_cm, ' (unperturbed)'
        iseed(1)=2
        iseed(2)=4
        iseed(3)=23 + task_id
        iseed(4)=43
        CALL perturb_lognormal(P_cm, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'P_cm  = ', P_cm, ' on task ', task_id
        
        ! P_cm_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'P_cm_d  = ', P_cm_d, ' (unperturbed)'
        iseed(1)=7
        iseed(2)=4
        iseed(3)=19 + task_id
        iseed(4)=47
        CALL perturb_lognormal(P_cm_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'P_cm_d  = ', P_cm_d, ' on task ', task_id
        
        ! Chl2N_max
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'Chl2N_max  = ', Chl2N_max, ' (unperturbed)'
        iseed(1)=14
        iseed(2)=15
        iseed(3)=17 + task_id
        iseed(4)=53
        CALL perturb_lognormal(Chl2N_max, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'Chl2N_max  = ', Chl2N_max, ' on task ', task_id
        
        ! Chl2N_max_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'Chl2N_max_d  = ', Chl2N_max_d, ' (unperturbed)'
        iseed(1)=15
        iseed(2)=16
        iseed(3)=13 + task_id
        iseed(4)=59
        CALL perturb_lognormal(Chl2N_max_d,perturb_scale,iseed)
        IF (mype_model==0) WRITE(*,*) 'Chl2N_max_d  = ', Chl2N_max_d, ' on task ', task_id
        
        ! deg_Chl
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'deg_Chl  = ', deg_Chl, ' (unperturbed)'
        iseed(1)=8
        iseed(2)=6
        iseed(3)=11 + task_id
        iseed(4)=61
        CALL perturb_lognormal(deg_Chl, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'deg_Chl  = ', deg_Chl, ' on task ', task_id
        
        ! deg_Chl_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'deg_Chl_d  = ', deg_Chl_d, ' (unperturbed)'
        iseed(1)=104
        iseed(2)=97
        iseed(3)=7 + task_id
        iseed(4)=67
        CALL perturb_lognormal(deg_Chl_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'deg_Chl_d  = ', deg_Chl_d, ' on task ', task_id
        
        ! graz_max
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'graz_max  = ', graz_max, ' (unperturbed)'
        iseed(1)=33
        iseed(2)=16
        iseed(3)=3 + task_id
        iseed(4)=73
        CALL perturb_lognormal(graz_max, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'graz_max  = ', graz_max, ' on task ', task_id
        
        ! graz_max2
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'graz_max2  = ', graz_max2, ' (unperturbed)'
        iseed(1)=3
        iseed(2)=6
        iseed(3)=5 + task_id
        iseed(4)=71
        CALL perturb_lognormal(graz_max2, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'graz_max2  = ', graz_max2, ' on task ', task_id

        ! grazEff
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'grazEff  = ', grazEff, ' (unperturbed)'
        iseed(1)=2
        iseed(2)=1
        iseed(3)=6+3*task_id
        iseed(4)=13
        CALL perturb_lognormal(grazEff, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'grazEff  = ', grazEff, ' on task ', task_id
        
        ! grazEff2
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'grazEff2  = ', grazEff2, ' (unperturbed)'
        iseed(1)=6
        iseed(2)=3
        iseed(3)=97+101*task_id
        iseed(4)=11
        CALL perturb_lognormal(grazEff2, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'grazEff2  = ', grazEff2, ' on task ', task_id
        
        ! VDet
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'VDet  = ', VDet, ' (unperturbed)'
        iseed(1)=68
        iseed(2)=54
        iseed(3)=63+3*task_id
        iseed(4)=21
        CALL perturb_lognormal(VDet, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'VDet  = ', VDet, ' on task ', task_id
        
        ! VDet_zoo2
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'vdet_zoo2  = ', vdet_zoo2, ' (unperturbed)'
        iseed(1)=87
        iseed(2)=28
        iseed(3)=93+5*task_id
        iseed(4)=23
        CALL perturb_lognormal(vdet_zoo2, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'vdet_zoo2  = ', vdet_zoo2, ' on task ', task_id
        
        ! Vdet_a
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'vdet_a  = ', vdet_a, ' (unperturbed)'
        iseed(1)=45
        iseed(2)=65
        iseed(3)=41+7*task_id
        iseed(4)=27
        CALL perturb_lognormal(vdet_a, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'vdet_a  = ', vdet_a, ' on task ', task_id
        
        ! k_din
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'k_din  = ', k_din, ' (unperturbed)'
        iseed(1)=26
        iseed(2)=87
        iseed(3)=56+4*task_id
        iseed(4)=29
        CALL perturb_lognormal(k_din, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'k_din  = ', k_din, ' on task ', task_id
        
        ! k_din_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'k_din_d  = ', k_din_d, ' (unperturbed)'
        iseed(1)=61
        iseed(2)=15
        iseed(3)=14+3*task_id
        iseed(4)=33
        CALL perturb_lognormal(k_din_d, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'k_din_d  = ', k_din_d, ' on task ', task_id
        
        ! res_phy
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'res_phy  = ', res_phy, ' (unperturbed)'
        iseed(1)=67
        iseed(2)=50
        iseed(3)=31+2*task_id
        iseed(4)=35
        CALL perturb_lognormal(res_phy, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'res_phy  = ', res_phy, ' on task ', task_id
        
        ! res_phy_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'res_phy_d  = ', res_phy_d, ' (unperturbed)'
        iseed(1)=20
        iseed(2)=50
        iseed(3)=30+9*task_id
        iseed(4)=37
        CALL perturb_lognormal(res_phy_d, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'res_phy_d  = ', res_phy_d, ' on task ', task_id
        
        ! rho_N
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'rho_N  = ', rho_N, ' (unperturbed)'
        iseed(1)=12
        iseed(2)=24
        iseed(3)=36+11*task_id
        iseed(4)=45
        CALL perturb_lognormal(rho_N, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'rho_N  = ', rho_N, ' on task ', task_id
        
        ! rho_C1
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'rho_c1  = ', rho_c1, ' (unperturbed)'
        iseed(1)=9
        iseed(2)=38
        iseed(3)=44+3*task_id
        iseed(4)=47
        CALL perturb_lognormal(rho_c1, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'rho_c1  = ', rho_c1, ' on task ', task_id
        
        ! lossN
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'lossN  = ', lossn, ' (unperturbed)'
        iseed(1)=18
        iseed(2)=51
        iseed(3)=54+8*task_id
        iseed(4)=53
        CALL perturb_lognormal(lossn, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'lossN  = ', lossn, ' on task ', task_id
        
        ! lossN_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'lossN_d  = ', lossN_d, ' (unperturbed)'
        iseed(1)=48
        iseed(2)=73
        iseed(3)=35+5*task_id
        iseed(4)=55
        CALL perturb_lognormal(lossN_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'lossN_d  = ', lossN_d, ' on task ', task_id
        
        ! lossC
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'lossC  = ', lossC, ' (unperturbed)'
        iseed(1)=94
        iseed(2)=2
        iseed(3)=22+12*task_id
        iseed(4)=63
        CALL perturb_lognormal(lossC, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'lossC  = ', lossc, ' on task ', task_id
        
        ! lossC_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'lossC_d  = ', lossC_d, ' (unperturbed)'
        iseed(1)=41
        iseed(2)=74
        iseed(3)=31+8*task_id
        iseed(4)=67
        CALL perturb_lognormal(lossC_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,*) 'lossC_d  = ', lossc_d, ' on task ', task_id
        
        ! reminN
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'reminN  = ', reminn, ' (unperturbed)'
        iseed(1)=80
        iseed(2)=25
        iseed(3)=97+101*task_id
        iseed(4)=69
        CALL perturb_lognormal(reminn, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'reminN  = ', reminn, ' on task ', task_id
        
        ! reminC    
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'reminC  = ', reminC, ' (unperturbed)'
        iseed(1)=13
        iseed(2)=15
        iseed(3)=66+17*task_id
        iseed(4)=71
        CALL perturb_lognormal(reminC, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'reminC  = ', reminc, ' on task ', task_id
        
        ! calc_prod_ratio   
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'calc_prod_ratio  = ', calc_prod_ratio, ' (unperturbed)'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(calc_prod_ratio, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'calc_prod_ratio  = ', calc_prod_ratio, ' on task ', task_id
        
        ! res_het
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'res_het  = ', res_het, ' (unperturbed)'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(res_het, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'res_het  = ', res_het, ' on task ', task_id
        
        ! res_zoo2
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'res_zoo2  = ', res_zoo2, ' (unperturbed)'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(res_zoo2, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'res_zoo2  = ', res_zoo2, ' on task ', task_id
        
        ! biosynth
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'biosynth  = ', biosynth, ' (unperturbed)'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(biosynth, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'biosynth  = ', biosynth, ' on task ', task_id
        
        ! calc_diss_rate
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'calc_diss_rate  = ', calc_diss_rate, ' (unperturbed)'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(calc_diss_rate, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'calc_diss_rate  = ', calc_diss_rate, ' on task ', task_id
        
        ! calc_diss_rate2
        IF (mype_model==0 .and. task_id==1) WRITE(*,*) 'calc_diss_rate2  = ', calc_diss_rate2, ' (unperturbed)'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(calc_diss_rate2, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,*) 'calc_diss_rate2  = ', calc_diss_rate2, ' on task ', task_id

ENDIF


! *****************************************************
! *** Call PDAF initialization routine on all PEs.  ***
! ***                                               ***
! *** For all filters, first the arrays of integer  ***
! *** and real number parameters are initialized.   ***
! *** Subsequently, PDAF_init is called.            ***
! *****************************************************

  ! *** All other filters                       ***
  ! *** SEIK, LSEIK, ETKF, LETKF, ESTKF, LESTKF ***
  filter_param_i(1) = dim_state_p ! State dimension
  filter_param_i(2) = dim_ens     ! Size of ensemble
  filter_param_i(3) = 0           ! Smoother lag (not implemented here)
  filter_param_i(4) = incremental ! Whether to perform incremental analysis
  filter_param_i(5) = type_forget ! Type of forgetting factor
  filter_param_i(6) = type_trans  ! Type of ensemble transformation
  filter_param_i(7) = type_sqrt   ! Type of transform square-root (SEIK-sub4/ESTKF)
  filter_param_r(1) = forget      ! Forgetting factor
     
  CALL PDAF_init(filtertype, subtype, step_null, &
       filter_param_i, 7,&
       filter_param_r, 2, &
       COMM_model, COMM_filter, COMM_couple, &
       task_id, n_modeltasks, filterpe, init_ens_pdaf, &
       screen, status_pdaf)


! *** Check whether initialization of PDAF was successful ***
  IF (status_pdaf /= 0) THEN
     WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
          'ERROR ', status_pdaf, &
          ' in initialization of PDAF - stopping! (PE ', mype_world,')'
     CALL abort_parallel()
  END IF

  
! ******************************
! *** Initialize file output ***
! ******************************

  writepe = .FALSE.
  IF (filterpe) THEN
     IF (mype_filter==0) writepe = .TRUE.
  ENDIF
  ! Initialize PDAF netCDF output file:
  !    - at beginning of each new year
  !    - at start of assimililation experiment
  IF ((.not. this_is_pdaf_restart) .or. (daynew==1)) THEN
    CALL netCDF_init('mean')
    IF (write_ens) CALL netCDF_init('memb')
    CALL netCDF_STD_init()
    CALL init_carbonfluxes_diags_out()
  ENDIF


! ******************************'***
! *** Prepare ensemble forecasts ***
! ******************************'***

  IF (mype_submodel==0) THEN
     WRITE (*,'(1x,a,i5)') 'FESOM-PDAF: INITIALIZE PDAF before barrier, task: ', task_id
  END IF


  call timeit(6, 'new')
  CALL MPI_BARRIER(MPI_COMM_WORLD, MPIerr)
  call timeit(6, 'old')

  CALL PDAF_get_state(steps, timenow, doexit, next_observation_pdaf, &
       distribute_state_pdaf, prepoststep_pdaf, status_pdaf)


! ***********************************************************************************
! *** Allocate arrays for effective observation dimension and localization radius ***
! ***********************************************************************************
  
  ALLOCATE(eff_dim_obs(mydim_nod2d))
  ALLOCATE(loc_radius(mydim_nod2d))

! **************************************************
! *** Initialize file for synthetic observations ***
! **************************************************

  IF (filtertype==100 .and. mype_world==0) THEN
        WRITE(*,*) "Synthetic observations not yet implemented in this version - stopping!"
        CALL abort_parallel
  END IF


! ***************************************
! *** Get domain limiting coordinates ***
! ***************************************

!~   IF (filterpe) CALL ignore_nod_pdaf() ! Seems to cause problems in FESOM2.0 (SigSegV)
                                          ! Not sure if needed

    CALL PDAFomi_get_domain_limits_unstr(myDim_nod2d, mesh_fesom%geo_coord_nod2D)

! ***********************************
! **** Atmospheric stochasticity  ***
! ***********************************
  
  IF (atmos_stochasticity_ON) THEN
    ! initialize atmospheric stochasticity at (re)start
    call init_atmos_ens_stochasticity()
    ! create stochasticity file
    !    - at start of assimilation experiment
    !    - at beginning of every new year
    IF ((.not. this_is_pdaf_restart) .or. (yearnew .ne. yearold)) THEN
       call init_atmos_stochasticity_output()
    ENDIF
  ENDIF
  
  ! in case of restart, reset forget to 0.99 or 1.00
  ! note: at restarts, forgetting factor is saved and read with atmospheric stochasticity
  IF (this_is_pdaf_restart) THEN
    CALL PDAF_reset_forget(forget)
  ENDIF
    
! *****************
! *** Debugging ***
! *****************

!   At this node and depth, we see a problem in the model output:
    debug_id_depth = 1
    debug_id_nod2  = 85615 !96487 !3833
    ens_member_debug = 0


END SUBROUTINE init_pdaf
