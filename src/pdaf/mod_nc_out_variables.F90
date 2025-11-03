MODULE mod_nc_out_variables

! USES:
USE mod_assim_pdaf, &
   ONLY: id, nfields, assimilateBGC, assimilatePHY, &
         phymin, phymax, bgcmin, bgcmax, &
         cda_phy, cda_bio
USE mod_parallel_pdaf, &
   ONLY: writepe, mype_world
USE obs_chl_cci_pdafomi,    ONLY: assim_o_chl_cci
USE obs_DIC_glodap_pdafomi, ONLY: assim_o_DIC_glodap
USE obs_Alk_glodap_pdafomi, ONLY: assim_o_Alk_glodap
USE obs_pCO2_SOCAT_pdafomi, ONLY: assim_o_pCO2_SOCAT
USE obs_o2_comf_pdafomi,    ONLY: assim_o_o2_comf
USE obs_n_comf_pdafomi,     ONLY: assim_o_n_comf
USE obs_o2_argo_pdafomi,    ONLY: assim_o_o2_argo
USE obs_n_argo_pdafomi,     ONLY: assim_o_n_argo
USE obs_o2_merged_pdafomi,  ONLY: assim_o_o2_merged
USE obs_sss_smos_pdafomi,   ONLY: assim_o_sss
USE obs_sss_cci_pdafomi,    ONLY: assim_o_sss_cci
USE obs_ssh_cmems_pdafomi,  ONLY: assim_o_ssh 
USE obs_sst_pdafomi,        ONLY: assim_o_sst 
USE obs_TSprof_EN4_pdafomi, ONLY: assim_o_en4_s
USE obs_TSprof_EN4_pdafomi, ONLY: assim_o_en4_t

IMPLICIT NONE

character(len=200) :: filename_std = ''       ! Full name of output file

LOGICAL :: w_daymemb    = .false.       ! whether to write any daily ensemble member states
LOGICAL :: w_dayensm    = .false.       ! whether to write any daily ensemble mean states
LOGICAL :: w_monmemb    = .false.       ! whether to write any monthly ensemble member states
LOGICAL :: w_monensm    = .false.       ! whether to write any monthly ensemble mean states
LOGICAL :: w_mm         = .false.       ! whether to write any m-fields (day-averages)
LOGICAL :: w_sm         = .false.       ! whether to write any m-fields (day-averages) of standard deviation

! Field description:

integer :: id_var                         ! Index of a variable in state vector

type state_field
   integer :: ndims = 0                   ! Number of field dimensions (1 or 2)
   logical :: nz1 = .true.                ! Vertical coordinates (on levels / on layers)
   character(len=10) :: variable = ''     ! Name of field
   character(len=50) :: long_name = ''    ! Long name of field
   character(len=20) :: units = ''        ! Unit of variable
   integer :: varid(9)                    ! To write to netCDF file
   logical :: updated = .true.            ! Whether variable is updated through assimilation
   logical :: output(8,3) = .false.       ! How frequently output is written
   logical :: bgc = .false.               ! Whether variable is biogeochemistry (or physics)
   integer :: id_state                    ! Field index in full state vector
   integer :: id_dim                      ! Field index in list of 2D/3D-fields
   integer :: id_type                     ! Field index in list of phy/bgc-fields
   integer :: trnumfesom = -999           ! Tracer index in FESOM-REcoM
   integer :: tridfesom = -999            ! Tracer ID in FESOM-REcoM
   integer :: id_tr = -999                ! Field index in list of 3D model tracer fields
end type state_field

type(state_field), allocatable :: sfields(:) ! Type variable holding the
                                             ! definitions of model fields
                                             
INTEGER :: nfields_3D                     ! Number of 3D fields in state vector
INTEGER :: nfields_2D                     ! """       2D fields """
INTEGER :: nfields_phy                    ! """       physics fields """
INTEGER :: nfields_bgc                    ! """       biogeochem. fields """
INTEGER :: nfields_tr3D                   ! """       3D model tracer """

INTEGER, ALLOCATABLE :: ids_3D(:)         ! List of 3D-field IDs
INTEGER, ALLOCATABLE :: ids_2D(:)         ! """       2D fields """
INTEGER, ALLOCATABLE :: ids_phy(:)        ! """       physics fields """
INTEGER, ALLOCATABLE :: ids_bgc(:)        ! """       biogeochem. fields """
INTEGER, ALLOCATABLE :: ids_tr3D(:)       ! """       3D model tracer """

! Indeces of whether to write: 
INTEGER, PARAMETER :: ff=1, aa=2, mm=3, ii=4 ! forecast (ff), analysis (aa), mean (mm) and initial (ii)
INTEGER, PARAMETER :: sf=5, sa=6, si=7, sm=8 ! ensemble standard deviation snapshots: forecast (sf), analysis (sa), initial (si) and mean (sm)
INTEGER, PARAMETER :: oo=1, ee=2, dd=3       ! any output (oo), ensemble members (ee) and daily values (dd)

LOGICAL :: setoutput(18)
                                             
CONTAINS

SUBROUTINE init_sfields()

! Local variables:
CHARACTER(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file
CHARACTER(len=12)  :: outputmessage(8,3)

LOGICAL            :: upd_ssh , &           ! physics
                      upd_u   , &
                      upd_v   , &
                      upd_w   , &
                      upd_temp, &
                      upd_salt, &
                      upd_ice , &
                      upd_MLD1, &           ! physics diagnostics
                      upd_MLD2, &
                      upd_PhyChl   , &      ! chlorophyll
                      upd_DiaChl   , &
                      upd_DIC      , &      ! dissolved tracers
                      upd_DOC      , &
                      upd_Alk      , &
                      upd_DIN      , &
                      upd_DON      , &
                      upd_O2       , &
                      upd_pCO2s    , &      ! surface carbon diags
                      upd_CO2f     , &
                      upd_alphaCO2 , &
                      upd_PistonVel, &
                      upd_DiaN     , &      ! diatoms
                      upd_DiaC     , &
                      upd_DiaSi    , &
                      upd_PhyCalc  , &      ! small phyto
                      upd_PhyC     , &
                      upd_PhyN     , &
                      upd_Zo1C     , &      ! het
                      upd_Zo1N     , &
                      upd_Zo2C     , &      ! zoo 2
                      upd_Zo2N     , &
                      upd_DetC     , &      ! small det
                      upd_DetCalc  , &
                      upd_DetSi    , &
                      upd_DetN     , &
                      upd_Det2C    , &      ! large det
                      upd_Det2N    , &
                      upd_Det2Si   , &
                      upd_Det2Calc , &
                      upd_export   , &      ! diags
                      upd_PAR      , &
                      upd_NPPn     , &
                      upd_NPPd     , &
                      upd_sigma
                      
INTEGER            :: b,p,s,j ! counters

ALLOCATE(sfields(nfields))

! ****************
! *** Physics ****
! ****************

! SSH
sfields(id% ssh) % ndims = 1
sfields(id% ssh) % variable = 'SSH'
sfields(id% ssh) % long_name = 'Sea surface height'
sfields(id% ssh) % units = 'm'
sfields(id% ssh) % updated = .true.
sfields(id% ssh) % bgc = .false.

! u
sfields(id% u) % ndims = 2
sfields(id% u) % nz1 = .true.
sfields(id% u) % variable = 'u'
sfields(id% u) % long_name = 'Zonal velocity (interpolated on nodes)'
sfields(id% u) % units = 'm/s'
sfields(id% u) % updated = .true.
sfields(id% u) % bgc = .false.

! v
sfields(id% v) % ndims = 2
sfields(id% v) % nz1 = .true.
sfields(id% v) % variable = 'v'
sfields(id% v) % long_name = 'Meridional velocity (interpolated on nodes)'
sfields(id% v) % units = 'm/s'
sfields(id% v) % updated = .true.
sfields(id% v) % bgc = .false.

! w
sfields(id% w) % ndims = 2
sfields(id% w) % nz1 = .false.
sfields(id% w) % variable = 'w'
sfields(id% w) % long_name = 'Vertical velocity'
sfields(id% w) % units = 'm/s'
sfields(id% w) % updated = .false.
sfields(id% w) % bgc = .false.


! temp
sfields(id% temp) % ndims = 2
sfields(id% temp) % nz1 = .true.
sfields(id% temp) % variable = 'T'
sfields(id% temp) % long_name = 'Temperature'
sfields(id% temp) % units = 'degC'
sfields(id% temp) % updated = .true.
sfields(id% temp) % bgc = .false.


! salt
sfields(id% salt) % ndims = 2
sfields(id% salt) % nz1 = .true.
sfields(id% salt) % variable = 'S'
sfields(id% salt) % long_name = 'Salinity'
sfields(id% salt) % units = 'psu'
sfields(id% salt) % updated = .true.
sfields(id% salt) % bgc = .false.


! **********************
! *** Physics Diags ****
! **********************

! ice
sfields(id% a_ice) % ndims = 1
sfields(id% a_ice) % variable = 'ice'
sfields(id% a_ice) % long_name = 'Sea-ice concentration'
sfields(id% a_ice) % units = '1'
sfields(id% a_ice) % updated = .false.
sfields(id% a_ice) % bgc = .false.

! MLD1
sfields(id% MLD1) % ndims = 1
sfields(id% MLD1) % variable = 'MLD1'
sfields(id% MLD1) % long_name = 'Boundary layer depth (Large et al. 1997)'
sfields(id% MLD1) % units = 'm'
sfields(id% MLD1) % updated = .false.
sfields(id% MLD1) % bgc = .false.

! MLD2
sfields(id% MLD2) % ndims = 1
sfields(id% MLD2) % variable = 'MLD2'
sfields(id% MLD2) % long_name = 'Mixed layer depth (sigma 0.03; Boyer-Montegut et al. 2004)'
sfields(id% MLD2) % units = 'm'
sfields(id% MLD2) % updated = .false.
sfields(id% MLD2) % bgc = .false.

! **********************
! *** Chlorophyll   ****
! **********************

! chlorophyll small phytoplankton
sfields(id% PhyChl) % ndims = 2
sfields(id% PhyChl) % nz1 = .true.
sfields(id% PhyChl) % variable = 'PhyChl'
sfields(id% PhyChl) % long_name = 'Chlorophyll-a small phytoplankton'
sfields(id% PhyChl) % units = 'mg chl m-3'
sfields(id% PhyChl) % updated = .false.
sfields(id% PhyChl) % bgc = .true.

! chlorophyll diatoms
sfields(id% DiaChl) % ndims = 2
sfields(id% DiaChl) % nz1 = .true.
sfields(id% DiaChl) % variable = 'DiaChl'
sfields(id% DiaChl) % long_name = 'Chlorophyll-a diatoms'
sfields(id% DiaChl) % units = 'mg chl m-3'
sfields(id% DiaChl) % updated = .false.
sfields(id% DiaChl) % bgc = .true.

! **************************
! *** Dissolved tracers ****
! **************************

! DIC
sfields(id% DIC) % ndims = 2
sfields(id% DIC) % nz1 = .true.
sfields(id% DIC) % variable = 'DIC'
sfields(id% DIC) % long_name = 'Dissolved inorganic carbon'
sfields(id% DIC) % units = 'mmol C m-3'
sfields(id% DIC) % updated = .false.
sfields(id% DIC) % bgc = .true.

! DOC
sfields(id% DOC) % ndims = 2
sfields(id% DOC) % nz1 = .true.
sfields(id% DOC) % variable = 'DOC'
sfields(id% DOC) % long_name = 'Dissolved organic carbon'
sfields(id% DOC) % units = 'mmol C m-3'
sfields(id% DOC) % updated = .false.
sfields(id% DOC) % bgc = .true.

! Alkalinity
sfields(id% Alk) % ndims = 2
sfields(id% Alk) % nz1 = .true.
sfields(id% Alk) % variable = 'Alk'
sfields(id% Alk) % long_name = 'Alkalinity'
sfields(id% Alk) % units = 'mmol m-3'
sfields(id% Alk) % updated = .false.
sfields(id% Alk) % bgc = .true.

! DIN
sfields(id% DIN) % ndims = 2
sfields(id% DIN) % nz1 = .true.
sfields(id% DIN) % variable = 'DIN'
sfields(id% DIN) % long_name = 'Dissolved inorganic nitrogen'
sfields(id% DIN) % units = 'mmol m-3'
sfields(id% DIN) % updated = .false.
sfields(id% DIN) % bgc = .true.

! DON
sfields(id% DON) % ndims = 2
sfields(id% DON) % nz1 = .true.
sfields(id% DON) % variable = 'DON'
sfields(id% DON) % long_name = 'Dissolved organic nitrogen'
sfields(id% DON) % units = 'mmol m-3'
sfields(id% DON) % updated = .false.
sfields(id% DON) % bgc = .true.

! Oxygen
sfields(id% O2) % ndims = 2
sfields(id% O2) % nz1 = .true.
sfields(id% O2) % variable = 'O2'
sfields(id% O2) % long_name = 'Oxygen'
sfields(id% O2) % units = 'mmol m-3'
sfields(id% O2) % updated = .false.
sfields(id% O2) % bgc = .true.

! *****************************
! *** Surface Carbon Diags ****
! *****************************

! pCO2
sfields(id% pCO2s) % ndims = 1
sfields(id% pCO2s) % variable = 'pCO2s'
sfields(id% pCO2s) % long_name = 'Partial pressure CO2 surface ocean'
sfields(id% pCO2s) % units = 'micro atm'
sfields(id% pCO2s) % updated = .false.
sfields(id% pCO2s) % bgc = .true.

! CO2f
sfields(id% CO2f) % ndims = 1
sfields(id% CO2f) % variable = 'CO2f'
sfields(id% CO2f) % long_name = 'CO2 flux from atmosphere into ocean'
sfields(id% CO2f) % units = 'mmol C m-2 d-1'
sfields(id% CO2f) % updated = .false.
sfields(id% CO2f) % bgc = .true.

! Solubility of CO2
sfields(id% alphaCO2) % ndims = 1
sfields(id% alphaCO2) % variable = 'alphaCO2'
sfields(id% alphaCO2) % long_name = 'solubility of surface CO2'
sfields(id% alphaCO2) % units = 'mol kg-1 atm-1'
sfields(id% alphaCO2) % updated = .false.
sfields(id% alphaCO2) % bgc = .true.

! Piston velocity
sfields(id% PistonVel) % ndims = 1
sfields(id% PistonVel) % variable = 'Kw660'
sfields(id% PistonVel) % long_name = 'air-sea piston velocity'
sfields(id% PistonVel) % units = 'm/s'
sfields(id% PistonVel) % updated = .false.
sfields(id% PistonVel) % bgc = .true.

! *****************************
! *** Small Phyto          ****
! *****************************

! PhyN
sfields(id% PhyN) % ndims = 2
sfields(id% PhyN) % variable = 'PhyN'
sfields(id% PhyN) % long_name = 'intracell nitrogen small phytoplankton'
sfields(id% PhyN) % units = 'mmol m-3'
sfields(id% PhyN) % updated = .false.
sfields(id% PhyN) % bgc = .true.

! PhyC
sfields(id% PhyC) % ndims = 2
sfields(id% PhyC) % variable = 'PhyC'
sfields(id% PhyC) % long_name = 'intracell carbon small phytoplankton'
sfields(id% PhyC) % units = 'mmol C m-3'
sfields(id% PhyC) % updated = .false.
sfields(id% PhyC) % bgc = .true.

! PhyCalc
sfields(id% PhyCalc) % ndims = 2
sfields(id% PhyCalc) % variable = 'PhyCalc'
sfields(id% PhyCalc) % long_name = 'calcium carbonate small phytoplankton'
sfields(id% PhyCalc) % units = 'mmol m-3'
sfields(id% PhyCalc) % updated = .false.
sfields(id% PhyCalc) % bgc = .true.

! *****************************
! *** diatoms              ****
! *****************************

! DiaN
sfields(id% DiaN) % ndims = 2
sfields(id% DiaN) % variable = 'DiaN'
sfields(id% DiaN) % long_name = 'intracell nitrogen diatoms'
sfields(id% DiaN) % units = 'mmol m-3'
sfields(id% DiaN) % updated = .false.
sfields(id% DiaN) % bgc = .true.

! DiaC
sfields(id% DiaC) % ndims = 2
sfields(id% DiaC) % variable = 'DiaC'
sfields(id% DiaC) % long_name = 'intracell carbon diatom'
sfields(id% DiaC) % units = 'mmol C m-3'
sfields(id% DiaC) % updated = .false.
sfields(id% DiaC) % bgc = .true.

! DiaSi
sfields(id% DiaSi) % ndims = 2
sfields(id% DiaSi) % variable = 'DiaSi'
sfields(id% DiaSi) % long_name = 'intracell Si diatom'
sfields(id% DiaSi) % units = 'mmol m-3'
sfields(id% DiaSi) % updated = .false.
sfields(id% DiaSi) % bgc = .true.

! *****************************
! *** zooplankton          ****
! *****************************

! Zo1C
sfields(id% Zo1C) % ndims = 2
sfields(id% Zo1C) % variable = 'Zo1C'
sfields(id% Zo1C) % long_name = 'carbon in small zooplankton'
sfields(id% Zo1C) % units = 'mmol C m-3'
sfields(id% Zo1C) % updated = .false.
sfields(id% Zo1C) % bgc = .true.

! Zo1N
sfields(id% Zo1N) % ndims = 2
sfields(id% Zo1N) % variable = 'Zo1N'
sfields(id% Zo1N) % long_name = 'nitrogen in small zooplankton'
sfields(id% Zo1N) % units = 'mmol C m-3'
sfields(id% Zo1N) % updated = .false.
sfields(id% Zo1N) % bgc = .true.

! Zo2C
sfields(id% Zo2C) % ndims = 2
sfields(id% Zo2C) % variable = 'Zo2C'
sfields(id% Zo2C) % long_name = 'carbon in macrozooplankton'
sfields(id% Zo2C) % units = 'mmol C m-3'
sfields(id% Zo2C) % updated = .false.
sfields(id% Zo2C) % bgc = .true.

! Zo2N
sfields(id% Zo2N) % ndims = 2
sfields(id% Zo2N) % variable = 'Zo2N'
sfields(id% Zo2N) % long_name = 'nitrogen in macrozooplankton'
sfields(id% Zo2N) % units = 'mmol C m-3'
sfields(id% Zo2N) % updated = .false.
sfields(id% Zo2N) % bgc = .true.

! *****************************
! *** detritus             ****
! *****************************

! DetC
sfields(id% DetC) % ndims = 2
sfields(id% DetC) % variable = 'DetC'
sfields(id% DetC) % long_name = 'carbon in small detritus'
sfields(id% DetC) % units = 'mmol C m-3'
sfields(id% DetC) % updated = .false.
sfields(id% DetC) % bgc = .true.

! DetCalc
sfields(id% DetCalc) % ndims = 2
sfields(id% DetCalc) % variable = 'DetCalc'
sfields(id% DetCalc) % long_name = 'calcite in small detritus'
sfields(id% DetCalc) % units = 'mmol C m-3'
sfields(id% DetCalc) % updated = .false.
sfields(id% DetCalc) % bgc = .true.

! DetN
sfields(id% DetN) % ndims = 2
sfields(id% DetN) % variable = 'DetN'
sfields(id% DetN) % long_name = 'nitrogen in small detritus'
sfields(id% DetN) % units = 'mmol C m-3'
sfields(id% DetN) % updated = .false.
sfields(id% DetN) % bgc = .true.

! DetSi
sfields(id% DetSi) % ndims = 2
sfields(id% DetSi) % variable = 'DetSi'
sfields(id% DetSi) % long_name = 'silicate in small detritus'
sfields(id% DetSi) % units = 'mmol C m-3'
sfields(id% DetSi) % updated = .false.
sfields(id% DetSi) % bgc = .true.

! Det2 C
sfields(id% Det2C) % ndims = 2
sfields(id% Det2C) % variable = 'Det2C'
sfields(id% Det2C) % long_name = 'carbon in large detritus'
sfields(id% Det2C) % units = 'mmol C m-3'
sfields(id% Det2C) % updated = .false.
sfields(id% Det2C) % bgc = .true.

! Det2 Calc
sfields(id% Det2Calc) % ndims = 2
sfields(id% Det2Calc) % variable = 'Det2Calc'
sfields(id% Det2Calc) % long_name = 'calcite in large detritus'
sfields(id% Det2Calc) % units = 'mmol C m-3'
sfields(id% Det2Calc) % updated = .false.
sfields(id% Det2Calc) % bgc = .true.

! Det2 N
sfields(id% Det2N) % ndims = 2
sfields(id% Det2N) % variable = 'Det2N'
sfields(id% Det2N) % long_name = 'nitrogen in large detritus'
sfields(id% Det2N) % units = 'mmol C m-3'
sfields(id% Det2N) % updated = .false.
sfields(id% Det2N) % bgc = .true.

! Det2 Si
sfields(id% Det2Si) % ndims = 2
sfields(id% Det2Si) % variable = 'Det2Si'
sfields(id% Det2Si) % long_name = 'silicate in large detritus'
sfields(id% Det2Si) % units = 'mmol C m-3'
sfields(id% Det2Si) % updated = .false.
sfields(id% Det2Si) % bgc = .true.

! *****************************
! *** diagnostics          ****
! *****************************

! PAR
sfields(id% PAR) % ndims = 2
sfields(id% PAR) % variable = 'PAR'
sfields(id% PAR) % long_name = 'photosynthetically active radiation'
sfields(id% PAR) % units = 'W m-2'
sfields(id% PAR) % updated = .false.
sfields(id% PAR) % bgc = .true.

! NPPn
sfields(id% NPPn) % ndims = 2
sfields(id% NPPn) % variable = 'NPPn'
sfields(id% NPPn) % long_name = 'mean net primary production small phytoplankton'
sfields(id% NPPn) % units = 'mmol C m-2 d-1'
sfields(id% NPPn) % updated = .false.
sfields(id% NPPn) % bgc = .true.

! NPPd
sfields(id% NPPd) % ndims = 2
sfields(id% NPPd) % variable = 'NPPd'
sfields(id% NPPd) % long_name = 'mean net primary production diatoms'
sfields(id% NPPd) % units = 'mmol C m-2 d-1'
sfields(id% NPPd) % updated = .false.
sfields(id% NPPd) % bgc = .true.

! Export production
sfields(id% export) % ndims = 1
sfields(id% export) % variable = 'export'
sfields(id% export) % long_name = 'export through particle sinking at 190m'
sfields(id% export) % units = 'mmol m-2 day-1'
sfields(id% export) % updated = .false.
sfields(id% export) % bgc = .true.

! Potential density
sfields(id% sigma) % ndims = 2
sfields(id% sigma) % variable = 'sigma'
sfields(id% sigma) % long_name = 'potential density'
sfields(id% sigma) % units = 'kg liter-1'
sfields(id% sigma) % updated = .false.
sfields(id% sigma) % bgc = .false.



!~ ! TChl
!~ sfields(id% TChl) % ndims = 2
!~ sfields(id% TChl) % variable = 'TChl'
!~ sfields(id% TChl) % long_name = 'Total chlorophyll (PhyChl+DiaChl)'
!~ sfields(id% TChl) % units = 'mg chl m-3'
!~ sfields(id% TChl) % updated = .false.
!~ sfields(id% TChl) % bgc = .true.

!~ ! TDN
!~ sfields(id% TDN) % ndims = 2
!~ sfields(id% TDN) % variable = 'TDN'
!~ sfields(id% TDN) % long_name = 'Total dissolved nitrogen (DIN+DON)'
!~ sfields(id% TDN) % units = 'mmol m-3'
!~ sfields(id% TDN) % updated = .false.
!~ sfields(id% TDN) % bgc = .true.

!~ ! TOC
!~ sfields(id% TOC) % ndims = 2
!~ sfields(id% TOC) % variable = 'TOC'
!~ sfields(id% TOC) % long_name = 'Total Organic Carbon (PhyC+DiaC+DetC+DOC+HetC)'
!~ sfields(id% TOC) % units = 'mmol C m-3'
!~ sfields(id% TOC) % updated = .false.
!~ sfields(id% TOC) % bgc = .true.

! ************************************************
! ***   Tracer index and ID from FESOM-REcoM   ***
! ************************************************

 sfields(id% temp    ) % trnumfesom = 1  ! temperature
 sfields(id% salt    ) % trnumfesom = 2  ! salinity 
 sfields(id% DIN     ) % trnumfesom = 3  ! DIN
 sfields(id% DIC     ) % trnumfesom = 4  ! DIC
 sfields(id% Alk     ) % trnumfesom = 5  ! Alk
 sfields(id% PhyN    ) % trnumfesom = 6  ! PhyN 
 sfields(id% PhyC    ) % trnumfesom = 7  ! PhyC   
 sfields(id% PhyChl  ) % trnumfesom = 8  ! PhyChl
 sfields(id% DetN    ) % trnumfesom = 9  ! DetN
 sfields(id% DetC    ) % trnumfesom = 10 ! DetC
 sfields(id% Zo1N    ) % trnumfesom = 11 ! HetN
 sfields(id% Zo1C    ) % trnumfesom = 12 ! HetC
 sfields(id% DON     ) % trnumfesom = 13 ! DON
 sfields(id% DOC     ) % trnumfesom = 14 ! DOC
 sfields(id% DiaN    ) % trnumfesom = 15 ! DiaN
 sfields(id% DiaC    ) % trnumfesom = 16 ! DiaC
 sfields(id% DiaChl  ) % trnumfesom = 17 ! DiaChl
 sfields(id% DiaSi   ) % trnumfesom = 18 ! DiaSi
 sfields(id% DetSi   ) % trnumfesom = 19 ! DetSi 
!sfields(id%         ) % trnumfesom = 20 ! DSi     
!sfields(id%         ) % trnumfesom = 21 ! Fe
 sfields(id% PhyCalc ) % trnumfesom = 22 ! PhyCalc
 sfields(id% DetCalc ) % trnumfesom = 23 ! DetCalc
 sfields(id% O2      ) % trnumfesom = 24 ! Oxy
 sfields(id% Zo2N    ) % trnumfesom = 25 ! Zoo2N
 sfields(id% Zo2C    ) % trnumfesom = 26 ! Zoo2C
 sfields(id% Det2N   ) % trnumfesom = 27 ! DetZ2N                              
 sfields(id% Det2C   ) % trnumfesom = 28 ! DetZ2C                                    
 sfields(id% Det2Si  ) % trnumfesom = 29 ! DetZ2Si                            
 sfields(id% Det2Calc) % trnumfesom = 30 ! DetZ2Calc
 
 sfields(id% temp    )% tridfesom =    0 ! temperature
 sfields(id% salt    )% tridfesom =    1 ! salinity
 sfields(id% DIN     )% tridfesom = 1001 ! DIN
 sfields(id% DIC     )% tridfesom = 1002 ! DIC
 sfields(id% Alk     )% tridfesom = 1003 ! Alk
 sfields(id% PhyN    )% tridfesom = 1004 ! PhyN 
 sfields(id% PhyC    )% tridfesom = 1005 ! PhyC   
 sfields(id% PhyChl  )% tridfesom = 1006 ! PhyChl
 sfields(id% DetN    )% tridfesom = 1007 ! DetN
 sfields(id% DetC    )% tridfesom = 1008 ! DetC
 sfields(id% Zo1N    )% tridfesom = 1009 ! HetN
 sfields(id% Zo1C    )% tridfesom = 1010 ! HetC
 sfields(id% DON     )% tridfesom = 1011 ! DON
 sfields(id% DOC     )% tridfesom = 1012 ! DOC
 sfields(id% DiaN    )% tridfesom = 1013 ! DiaN
 sfields(id% DiaC    )% tridfesom = 1014 ! DiaC
 sfields(id% DiaChl  )% tridfesom = 1015 ! DiaChl
 sfields(id% DiaSi   )% tridfesom = 1016 ! DiaSi
 sfields(id% DetSi   )% tridfesom = 1017 ! DetSi 
!sfields(id%         )% tridfesom = 1018 ! DSi     
!sfields(id%         )% tridfesom = 1019 ! Fe
 sfields(id% PhyCalc )% tridfesom = 1020 ! PhyCalc
 sfields(id% DetCalc )% tridfesom = 1021 ! DetCalc
 sfields(id% O2      )% tridfesom = 1022 ! Oxy
 sfields(id% Zo2N    )% tridfesom = 1023 ! Zoo2N
 sfields(id% Zo2C    )% tridfesom = 1024 ! Zoo2C
 sfields(id% Det2N   )% tridfesom = 1025 ! DetZ2N                              
 sfields(id% Det2C   )% tridfesom = 1026 ! DetZ2C                                    
 sfields(id% Det2Si  )% tridfesom = 1027 ! DetZ2Si                            
 sfields(id% Det2Calc)% tridfesom = 1028 ! DetZ2Calc

! ************************************************
! ***   Read updated variables from namelist   ***
! ************************************************

! The logical "updated" describes whether a variables is updated in at least one sweep

! In case of diagnostic variables, "updated" is False. Setting diagnostics variables to False and the others to True, is set in namelist.
! In case of weak coupling and only physics or BGC assimilation, "updated" is False for the other type of fields. This is done below.
! "updated" is used in init_dim_l_pdaf: only updated fields are included in local state
! "updated" is used in the output routine: option to write out only updated fields

! *** Read namelist file ***
  IF (mype_world==0) WRITE(*,*) 'Read namelist file for updated variables: ',nmlfile
  
  NAMELIST /updated/ &
     upd_ssh , &           ! physics
     upd_u   , &
     upd_v   , &
     upd_w   , &
     upd_temp, &
     upd_salt, &
     upd_ice , &
     upd_MLD1, &           ! physics diagnostics
     upd_MLD2, &
     upd_PhyChl   , &      ! chlorophyll
     upd_DiaChl   , &
     upd_DIC      , &      ! dissolved tracers
     upd_DOC      , &
     upd_Alk      , &
     upd_DIN      , &
     upd_DON      , &
     upd_O2       , &
     upd_pCO2s    , &      ! surface carbon diags
     upd_CO2f     , &
     upd_alphaCO2 , &
     upd_PistonVel, &
     upd_DiaN     , &      ! diatoms
     upd_DiaC     , &
     upd_DiaSi    , &
     upd_PhyCalc  , &      ! small phyto
     upd_PhyC     , &
     upd_PhyN     , &
     upd_Zo1C     , &      ! het
     upd_Zo1N     , &
     upd_Zo2C     , &      ! zoo 2
     upd_Zo2N     , &
     upd_DetC     , &      ! small det
     upd_DetCalc  , &
     upd_DetSi    , &
     upd_DetN     , &
     upd_Det2C    , &      ! large det
     upd_Det2N    , &
     upd_Det2Si   , &
     upd_Det2Calc , &
     upd_export   , &      ! diags
     upd_PAR      , &
     upd_NPPn     , &
     upd_NPPd     , &
     upd_sigma
     

  OPEN  (20,file=nmlfile)
  READ  (20,NML=updated)
  CLOSE (20)
  
  sfields(id% ssh      ) % updated = upd_ssh
  sfields(id% u        ) % updated = upd_u
  sfields(id% v        ) % updated = upd_v
  sfields(id% w        ) % updated = upd_w
  sfields(id% temp     ) % updated = upd_temp
  sfields(id% salt     ) % updated = upd_salt
  sfields(id% a_ice    ) % updated = upd_ice
  sfields(id% MLD1     ) % updated = upd_MLD1
  sfields(id% MLD2     ) % updated = upd_MLD2
  sfields(id% sigma    ) % updated = upd_sigma
  
  sfields(id% PhyChl   ) % updated = upd_PhyChl
  sfields(id% DiaChl   ) % updated = upd_DiaChl
  
  sfields(id% DIC      ) % updated = upd_DIC
  sfields(id% DOC      ) % updated = upd_DOC
  sfields(id% Alk      ) % updated = upd_Alk
  sfields(id% DIN      ) % updated = upd_DIN
  sfields(id% DON      ) % updated = upd_DON
  sfields(id% O2       ) % updated = upd_O2
  
  sfields(id% pCO2s    ) % updated = upd_pCO2s
  sfields(id% CO2f     ) % updated = upd_CO2f
  sfields(id% alphaCO2 ) % updated = upd_alphaCO2
  sfields(id% PistonVel) % updated = upd_PistonVel  
  
  sfields(id% DiaN     ) % updated = upd_DiaN
  sfields(id% DiaC     ) % updated = upd_DiaC
  sfields(id% DiaSi    ) % updated = upd_DiaSi
  
  sfields(id% PhyCalc  ) % updated = upd_PhyCalc
  sfields(id% PhyC     ) % updated = upd_PhyC
  sfields(id% PhyN     ) % updated = upd_PhyN
  
  sfields(id% Zo1C     ) % updated = upd_Zo1C
  sfields(id% Zo1N     ) % updated = upd_Zo1N
  sfields(id% Zo2C     ) % updated = upd_Zo2C
  sfields(id% Zo2N     ) % updated = upd_Zo2N
  
  sfields(id% DetC     ) % updated = upd_DetC
  sfields(id% DetCalc  ) % updated = upd_DetCalc
  sfields(id% DetSi    ) % updated = upd_DetSi
  sfields(id% DetN     ) % updated = upd_DetN
  
  sfields(id% Det2C     ) % updated = upd_Det2C
  sfields(id% Det2Calc  ) % updated = upd_Det2Calc
  sfields(id% Det2Si    ) % updated = upd_Det2Si
  sfields(id% Det2N     ) % updated = upd_Det2N
  
  sfields(id% PAR      ) % updated = upd_PAR
  sfields(id% NPPn     ) % updated = upd_NPPn
  sfields(id% NPPd     ) % updated = upd_NPPd
  sfields(id% export   ) % updated = upd_export
  
  
  ! Physics not assimilated and coupling weak: No update to physics
  IF ((.not. assimilatePHY) .and. (cda_phy=='weak')) THEN
     do p=phymin,phymax
       sfields(p) % updated = .false.
     enddo
  ENDIF
  ! BGC not assimilated and coupling weak: No update to BGC
  IF ((.not. assimilateBGC) .and. (cda_phy=='weak')) THEN
     do b=bgcmin,bgcmax
       sfields(b) % updated = .false.
     enddo
  ENDIF
  
  ! **************************************
  ! ***  Indeces of by type of field   ***
  ! **************************************
  
  ! ___________
  ! ___2D/3D___
  ! count number of 3D and 2D fields (tracers and diagnostics)
  nfields_3D = 0
  nfields_2D = 0
  DO b=1,nfields
   IF (sfields(b) % ndims == 2) nfields_3D = nfields_3D + 1
   IF (sfields(b) % ndims == 1) nfields_2D = nfields_2D + 1
  ENDDO
  
  ! indeces of 3D and 2D fields
  ALLOCATE(ids_3D(nfields_3D))
  p = 1
  DO b=1,nfields
    IF (sfields(b) % ndims == 2) THEN
      ids_3D(p) = b
      sfields(b) % id_dim = p
      p = p+1
    ENDIF
  ENDDO
  IF (mype_world==0) THEN
     DO p=1,nfields_3D
       WRITE (*,'(a, 10x,3a,1x,7a)') &
          'FESOM-PDAF', '3D fields in state vector: ', sfields(ids_3D(p)) % variable
     ENDDO
  ENDIF
  
  ALLOCATE(ids_2D(nfields_2D))
  p = 1
  DO b=1,nfields
    IF (sfields(b) % ndims == 1) THEN
      ids_2D(p) = b
      sfields(b) % id_dim = p
      p = p+1
    ENDIF
  ENDDO
  IF (mype_world==0) THEN
     DO p=1,nfields_2D
       WRITE (*,'(a, 10x,3a,1x,7a)') &
          'FESOM-PDAF', '2D fields in state vector: ', sfields(ids_2D(p)) % variable
     ENDDO
  ENDIF
  
  ! ________________________
  ! ___ model 3D tracers ___
  ! count number of 3D model tracer fields in state vector (only tracers)
  nfields_tr3D = 0
  DO b=1,nfields
   IF (sfields(b)% trnumfesom > 0) nfields_tr3D = nfields_tr3D + 1
  ENDDO
  
  ! indeces of 3D model tracer fields in state vector
  ALLOCATE(ids_tr3D(nfields_tr3D))
  p = 1
  DO b=1,nfields
    IF (sfields(b)% trnumfesom > 0) THEN
      ids_tr3D(p)=b
      sfields(b) % id_tr = p
      p = p+1
    ENDIF
  ENDDO
  IF (mype_world==0) THEN
     DO p=1,nfields_tr3D
       WRITE (*,'(a, 10x,3a,1x,7a)') &
          'FESOM-PDAF', '3D model tracer fields in state vector: ', sfields(ids_tr3D(p)) % variable
     ENDDO
  ENDIF
  
  ! ________________________
  ! ___ phy/bgc fields  ____
  ! count number of phy/bgc fields
  nfields_phy = 0
  nfields_bgc = 0
  DO b=1,nfields
   IF (      sfields(b) % bgc) nfields_bgc = nfields_bgc + 1
   IF (.not. sfields(b) % bgc) nfields_phy = nfields_phy + 1
  ENDDO
  
  ALLOCATE(ids_bgc(nfields_bgc))
  p = 1
  DO b=1,nfields
    IF (sfields(b) % bgc) THEN
      ids_bgc(p) = b
      sfields(b) % id_type = p
      p = p+1
    ENDIF
  ENDDO
  IF (mype_world==0) THEN
     DO p=1,nfields_bgc
       WRITE (*,'(a, 10x,3a,1x,7a)') &
          'FESOM-PDAF', 'bgc fields in state vector: ', sfields(ids_bgc(p)) % variable
     ENDDO
  ENDIF
  
  ALLOCATE(ids_phy(nfields_phy))
  p = 1
  DO b=1,nfields
    IF (.not. sfields(b) % bgc) THEN
      ids_phy(p) = b
      sfields(b) % id_type = p
      p = p+1
    ENDIF
  ENDDO
  IF (mype_world==0) THEN
     DO p=1,nfields_phy
       WRITE (*,'(a, 10x,3a,1x,7a)') &
          'FESOM-PDAF', 'phy fields in state vector: ', sfields(ids_phy(p)) % variable
     ENDDO
  ENDIF  
  
  
! *********************************************************
! ***   Set field-specific output type and frequency    ***
! *********************************************************

! [oo] True   - write output
!      False  - no output
! [ee] True   - write ensemble members
!      False  - write ensemble mean
! [dd] True   - daily
!      False  - monthly

! Defaults: False
! activate one (or multiple if not contradictory - think!) of the following predefined output schemes:

! ___________________________________________________________
! ___ write daily forecast and analysis ensemble members  ___
IF (setoutput(1)) THEN
  DO s=1, nfields
    ! forecast
    sfields(s)% output(ff,oo) = .True.
    sfields(s)% output(ff,ee) = .True.
    sfields(s)% output(ff,dd) = .True.
    ! analysis
    sfields(s)% output(aa,oo) = .True.
    sfields(s)% output(aa,ee) = .True.
    sfields(s)% output(aa,dd) = .True.
  ENDDO
ENDIF

! ________________________________________________________
! ___ write daily forecast and analysis ensemble mean  ___
IF (setoutput(2)) THEN
  DO s=1, nfields
    ! forecast
    sfields(s)% output(ff,oo) = .True.
    sfields(s)% output(ff,dd) = .True.
    ! analysis
    sfields(s)% output(aa,oo) = .True.
    sfields(s)% output(aa,dd) = .True.
  ENDDO
ENDIF

! ______________________________________________________________________________
! ___ write monthly forecast and analysis ensemble mean of updated variables ___
IF (setoutput(3)) THEN
  DO s=1, nfields
    ! forecast
    IF (sfields(s)%updated)   sfields(s)% output(ff,oo) = .True.
    ! analysis
    IF (sfields(s)%updated)   sfields(s)% output(aa,oo) = .True.
  ENDDO
ENDIF

! ___________________________________________
! ___ write daily m-fields ensemble mean  ___
IF (setoutput(4)) THEN
  DO s=1, nfields
    sfields(s)% output(mm,oo) = .True.
    sfields(s)% output(mm,dd) = .True.
  ENDDO
ENDIF

! ___________________________________________
! ___ write initial fields ensemble mean  ___
IF (setoutput(5)) THEN
  DO s=1, nfields
    sfields(s)% output(ii,oo) = .True.
  ENDDO
ENDIF

! ______________________________________________
! ___ write initial fields ensemble members  ___
IF (setoutput(6)) THEN
  DO s=1, nfields
    sfields(s)% output(ii,oo) = .True.
    sfields(s)% output(ii,ee) = .True.
  ENDDO
ENDIF

! ____________________________________________
! ___ write monthly m-fields ensemble mean ___
IF (setoutput(7)) THEN
  DO s=1, nfields
    sfields(s)% output(mm,oo) = .True.
  ENDDO
ENDIF

! ___________________________________________________________________
! ___ write daily m-fields of assimilated variables ensemble mean ___
IF (setoutput(8)) THEN
  ! activate m-field output
  IF (assim_o_sst)          sfields(id% temp)   % output(mm,oo) = .True.
  IF (assim_o_sss)          sfields(id% salt)   % output(mm,oo) = .True.
  IF (assim_o_sss_cci)      sfields(id% salt)   % output(mm,oo) = .True.
  IF (assim_o_en4_t)        sfields(id% temp)   % output(mm,oo) = .True.
  IF (assim_o_en4_s)        sfields(id% salt)   % output(mm,oo) = .True.
  IF (assim_o_ssh)          sfields(id% SSH)    % output(mm,oo) = .True.
  IF (assim_o_chl_cci)      sfields(id% PhyChl) % output(mm,oo) = .True.
  IF (assim_o_chl_cci)      sfields(id% DiaChl) % output(mm,oo) = .True.
  IF (assim_o_DIC_glodap)   sfields(id% DIC)    % output(mm,oo) = .True.
  IF (assim_o_Alk_glodap)   sfields(id% Alk)    % output(mm,oo) = .True.
  IF (assim_o_pCO2_SOCAT)   sfields(id% pCO2s)  % output(mm,oo) = .True.
  IF (assim_o_o2_comf)      sfields(id% O2)     % output(mm,oo) = .True.
  IF (assim_o_n_comf)       sfields(id% DIN)    % output(mm,oo) = .True.
  IF (assim_o_o2_argo)      sfields(id% O2)     % output(mm,oo) = .True.
  IF (assim_o_n_argo)       sfields(id% DIN)    % output(mm,oo) = .True.
  IF (assim_o_o2_merged)    sfields(id% O2)     % output(mm,oo) = .True.
  ! set to daily
  IF (assim_o_sst)          sfields(id% temp)   % output(mm,dd) = .True.
  IF (assim_o_sss)          sfields(id% salt)   % output(mm,dd) = .True.
  IF (assim_o_sss_cci)      sfields(id% salt)   % output(mm,dd) = .True.
  IF (assim_o_en4_t)        sfields(id% temp)   % output(mm,dd) = .True.
  IF (assim_o_en4_s)        sfields(id% salt)   % output(mm,dd) = .True.
  IF (assim_o_ssh)          sfields(id% SSH)    % output(mm,dd) = .True.
  IF (assim_o_chl_cci)      sfields(id% PhyChl) % output(mm,dd) = .True.
  IF (assim_o_chl_cci)      sfields(id% DiaChl) % output(mm,dd) = .True.
  IF (assim_o_DIC_glodap)   sfields(id% DIC)    % output(mm,dd) = .True.
  IF (assim_o_Alk_glodap)   sfields(id% Alk)    % output(mm,dd) = .True.
  IF (assim_o_pCO2_SOCAT)   sfields(id% pCO2s)  % output(mm,dd) = .True.
  IF (assim_o_o2_comf)      sfields(id% O2)     % output(mm,dd) = .True.
  IF (assim_o_n_comf)       sfields(id% DIN)    % output(mm,dd) = .True.
  IF (assim_o_o2_argo)      sfields(id% O2)     % output(mm,dd) = .True.
  IF (assim_o_n_argo)       sfields(id% DIN)    % output(mm,dd) = .True.
  IF (assim_o_o2_merged)    sfields(id% O2)     % output(mm,dd) = .True.
ENDIF

! ___________________________________________________________________________
! ___ write daily m-fields of assimilate-able BGC variables ensemble mean ___
IF (setoutput(9)) THEN
  ! activate m-field output
  sfields(id% PhyChl) % output(mm,oo) = .True.
  sfields(id% DiaChl) % output(mm,oo) = .True.
  sfields(id% DIC)    % output(mm,oo) = .True.
  sfields(id% Alk)    % output(mm,oo) = .True.
  sfields(id% pCO2s)  % output(mm,oo) = .True.
  sfields(id% O2)     % output(mm,oo) = .True.
  sfields(id% DIN)    % output(mm,oo) = .True.
  ! set to daily
  sfields(id% PhyChl) % output(mm,dd) = .True.
  sfields(id% DiaChl) % output(mm,dd) = .True.
  sfields(id% DIC)    % output(mm,dd) = .True.
  sfields(id% Alk)    % output(mm,dd) = .True.
  sfields(id% pCO2s)  % output(mm,dd) = .True.
  sfields(id% O2)     % output(mm,dd) = .True.
  sfields(id% DIN)    % output(mm,dd) = .True.
ENDIF
  
! _______________________________________________________________
! ___ write daily m-fields of CO2 flux and pCO2 ensemble mean ___
IF (setoutput(10)) THEN
  sfields(id% CO2f)  % output(mm,oo) = .True.
  sfields(id% pCO2s) % output(mm,oo) = .True.
  sfields(id% CO2f)  % output(mm,dd) = .True.
  sfields(id% pCO2s) % output(mm,dd) = .True.
ENDIF

! _______________________________________________________________
! ___ write daily m-fields of variables defining the CO2 flux ___
IF (setoutput(11)) THEN
  sfields(id% CO2f     ) % output(mm,oo) = .True.
  sfields(id% pCO2s    ) % output(mm,oo) = .True.
  sfields(id% a_ice    ) % output(mm,oo) = .True.
  sfields(id% alphaCO2 ) % output(mm,oo) = .True.
  sfields(id% PistonVel) % output(mm,oo) = .True.
  
  sfields(id% CO2f     ) % output(mm,dd) = .True.
  sfields(id% pCO2s    ) % output(mm,dd) = .True.
  sfields(id% a_ice    ) % output(mm,dd) = .True.
  sfields(id% alphaCO2 ) % output(mm,dd) = .True.
  sfields(id% PistonVel) % output(mm,dd) = .True.
ENDIF

! ________________________________________________
! ___ write initial fields standard deviation  ___
IF (setoutput(12)) THEN
  DO s=1, nfields
    sfields(s)% output(si,oo) = .True.
  ENDDO
ENDIF

! ____________________________________________
! ___ write daily forecast of standard deviation ___
IF (setoutput(13)) THEN
  DO s=1, nfields
    sfields(s)% output(sf,oo) = .True.
    sfields(s)% output(sf,dd) = .True.
  ENDDO
ENDIF

! ____________________________________________
! ___ write monthly forecast of standard deviation ___
IF (setoutput(14)) THEN
  DO s=1, nfields
    sfields(s)% output(sf,oo) = .True.
  ENDDO
ENDIF

! ____________________________________________
! ___ write monthly forecast and analysis of standard deviation for updated fields ___
IF (setoutput(15)) THEN
  DO s=1, nfields
    IF (sfields(s)%updated)   sfields(s)% output(sf,oo) = .True.
    IF (sfields(s)%updated)   sfields(s)% output(sa,oo) = .True.
  ENDDO
ENDIF

! ____________________________________________
! ___ write monthly mm-fields of standard deviation ___
IF (setoutput(16)) THEN
  DO s=1, nfields
    sfields(s)% output(sm,oo) = .True.
  ENDDO
ENDIF

! ______________________________________________________________
! ___ write daily forecast fields for oxygen and DIN         ___
! (to reconstruct inno_omit)
IF (setoutput(17)) THEN
   ! Oxygen
      sfields(id% O2) % output(ff,oo) = .True. ! activate
      sfields(id% O2) % output(ff,dd) = .True. ! daily
   ! DIN
      sfields(id% DIN) % output(ff,oo) = .True. ! activate
      sfields(id% DIN) % output(ff,dd) = .True. ! daily
ENDIF

! ___________________________________________________________________________
! ___ write daily analysis of assimilate-able BGC variables ensemble mean ___
IF (setoutput(18)) THEN
  ! activate m-field output
  sfields(id% PhyChl) % output(aa,oo) = .True.
  sfields(id% DiaChl) % output(aa,oo) = .True.
  sfields(id% DIC)    % output(aa,oo) = .True.
  sfields(id% Alk)    % output(aa,oo) = .True.
  sfields(id% pCO2s)  % output(aa,oo) = .True.
  sfields(id% O2)     % output(aa,oo) = .True.
  sfields(id% DIN)    % output(aa,oo) = .True.
  ! set to daily
  sfields(id% PhyChl) % output(aa,dd) = .True.
  sfields(id% DiaChl) % output(aa,dd) = .True.
  sfields(id% DIC)    % output(aa,dd) = .True.
  sfields(id% Alk)    % output(aa,dd) = .True.
  sfields(id% pCO2s)  % output(aa,dd) = .True.
  sfields(id% O2)     % output(aa,dd) = .True.
  sfields(id% DIN)    % output(aa,dd) = .True.
ENDIF


! ________________________
! ___ FINALIZE        ____


DO s=1, nfields
  ! no monthly initial fields. monthly fields are written at last day of month, but initial fields at first day.
  ! because monthly initial fields make no sense, we set "daily" (dd=True) for all initial fields.
  sfields(s)% output(ii,dd) = .True.
  sfields(s)% output(si,dd) = .True.
  
  ! do not compute full ensemble state for m-fields! this takes memory and time. simply use fesom-output instead.
  ! whatever settings made be before, we reset m-fields ens-member output to False, in the end.
  sfields(s)% output(mm,ee) = .False.
  
  ! the standard deviation is written to ensemble-mean file
  sfields(s)% output(si,ee) = .False.
  sfields(s)% output(sf,ee) = .False.
  sfields(s)% output(sa,ee) = .False.
  sfields(s)% output(sm,ee) = .False.
  
ENDDO

! __________________________________________
! ___ which files are needed?           ____

w_daymemb = .false.          
w_dayensm = .false.          
w_monmemb = .false.          
w_monensm = .false.          

DO s=1, nfields
  !                               -- output? -----------        -- ensemble members? ----------       -- daily? ----------------------
  w_daymemb = w_daymemb .or. any( sfields(s)%output(:,oo) .and. (      sfields(s)%output(:,ee)) .and. (      sfields(s)%output(:,dd)))  ! have any daily member states to write?
  w_dayensm = w_dayensm .or. any( sfields(s)%output(:,oo) .and. (.not. sfields(s)%output(:,ee)) .and. (      sfields(s)%output(:,dd)))  ! have any daily ensemble mean states to write?
  w_monmemb = w_monmemb .or. any( sfields(s)%output(:,oo) .and. (      sfields(s)%output(:,ee)) .and. (.not. sfields(s)%output(:,dd)))  ! have any monthly member states to write?
  w_monensm = w_monensm .or. any( sfields(s)%output(:,oo) .and. (.not. sfields(s)%output(:,ee)) .and. (.not. sfields(s)%output(:,dd)))  ! have any monthly ensemble mean states to write?
ENDDO

w_dayensm = .true. ! daily file to protocol forgetting factor and global standard deviation

! _____________________________________
w_mm = any(sfields(:)%output(mm,oo))    ! any m-fields (day-mean) to be written?
w_sm = any(sfields(:)%output(sm,oo))    ! any m-fields (day-mean) of standard deviation to be written?

! output message
IF (mype_world==0) THEN
  DO s=1, nfields
  
    outputmessage(aa,1) = 'Analysis'
    outputmessage(ff,1) = 'Forecast'
    outputmessage(mm,1) = 'M-Field'
    outputmessage(ii,1) = 'Initial'
    outputmessage(sf,1) = 'STD Forcast'
    outputmessage(sa,1) = 'STD Analysis'
    outputmessage(si,1) = 'STD Initial'
    outputmessage(sm,1) = 'STD M-Field'
    
    DO j=1,8
      IF (sfields(s)%output(j,ee)) then
        outputmessage(j,ee) = 'Ens-Memb'
      ELSE
        outputmessage(j,ee) = 'Ens-Mean'
      ENDIF
      IF (sfields(s)%output(j,dd)) then
        outputmessage(j,dd) = 'Daily'
      ELSE
        outputmessage(j,dd) = 'Monthly'
      ENDIF
    ENDDO ! j=1,4
    
    outputmessage(ii,dd) = '-'
    outputmessage(si,dd) = '-'
  
    if (.not. any(sfields(s)%output(:,oo))) then ! this field any output?
      write (*, '(a,4x,a,1x,a10,1x,a9)') 'FESOM-PDAF', 'Field', sfields(s)%variable, 'No Output'
    else
      DO j=1,4
        if (sfields(s)%output(j,oo)) write (*, '(a,4x,a,1x,a10,1x,a12,1x,a10,1x,a10,1x,a10)') 'FESOM-PDAF', 'Field', &
                                                                                               sfields(s)%variable, &
                                                                                               outputmessage(j, 1), &
                                                                                               outputmessage(j,ee), &
                                                                                               outputmessage(j,dd)
      ENDDO ! j=1,8
    endif ! this field any output?
  ENDDO ! s=1, nfields
ENDIF ! writepe

END SUBROUTINE init_sfields
  
  
END MODULE mod_nc_out_variables
