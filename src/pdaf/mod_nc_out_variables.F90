MODULE mod_nc_out_variables

! USES:
USE mod_assim_pdaf, &
   ONLY: id, nfields

IMPLICIT NONE

character(len=200) :: filename_phy = ''       ! Full name of output file
character(len=200) :: filename_bgc = ''       ! Full name of output file

LOGICAL :: write_ens    = .true.          ! Whether to write ensemble states
INTEGER :: write_pos_da = 1


! Field description:

integer :: id_var                         ! Index of a variable in state vector

type state_field
   integer :: ndims = 0                   ! Number of field dimensions (1 or 2)
   logical :: nz1 = .true.                ! Vertical coordinates (on levels / on layers)
   character(len=10) :: variable = ''     ! Name of field
   character(len=40) :: long_name = ''    ! Long name of field
   character(len=20) :: units = ''        ! Unit of variable
   integer :: varid(3)                    ! To write to netCDF file
   logical :: updated = .true.            ! Whether variable is updated through assimilation
   logical :: bgc = .false.               ! Whether variable is biogeochemistry (or physics)
end type state_field

type(state_field), allocatable :: sfields(:) ! Type variable holding the
                                             ! definitions of model fields
                                             
CONTAINS

SUBROUTINE init_sfields()
ALLOCATE(sfields(nfields))

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
sfields(id% temp) % bgc = .false.


! salt
sfields(id% salt) % ndims = 2
sfields(id% salt) % nz1 = .true.
sfields(id% salt) % variable = 'S'
sfields(id% salt) % long_name = 'Salinity'
sfields(id% salt) % units = 'psu'
sfields(id% salt) % updated = .true.
sfields(id% salt) % bgc = .false.


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
sfields(id% MLD1) % long_name = 'Mixed layer depth (Large et al. 1997)'
sfields(id% MLD1) % units = 'm'
sfields(id% MLD1) % updated = .false.
sfields(id% MLD1) % bgc = .false.

! chlorophyll small phytoplankton
sfields(id% PhyChl) % ndims = 2
sfields(id% PhyChl) % nz1 = .true.
sfields(id% PhyChl) % variable = 'PhyChl'
sfields(id% PhyChl) % long_name = 'Chlorophyll-a small phytoplankton'
sfields(id% PhyChl) % units = 'mg chl m-3'
sfields(id% PhyChl) % updated = .false.
sfields(id% PhyChl) % bgc = .true.

! chlorophyll diatoms
!~ sfields(id% DiaChl) % ndims = 2
!~ sfields(id% DiaChl) % nz1 = .true.
!~ sfields(id% DiaChl) % variable = 'DiaChl'
!~ sfields(id% DiaChl) % long_name = 'Chlorophyll-a diatoms'
!~ sfields(id% DiaChl) % units = 'mg chl m-3'
!~ sfields(id% DiaChl) % updated = .false.
!~ sfields(id% DiaChl) % bgc = .true.

! DIC
!~ sfields(id% DIC) % ndims = 2
!~ sfields(id% DIC) % nz1 = .true.
!~ sfields(id% DIC) % variable = 'DIC'
!~ sfields(id% DIC) % long_name = 'Dissolved inorganic carbon'
!~ sfields(id% DIC) % units = 'mmol m-3'
!~ sfields(id% DIC) % updated = .false.
!~ sfields(id% DIC) % bgc = .true.

! DOC
!~ sfields(id% DOC) % ndims = 2
!~ sfields(id% DOC) % nz1 = .true.
!~ sfields(id% DOC) % variable = 'DOC'
!~ sfields(id% DOC) % long_name = 'Dissolved organic carbon'
!~ sfields(id% DOC) % units = 'mmol m-3'
!~ sfields(id% DOC) % updated = .false.
!~ sfields(id% DOC) % bgc = .true.

! Alkalinity
!~ sfields(id% Alk) % ndims = 2
!~ sfields(id% Alk) % nz1 = .true.
!~ sfields(id% Alk) % variable = 'Alk'
!~ sfields(id% Alk) % long_name = 'Alkalinity'
!~ sfields(id% Alk) % units = 'mmol m-3'
!~ sfields(id% Alk) % updated = .false.
!~ sfields(id% Alk) % bgc = .true.

! DIN
!~ sfields(id% DIN) % ndims = 2
!~ sfields(id% DIN) % nz1 = .true.
!~ sfields(id% DIN) % variable = 'DIN'
!~ sfields(id% DIN) % long_name = 'Dissolved inorganic nitrogen'
!~ sfields(id% DIN) % units = 'mmol m-3'
!~ sfields(id% DIN) % updated = .false.
!~ sfields(id% DIN) % bgc = .true.

! DON
!~ sfields(id% DON) % ndims = 2
!~ sfields(id% DON) % nz1 = .true.
!~ sfields(id% DON) % variable = 'DON'
!~ sfields(id% DON) % long_name = 'Dissolved organic nitrogen'
!~ sfields(id% DON) % units = 'mmol m-3'
!~ sfields(id% DON) % updated = .false.
!~ sfields(id% DON) % bgc = .true.

! Oxygen
!~ sfields(id% O2) % ndims = 2
!~ sfields(id% O2) % nz1 = .true.
!~ sfields(id% O2) % variable = 'O2'
!~ sfields(id% O2) % long_name = 'Oxygen'
!~ sfields(id% O2) % units = 'mmol m-3'
!~ sfields(id% O2) % updated = .false.
!~ sfields(id% O2) % bgc = .true.

! pCO2
!~ sfields(id% pCO2s) % ndims = 1
!~ sfields(id% pCO2s) % variable = 'pCO2s'
!~ sfields(id% pCO2s) % long_name = 'Partial pressure CO2 surface ocean'
!~ sfields(id% pCO2s) % units = 'micro atm'
!~ sfields(id% pCO2s) % updated = .false.
!~ sfields(id% pCO2s) % bgc = .true.

! CO2f
!~ sfields(id% CO2f) % ndims = 1
!~ sfields(id% CO2f) % variable = 'CO2f'
!~ sfields(id% CO2f) % long_name = 'CO2 flux from atmosphere into ocean'
!~ sfields(id% CO2f) % units = 'mmolC m-2 d-1'
!~ sfields(id% CO2f) % updated = .false.
!~ sfields(id% CO2f) % bgc = .true.


END SUBROUTINE init_sfields
  
  
END MODULE mod_nc_out_variables
