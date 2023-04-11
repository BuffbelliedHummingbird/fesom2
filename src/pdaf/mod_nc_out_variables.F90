MODULE mod_nc_out_variables

! USES:
USE mod_assim_pdaf, &
   ONLY: id, nfields

IMPLICIT NONE

character(len=200) :: filename_phy = ''       ! Full name of output file
character(len=200) :: filename_bio = ''       ! Full name of output file

! TO-DO: INITIALIZE THIS FROM NAMELIST!
LOGICAL :: write_ens_snapshot = .true. ! Whether to write ensemble states


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

! u
sfields(id% u) % ndims = 2
sfields(id% u) % nz1 = .true.
sfields(id% u) % variable = 'u'
sfields(id% u) % long_name = 'Zonal velocity (interpolated on nodes)'
sfields(id% u) % units = 'm/s'
sfields(id% u) % updated = .true.


! v
sfields(id% v) % ndims = 2
sfields(id% v) % nz1 = .true.
sfields(id% v) % variable = 'v'
sfields(id% v) % long_name = 'Meridional velocity (interpolated on nodes)'
sfields(id% v) % units = 'm/s'
sfields(id% v) % updated = .true.


! w
sfields(id% w) % ndims = 2
sfields(id% w) % nz1 = .false.
sfields(id% w) % variable = 'w'
sfields(id% w) % long_name = 'Vertical velocity'
sfields(id% w) % units = 'm/s'
sfields(id% w) % updated = .false.


! temp
sfields(id% temp) % ndims = 2
sfields(id% temp) % nz1 = .true.
sfields(id% temp) % variable = 'T'
sfields(id% temp) % long_name = 'Temperature'
sfields(id% temp) % units = 'degC'

! salt
sfields(id% salt) % ndims = 2
sfields(id% salt) % nz1 = .true.
sfields(id% salt) % variable = 'S'
sfields(id% salt) % long_name = 'Salinity'
sfields(id% salt) % units = 'psu'
sfields(id% salt) % updated = .true.


! ice
sfields(id% a_ice) % ndims = 1
sfields(id% a_ice) % variable = 'ice'
sfields(id% a_ice) % long_name = 'Sea-ice concentration'
sfields(id% a_ice) % units = '1'
sfields(id% a_ice) % updated = .false.


END SUBROUTINE init_sfields
  
  
END MODULE mod_nc_out_variables
