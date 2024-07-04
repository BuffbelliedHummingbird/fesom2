MODULE mod_nc_out_variables

! USES:
USE mod_assim_pdaf, &
   ONLY: id, nfields, assimilateBGC, assimilatePHY, &
         phymin, phymax, bgcmin, bgcmax
USE mod_parallel_pdaf, &
   ONLY: writepe, mype_world

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
   character(len=50) :: long_name = ''    ! Long name of field
   character(len=20) :: units = ''        ! Unit of variable
   integer :: varid(4)                    ! To write to netCDF file
   logical :: updated = .true.            ! Whether variable is updated through assimilation
   logical :: bgc = .false.               ! Whether variable is biogeochemistry (or physics)
end type state_field

type(state_field), allocatable :: sfields(:) ! Type variable holding the
                                             ! definitions of model fields
                                             
INTEGER :: nfields_3D                     ! Number of 3D fields in state vector
INTEGER, ALLOCATABLE :: ids_3D(:)         ! List of 3D-field IDs
                                             
CONTAINS

SUBROUTINE init_sfields()

! Local variables:
CHARACTER(len=100) :: nmlfile ='namelist.fesom.pdaf'    ! name of namelist file
LOGICAL            :: upd_ssh , &
                      upd_u   , &
                      upd_v   , &
                      upd_w   , &
                      upd_temp, &
                      upd_salt, &
                      upd_ice , &
                      upd_MLD1, &
                      upd_PhyChl   , &
                      upd_DiaChl   , &
                      upd_DIC      , &
                      upd_DOC      , &
                      upd_Alk      , &
                      upd_DIN      , &
                      upd_DON      , &
                      upd_O2       , &
                      upd_pCO2s    , &
                      upd_CO2f     , &
                      upd_DiaN     , &
                      upd_DiaC     , &
                      upd_PAR      , &
                      upd_NPPn     , &
                      upd_NPPd     , &
                      upd_DetC     , &
                      upd_PhyCalc  , &
                      upd_export   , &
                      upd_alphaCO2 , &
                      upd_PistonVel, &
                      upd_Zo1C     , &
                      upd_Zo1N     , &
                      upd_Zo2C     , &
                      upd_Zo2N
INTEGER            :: b,p ! counters

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
sfields(id% DiaChl) % ndims = 2
sfields(id% DiaChl) % nz1 = .true.
sfields(id% DiaChl) % variable = 'DiaChl'
sfields(id% DiaChl) % long_name = 'Chlorophyll-a diatoms'
sfields(id% DiaChl) % units = 'mg chl m-3'
sfields(id% DiaChl) % updated = .false.
sfields(id% DiaChl) % bgc = .true.

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

! DetC
sfields(id% DetC) % ndims = 2
sfields(id% DetC) % variable = 'DetC'
sfields(id% DetC) % long_name = 'carbon in detritus'
sfields(id% DetC) % units = 'mmol C m-3'
sfields(id% DetC) % updated = .false.
sfields(id% DetC) % bgc = .true.

!~ ! TOC
!~ sfields(id% TOC) % ndims = 2
!~ sfields(id% TOC) % variable = 'TOC'
!~ sfields(id% TOC) % long_name = 'Total Organic Carbon (PhyC+DiaC+DetC+DOC+HetC)'
!~ sfields(id% TOC) % units = 'mmol C m-3'
!~ sfields(id% TOC) % updated = .false.
!~ sfields(id% TOC) % bgc = .true.

! PhyCalc
sfields(id% PhyCalc) % ndims = 2
sfields(id% PhyCalc) % variable = 'PhyCalc'
sfields(id% PhyCalc) % long_name = 'calcium carbonate small phytoplankton'
sfields(id% PhyCalc) % units = 'mmol m-3'
sfields(id% PhyCalc) % updated = .false.
sfields(id% PhyCalc) % bgc = .true.

! Export production
sfields(id% export) % ndims = 1
sfields(id% export) % variable = 'export'
sfields(id% export) % long_name = 'export through particle sinking at 190m'
sfields(id% export) % units = 'mmol m-2 day-1'
sfields(id% export) % updated = .false.
sfields(id% export) % bgc = .true.

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

! ************************************************
! ***   Read updated variables from namelist   ***
! ************************************************

! *** Read namelist file ***
  IF (mype_world==0) WRITE(*,*) 'Read namelist file for updated variables: ',nmlfile
  
  NAMELIST /updated/ &
     upd_ssh , &
     upd_ssh , &
     upd_u   , &
     upd_v   , &
     upd_w   , &
     upd_temp, &
     upd_salt, &
     upd_ice , &
     upd_MLD1, &
     upd_PhyChl   , &
     upd_DiaChl   , &
     upd_DIC      , &
     upd_DOC      , &
     upd_Alk      , &
     upd_DIN      , &
     upd_DON      , &
     upd_O2       , &
     upd_pCO2s    , &
     upd_CO2f     , &
     upd_DiaN     , &
     upd_DiaC     , &
     upd_PAR      , &
     upd_NPPn     , &
     upd_NPPd     , &
     upd_DetC     , &
     upd_PhyCalc  , &
     upd_export   , &
     upd_alphaCO2 , &
     upd_PistonVel, &
     upd_Zo1C     , &
     upd_Zo1N     , &
     upd_Zo2C     , &
     upd_Zo2N
     

  OPEN  (20,file=nmlfile)
  READ  (20,NML=updated)
  CLOSE (20)
  
  sfields(id% ssh      ) % updated = upd_ssh
  sfields(id% ssh      ) % updated = upd_ssh
  sfields(id% u        ) % updated = upd_u
  sfields(id% v        ) % updated = upd_v
  sfields(id% w        ) % updated = upd_w
  sfields(id% temp     ) % updated = upd_temp
  sfields(id% salt     ) % updated = upd_salt
  sfields(id% a_ice    ) % updated = upd_ice
  sfields(id% MLD1     ) % updated = upd_MLD1
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
  sfields(id% DiaN     ) % updated = upd_DiaN
  sfields(id% DiaC     ) % updated = upd_DiaC
  sfields(id% PAR      ) % updated = upd_PAR
  sfields(id% NPPn     ) % updated = upd_NPPn
  sfields(id% NPPd     ) % updated = upd_NPPd
  sfields(id% DetC     ) % updated = upd_DetC
  sfields(id% PhyCalc  ) % updated = upd_PhyCalc
  sfields(id% export   ) % updated = upd_export
  sfields(id% alphaCO2 ) % updated = upd_alphaCO2
  sfields(id% PistonVel) % updated = upd_PistonVel
  sfields(id% Zo1C     ) % updated = upd_Zo1C
  sfields(id% Zo1N     ) % updated = upd_Zo1N
  sfields(id% Zo2C     ) % updated = upd_Zo2C
  sfields(id% Zo2N     ) % updated = upd_Zo2N
  
  ! If BGC or physics are not assimilated,
  ! do not update the respective part of state vector.  
  IF (.not. assimilatePHY) THEN
     do p=phymin,phymax
       sfields(p) % updated = .false.
     enddo
  ELSEIF (.not. assimilateBGC) THEN
     do b=bgcmin,bgcmax
       sfields(b) % updated = .false.
     enddo
  ENDIF
   
  ! count number of 3D fields
  nfields_3D = 0
  DO b=1,nfields
   IF (sfields(b) % ndims == 2) nfields_3D = nfields_3D + 1
  ENDDO
  
  ! indeces of 3D fields
  ALLOCATE(ids_3D(nfields_3D))
  p = 1
  DO b=1,nfields
    IF (sfields(b) % ndims == 2) THEN
      ids_3D(p) = b
      p = p+1
    ENDIF
  ENDDO
  IF (mype_world==0) THEN
     DO p=1,nfields_3D
     WRITE (*,'(a, 10x,3a,1x,7a)') &
          'FESOM-PDAF', '3D fields in state vector: ', sfields % variable(ids_3D(p))
  ENDIF

END SUBROUTINE init_sfields
  
  
END MODULE mod_nc_out_variables
