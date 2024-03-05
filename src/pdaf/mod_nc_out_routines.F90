MODULE mod_nc_out_routines

! contains:
!   - check()
!   - netCDF_init()
!        -> calls netCDF_init_mm()
!   - netCDF_init_rmsens()
!   - netCDF_out()


! USES:
USE mod_nc_out_variables
USE mod_assim_pdaf, &
   ONLY: nfields, id, mesh_fesom, DAoutput_path, dim_ens, dim_state, &
         dim_state_p, offset, phymin, phymax, bgcmin, bgcmax
USE mod_parallel_pdaf, &
   ONLY: abort_parallel, writepe
USE mod_obs_f_pdaf, &
   ONLY: pi
USE g_config, &
   ONLY: runid
USE g_clock, &
   ONLY: cyearnew
USE g_parsup, &
   ONLY: myDim_nod2D
USE g_comm_auto, &
   ONLY: gather_nod
USE netcdf

IMPLICIT NONE

! Private variables:
INTEGER :: i                          ! counter (state variables)
INTEGER :: j                          ! counter (ini/forc/ana)
INTEGER :: member                     ! counter (ensemble member)
INTEGER :: n                          ! counter (nodes)
INTEGER :: l                          ! counter (levels/layers)

CHARACTER(len=2) :: memberstr         ! ensemble member
INTEGER :: fileid                     ! ID of netCDF file
CHARACTER, dimension(4) :: IFA = ['i', 'f', 'a', 'm']
                                      ! "i" (initial), 
                                      ! "f" (forecast),
                                      ! "a" (analysis),
                                      ! "m" (daily means)
CHARACTER(len=8), dimension(4) :: IFA_long = ['initial','forecast','analysis','daymean']
INTEGER :: dimID_nod2, dimID_iter, &
           dimID_nz, dimID_nz1
INTEGER :: dimIDs(3)
INTEGER :: varid_nod2, varid_iter, &
           varid_nz, varid_nz1, &
           varid_lon, varid_lat, &
           varid_forget
INTEGER :: ndims
INTEGER :: nlay

REAL, allocatable :: lon(:)
REAL, allocatable :: lat(:)

CONTAINS


! ********************************
! ***                          ***
! ***   netCDF check           ***
! ***                          ***
! ********************************
! Checks for errors during netCDF operations.

SUBROUTINE check(status)

! *** Arguments ***
integer, intent ( in) :: status   ! Reading status

    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       call abort_parallel()
    end if

END SUBROUTINE check

! ********************************
! ***                          ***
! ***   netCDF_init            ***
! ***                          ***
! ********************************
! Calls netCDF_init_mm on main PE depending on whether "m"ean or ensemble "m"ember
! file(s) are to be initialized.

SUBROUTINE netCDF_init(mm)

! *** Arguments ***
CHARACTER(len=4) :: mm ! Mean state or member?

! gather GEO coordinates (from all PEs)
allocate(lon(mesh_fesom% nod2D),lat(mesh_fesom% nod2D))
call gather_nod(mesh_fesom%geo_coord_nod2D(1, 1:myDim_nod2D), lon)
call gather_nod(mesh_fesom%geo_coord_nod2D(2, 1:myDim_nod2D), lat)

! initialize file (on main PE)
IF (writepe) THEN
  IF (mm=='mean') THEN
  
     filename_phy = TRIM(DAoutput_path)//'fesom-recom-pdaf.phy.'//cyearnew//'.nc'
     filename_bgc = TRIM(DAoutput_path)//'fesom-recom-pdaf.bgc.'//cyearnew//'.nc'
     
     call netCDF_init_mm    (filename_phy, phymin, phymax)
     call netCDF_init_rmsens(filename_phy, phymin, phymax)
     
     call netCDF_init_mm    (filename_bgc, bgcmin, bgcmax)
     call netCDF_init_rmsens(filename_bgc, bgcmin, bgcmax)
  
  ELSEIF (mm=='memb') THEN
     DO member=1, dim_ens
     
        WRITE(memberstr,'(i2.2)') member
        filename_phy = TRIM(DAoutput_path)//'fesom-recom-pdaf.phy.'//cyearnew//'.'//memberstr//'.nc'
        filename_bgc = TRIM(DAoutput_path)//'fesom-recom-pdaf.bgc.'//cyearnew//'.'//memberstr//'.nc'

        call netCDF_init_mm(filename_phy, phymin, phymax)
        call netCDF_init_mm(filename_bgc, bgcmin, bgcmax)
     
     ENDDO
  ELSE
     WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initializing mean or member netCDF file?'
  ENDIF ! mean/memb

ENDIF ! writepe
deallocate(lat,lon)


END SUBROUTINE netCDF_init

! ********************************
! ***                          ***
! ***   netCDF_init_mm         ***
! ***                          ***
! ********************************
! Initializes a netCDF file.

SUBROUTINE netCDF_init_mm(filename, istart, iend)

character(len=200), intent(in) :: filename ! Full name of output file
integer, intent(in) :: istart              ! First field in state vector
integer, intent(in) :: iend                ! Last field in state vector

WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize assimilation NetCDF file: '//TRIM(filename)

! open file
call check(NF90_CREATE(trim(filename),NF90_NETCDF4,fileid))

! define dimensions
call check( NF90_DEF_DIM(fileid, 'nod2', mesh_fesom%nod2d, dimID_nod2))
call check( NF90_DEF_DIM(fileid, 'nz',   mesh_fesom%nl,    dimID_nz))
call check( NF90_DEF_DIM(fileid, 'nz1',  mesh_fesom%nl-1,  dimID_nz1))
call check( NF90_DEF_DIM(fileid, 'time', NF90_UNLIMITED,   dimID_iter))

! dimension variables
call check( nf90_def_var(fileid, 'nod2', NF90_INT,   dimID_nod2, varid_nod2))
call check( nf90_def_var(fileid, 'nz',   NF90_FLOAT, dimID_nz,   varid_nz))
call check( nf90_def_var(fileid, 'nz1',  NF90_FLOAT, dimID_nz1,  varid_nz1))
call check( nf90_def_var(fileid, 'time', NF90_INT,   dimID_iter, varid_iter))

call check( nf90_def_var(fileid, 'lon', NF90_FLOAT,   dimID_nod2, varid_lon))
call check( nf90_def_var(fileid, 'lat', NF90_FLOAT,   dimID_nod2, varid_lat))

! dimension description
call check( nf90_put_att(fileid, varid_nod2, 'long_name', 'surface nodes'))
call check( nf90_put_att(fileid, varid_nz,   'long_name', 'vertical levels (upper layer bounds)'))
call check( nf90_put_att(fileid, varid_nz,   'units', 'm'))
call check( nf90_put_att(fileid, varid_nz1,  'long_name', 'vertical layers (mid-layer depths)'))
call check( nf90_put_att(fileid, varid_nz1,  'units', 'm'))
call check( nf90_put_att(fileid, varid_iter, 'long_name', 'iteration'))

call check( nf90_put_att(fileid, varid_lon,  'long_name', 'longitude'))
call check( nf90_put_att(fileid, varid_lat,  'long_name', 'latitude'))
call check( nf90_put_att(fileid, varid_lon,  'units', 'degE [-180;180]'))
call check( nf90_put_att(fileid, varid_lat,  'units', 'degN [-90;90]'))


! fill dimension variables
call check (nf90_put_var(fileid, varid_nz,   mesh_fesom% zbar))
call check (nf90_put_var(fileid, varid_nz1,  mesh_fesom% Z))
call check (nf90_put_var(fileid, varid_nod2, [(n,n=1,mesh_fesom% nod2D)]))

call check (nf90_put_var(fileid, varid_lon,  REAL(180./pi * lon, 4)))
call check (nf90_put_var(fileid, varid_lat,  REAL(180./pi * lat, 4)))


! define state variables
DO j = 1, 3 ! ini / forc / ana
  DO i = istart, iend ! state fields
    
    ! don't write analysis for not-updated variables
    IF ((j==3) .and. (.not. (sfields(i)% updated))) CYCLE
    
    ! set dimensions according to field before defining variables
    dimIDs(1) = dimID_nod2
    IF (sfields(i)% ndims == 1) THEN ! surface fields
      dimIDs(2) = dimID_iter
    ELSEIF (sfields(i)% ndims == 2) THEN ! 3D-fields
      dimIDs(2) = dimID_nz1
      IF (.not.(sfields(i)%nz1)) dimIDs(2) = dimID_nz ! variables on levels (vertical velocity)
      dimIDs(3) = dimID_iter
    ENDIF
    
    ! number of dimensions
    IF (IFA(j)=='i') THEN
      ndims = sfields(i)% ndims     ! no iteration
    ELSE
      ndims = sfields(i)% ndims +1  ! plus iteration
    ENDIF
    
    ! define variable
    call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
    ! variable description
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))

  ENDDO ! state fields
ENDDO ! ini / forc / ana

call check(NF90_ENDDEF(fileid))
call check(NF90_CLOSE(fileid))

END SUBROUTINE netCDF_init_mm

! ********************************
! ***                          ***
! ***   netCDF_out             ***
! ***                          ***
! ********************************

SUBROUTINE netCDF_out(writetype, writepos, iteration, state_p, ens_p, rms, forget)

! ARGUMENTS:
CHARACTER(len=1), intent(in) :: writetype          ! Write (i) initial, (a) assimilated, (f) forecast, (m) daily-average fields
INTEGER, INTENT(in) :: writepos                    ! Write position
INTEGER, INTENT(in) :: iteration                   ! Current model time step
REAL, target, INTENT(in)    :: state_p(dim_state_p)        ! Ensemble mean state vector
REAL, target, INTENT(in)    :: ens_p(dim_state_p, dim_ens) ! Ensemble states
REAL, INTENT(in) :: rms(nfields)                   ! RMS ensemble spread
REAL, INTENT(in) :: forget                         ! Forgetting factor

! Local variables:
REAL, pointer :: state_ptr(:)                      ! Points to state vector (mean or ensemble member state)
REAL, allocatable :: myData2(:)                    ! Temporary array for pe-local surface fields
REAL, allocatable :: myData3(:,:)                  ! Temporary array for pe-local 3D-fields
REAL, allocatable :: data2_g(:)                    ! Temporary array for global surface fields
REAL, allocatable :: data3_g(:,:)                  ! Temporary array for global 3D-fields

! character(len=200) :: filename                     ! Full name of output file
integer            :: fileid_bgc, fileid_phy        ! nc-file ID for biogeochemistry and physics output file

allocate(state_ptr(dim_state_p))

! Print screen information:
IF (writepe) THEN
  IF (writetype == 'i') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write initial ocean state to NetCDF at step ', &
          iteration, ' position ', writepos
  ELSE IF (writetype== 'f') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean forecast to NetCDF at step ', &
          iteration, ' position ', writepos
  ELSE IF (writetype== 'a') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean analysis to NetCDF at step ', &
          iteration, ' position ', writepos
  ELSE IF (writetype== 'm') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ocean daily mean to NetCDF at step ', &
          iteration, ' position ', writepos
  END IF
END IF ! writepe

DO member=0, dim_ens
  
  ! Set file name and ensemble mean or member data:
  IF (member==0) THEN ! mean state
    filename_phy = TRIM(DAoutput_path)//'fesom-recom-pdaf.phy.'//cyearnew//'.nc'
    filename_bgc = TRIM(DAoutput_path)//'fesom-recom-pdaf.bgc.'//cyearnew//'.nc'
    state_ptr => state_p
  ELSE ! ensemble member state
    IF (.not. (write_ens)) EXIT
    IF (writetype== 'm') EXIT
    WRITE(memberstr,'(i2.2)') member
    filename_phy = TRIM(DAoutput_path)//'fesom-recom-pdaf.phy.'//cyearnew//'.'//memberstr//'.nc'
    filename_bgc = TRIM(DAoutput_path)//'fesom-recom-pdaf.bgc.'//cyearnew//'.'//memberstr//'.nc'
    state_ptr => ens_p(:,member)
  ENDIF ! mean/member
  
  ! Open netCDF file:
  IF (writepe) THEN
    call check( nf90_open(TRIM(filename_phy), nf90_write, fileid_phy))
    call check( nf90_open(TRIM(filename_bgc), nf90_write, fileid_bgc))
  
    ! Non field-specific writing:
    IF (writetype=='a') THEN
      call check( nf90_inq_varid(fileid_phy, "time", varid_iter))
      call check( nf90_inq_varid(fileid_bgc, "time", varid_iter))
      
      call check( nf90_put_var  (fileid_phy, varid_iter, iteration, &
                                 start=(/ writepos /)))
      call check( nf90_put_var  (fileid_bgc, varid_iter, iteration, &
                                 start=(/ writepos /)))
                                 
      IF (member==0) THEN
        call check( nf90_inq_varid(fileid_phy, "forget", varid_forget))
        call check( nf90_put_var  (fileid_phy, varid_forget, REAL(forget,4), &
                                   start=(/ writepos /)))
      ENDIF ! member
    ENDIF ! writetype
  ENDIF ! writepe

  ! Field-specific:
  DO i=1, nfields
    
    ! don't write analysis of not-updated variables
    IF ((writetype=='a') .and. (.not. (sfields(i)% updated))) CYCLE
    
    ! physics or biogeochemistry?
    IF (sfields(i)% bgc) THEN
      fileid = fileid_bgc
    ELSE
      fileid = fileid_phy
    ENDIF
    
    ! --------------
    ! surface fields
    ! --------------
    IF (sfields(i)% ndims == 1) THEN
      ! select pe-local field from state vector
      allocate(myData2(myDim_nod2d))
      DO n = 1, myDim_nod2D
        myData2(n) = state_ptr(n + offset(i))
      END DO
      ! gather global field
      allocate(data2_g(mesh_fesom% nod2D))
      CALL gather_nod(myData2, data2_g)
      
      IF (writepe) THEN
        ! Inquire variable ID
        call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(1)))
        ! Write variable to netCDF:
        IF (writetype=='i') THEN ! initial field
          call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(data2_g,4)))
        ELSE                     ! forecast, analysis and mean fields
          call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(data2_g,4), &
                                   start=(/ 1, writepos /), &
                                   count=(/ mesh_fesom% nod2D, 1 /) ))
                                   ! (dims: 1-nod2, 2-iter)
        ENDIF ! writetype (i,f,a,m)
      ENDIF ! writepe
      deallocate(myData2, data2_g)
      
    ! ---------
    ! 3D fields
    ! ---------
    ELSEIF ((sfields(i)% ndims == 2)) THEN
      ! Levels or layers?
      IF (sfields(i)% nz1) THEN
        nlay = mesh_fesom% nl-1
      ELSE
        nlay = mesh_fesom% nl
      ENDIF
      ! select pe-local field from state vector
      allocate(myData3(nlay, myDim_nod2d))
      DO n = 1, myDim_nod2D
      DO l = 1, nlay
        myData3(l, n) = state_ptr((n-1) * nlay + l + offset(i))
      END DO
      END DO
      ! gather global field
      allocate(data3_g(nlay,mesh_fesom% nod2D))
      CALL gather_nod(myData3, data3_g)
      
      IF (writepe) THEN
        ! Inquire variable ID
        call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(1)))
        ! Write variable to netCDF:
        IF (writetype=='i') THEN ! initial field
          call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(TRANSPOSE(data3_g),4)))
        ELSE                     ! forecast, analysis and mean fields
          call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(TRANSPOSE(data3_g),4), &
                                   start=(/ 1, 1, writepos /), &
                                   count=(/ mesh_fesom% nod2D, nlay, 1 /) ))
                                   ! dims: 1-nod2, 2-nz / nz1, 3-iter
        ENDIF ! writetype
      ENDIF ! writepe
      deallocate(myData3, data3_g)
      
    ENDIF ! surface / 3D-fields
    
    ! ---------------------
    ! RMS (Ensemble spread)
    ! ---------------------
    IF (member/=0) CYCLE
    IF (writetype=='m') CYCLE
    
    IF (writepe) THEN
      ! Inquire variable ID
      call check( nf90_inq_varid(fileid, 'rms_'//TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(2)))
      ! Write RMS to netCDF
      call check( nf90_put_var(fileid, sfields(i)% varid(2), REAL(rms(i),4), &
                               start=(/ writepos /)))
    ENDIF ! writepe
                             
  ENDDO ! nfields
  
  ! Close file:
  IF (writepe) THEN
    call check (nf90_close(fileid_phy))
    call check (nf90_close(fileid_bgc))
  ENDIF ! writepe
  
ENDDO ! dim_ens
END SUBROUTINE netCDF_out

! ********************************
! ***                          ***
! ***   netCDF_init_rmsens     ***
! ***                          ***
! ********************************

! Inits netCDF output for ensemble spread (root-mean-square deviation)

SUBROUTINE netCDF_init_rmsens(filename, istart, iend)

character(len=200), intent(in) :: filename ! Full name of output file
integer, intent(in) :: istart              ! First field in state vector
integer, intent(in) :: iend                ! Last field in state vector

! open file
call check( nf90_open(TRIM(filename), nf90_write, fileid))

! inquire dimension ID
call check( nf90_inq_dimid( fileid, "time", dimID_iter))
call check( nf90_inq_dimid( fileid, "nod2", dimID_nod2))
call check( nf90_inq_dimid( fileid, "nz",   dimID_nz))
call check( nf90_inq_dimid( fileid, "nz1",  dimID_nz1))

! ********************
! define rms variables
! ********************
DO j = 1, 3 ! ini / forc / ana
  DO i = istart, iend ! state fields
    
    ! don't write analysis of not-updated variables
    IF ((j==3) .and. (.not. (sfields(i)% updated))) CYCLE
    ! set dimension
    IF (j==1) THEN
      ! initial RMS
      ! define dimensionless (1D) variable
    call check( NF90_DEF_VAR(fileid, 'rms_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, sfields(i)% varid(j)))
    ELSE
      ! forecast / analysis RMS
      ! define timeseries variable
    call check( NF90_DEF_VAR(fileid, 'rms_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimID_iter, sfields(i)% varid(j)))
    ENDIF
    
    ! variable description
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', 'Ensemble spread (RMS) for '//trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    
  ENDDO ! state fields
ENDDO ! ini / ana / forc

! ***************************
! define daily mean variables
! ***************************

j = 4 ! daily mean data
DO i = istart, iend ! state fields

    ! set dimensions according to field before defining variables
    dimIDs(1) = dimID_nod2
    IF (sfields(i)% ndims == 1) THEN ! surface fields
      dimIDs(2) = dimID_iter
    ELSEIF (sfields(i)% ndims == 2) THEN ! 3D-fields
      dimIDs(2) = dimID_nz1
      IF (.not.(sfields(i)%nz1)) dimIDs(2) = dimID_nz ! variables on levels (vertical velocity)
      dimIDs(3) = dimID_iter
    ENDIF
    
    ! number of dimensions
    ndims = sfields(i)% ndims +1  ! plus iteration
    
    ! define variable
    call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
    ! variable description
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))

  ENDDO ! state fields
  
! ***************************
! define forget
! ***************************

! define variable
call check( NF90_DEF_VAR(fileid,'forget', NF90_FLOAT, dimID_iter, varID_forget))
! variable description
call check( nf90_put_att(fileid, varID_forget, 'long_name', 'variable forgetting factor, snapshots'))


! close file
call check (nf90_close(fileid))
END SUBROUTINE netCDF_init_rmsens

END MODULE mod_nc_out_routines
