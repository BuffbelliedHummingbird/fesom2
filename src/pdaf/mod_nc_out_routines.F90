MODULE mod_nc_out_routines

! contains:
!   - check
!   - netCDF_init
!      + netCDF_deffile
!      + netCDF_defvar
!   - netCDF_out
!      + netCDF_openfile

! USES:
USE mod_nc_out_variables
USE mod_assim_pdaf, &
   ONLY: nfields, id, mesh_fesom, nlmax, DAoutput_path, dim_ens, dim_state, &
         dim_state_p, offset, phymin, phymax, bgcmin, bgcmax, &
         topography3D_g
USE mod_parallel_pdaf, &
   ONLY: abort_parallel, writepe
USE mod_obs_f_pdaf, &
   ONLY: pi
USE g_config, &
   ONLY: runid
USE g_clock, &
   ONLY: cyearnew, timeold, dayold, yearold, yearnew, num_day_in_month, fleapyear
USE g_parsup, &
   ONLY: myDim_nod2D
USE g_comm_auto, &
   ONLY: gather_nod
USE netcdf

IMPLICIT NONE
save

! Private variables:
INTEGER :: i                          ! counter (state variables)
INTEGER :: j                          ! counter (ini/forc/ana/mean)
INTEGER :: member                     ! counter (ensemble member)
INTEGER :: n                          ! counter (nodes)
INTEGER :: l                          ! counter (levels/layers)

CHARACTER(len=2)   :: memberstr       ! ensemble member
CHARACTER(len=200) :: filename        ! name of output file
INTEGER :: fileid                     ! IDs of netCDF files
INTEGER :: fidphyday, fidbgcday, &    ! IDs of netCDF files
           fidphymon, fidbgcmon
           
!~ ff=1, aa=2, mm=3, ii=4 ! forecast (ff), analysis (aa), mean (mm) and initial (ii)
!~ sf=5, sa=6, si=7, sm=8 ! ensemble standard deviation of forecast (sf), analysis (sa) and initial (si)
!~ oo=1, ee=2, dd=3       ! any output (oo), ensemble members (ee) and daily values (dd)

CHARACTER(len=2), dimension(8) :: IFA       ! Type character ('ff','aa','mm','ii','sf','sa','si')
CHARACTER(len=8), dimension(8) :: IFA_long  ! Type character
INTEGER :: dimID_nod2, dimID_time, &
           dimID_nz, dimID_iter
INTEGER :: dimIDs(3)
INTEGER :: varid_nod2, varid_nz, &    ! netCDF IDs for dimensions
           varid_time, varid_iter, &              
           varid_lon, varid_lat, &
           varid_forget               ! netCDF ID for forgetting factor
INTEGER :: ndims

REAL, allocatable :: lon(:)
REAL, allocatable :: lat(:)

REAL, parameter    :: fill_value = -999.0
INTEGER, parameter :: int0=0

LOGICAL :: debug = .false.

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
! This routine initializes netCDF output for monthly and daily physics and BGC fields.
! It is called from init_PDAF on all PEs.
! It is optionally called for the ensemble mean and/or for each ensemble member.

SUBROUTINE netCDF_init()

USE mod_assim_pdaf, only: dim_ens

INTEGER :: memb                   ! Zero:    Mean state
                                  ! Number:  Ensemble member
                                                  
! initialize:

!~ ff=1, aa=2, mm=3, ii=4 ! forecast (ff), analysis (aa), mean (mm) and initial (ii)
!~ sf=5, sa=6, si=7, sm=8 ! ensemble standard deviation of forecast (sf), analysis (sa) and initial (si)
!~ oo=1, ee=2, dd=3       ! any output (oo), ensemble members (ee) and daily values (dd)

IFA(ff)='ff' ! "f" (forecast)
IFA(aa)='aa' ! "a" (analysis)
IFA(mm)='mm' ! "m" (daily means)
IFA(ii)='ii' ! "i" (initial)
IFA(sf)='sf'
IFA(sa)='sa'
IFA(si)='si'
IFA(sm)='sm'

IFA_long(ff)='forecast'
IFA_long(aa)='analysis'
IFA_long(mm)='timemean'
IFA_long(ii)='initial'
IFA_long(sf)='STD forc'
IFA_long(sa)='STD anls'
IFA_long(si)='STD init'
IFA_long(sm)='STD mean'

! gather GEO coordinates
allocate(lon(mesh_fesom% nod2D),lat(mesh_fesom% nod2D))
call gather_nod(mesh_fesom%geo_coord_nod2D(1, 1:myDim_nod2D), lon)
call gather_nod(mesh_fesom%geo_coord_nod2D(2, 1:myDim_nod2D), lat)

! initialize file (on main PE)
IF (writepe) THEN
     
     DO memb=0,dim_ens
        ! _______________________________________
        ! (1) netCDF_deffile: create and open file, define dimensions and coordinates
        IF ((memb==0) .and. w_dayensm) THEN
           CALL netCDF_deffile('phy','day',int0,fidphyday)
           if (debug) write (*,'(a, 10x,a,1x,i0)') 'FESOM-PDAF', 'fidphyday is', fidphyday
           CALL netCDF_deffile('bgc','day',int0,fidbgcday)
           if (debug) write (*,'(a, 10x,a,1x,i0)') 'FESOM-PDAF', 'fidbgcday is', fidbgcday
        ENDIF
        
        IF ((memb >0) .and. w_daymemb) THEN
           CALL netCDF_deffile('phy','day',memb,fidphyday)
           CALL netCDF_deffile('bgc','day',memb,fidbgcday)
        ENDIF
        
        IF ((memb==0) .and. w_monensm) THEN
           CALL netCDF_deffile('phy','mon',int0,fidphymon)
           CALL netCDF_deffile('bgc','mon',int0,fidbgcmon)
        ENDIF
        
        IF ((memb >0) .and. w_monmemb) THEN
           CALL netCDF_deffile('phy','mon',memb,fidphymon)
           CALL netCDF_deffile('bgc','mon',memb,fidbgcmon)
        ENDIF
        
        ! _______________________________________
        ! (2) netCDF_defvar:  define variables and close file
        IF ((memb==0) .and. (w_dayensm .or. w_monensm)) CALL netCDF_defvar(int0)
        IF ((memb >0) .and. (w_daymemb .or. w_monmemb)) CALL netCDF_defvar(memb)

     ENDDO ! memb=0,dim_ens
     
ENDIF ! writepe
deallocate(lon,lat)
END SUBROUTINE netCDF_init

! ********************************
! ***   netCDF_deffile         ***
! ********************************
! This routine creates and opens the netCDF output file into put-var mode
! It is called from netCDF_init on writepe

SUBROUTINE netCDF_deffile(typ,freq,memb,fid)
! *** Arguments ***
CHARACTER(len=3), intent (in)  :: typ   ! phy or bgc
CHARACTER(len=3), intent (in)  :: freq  ! day or mon
INTEGER,          intent (in)  :: memb  ! Zero:    Mean state
                                        ! Number:  Ensemble member
INTEGER,          intent (out) :: fid   ! file ID
! *** Local variables ***
character(40) :: att_text
              
IF (memb==0) THEN
  ! mean state
  filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.'//typ//'.'//cyearnew//'.'//freq//'.nc'
ELSE
  ! ensemble member
  WRITE(memberstr,'(i2.2)') memb
  filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.'//typ//'.'//cyearnew//'.'//freq//'.'//memberstr//'.nc'
ENDIF

! open file
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'Init output to file '//trim(filename)

call check(NF90_CREATE(trim(filename),NF90_NETCDF4,fid))
if (debug) write (*,'(a, 10x,a,1x,i0)') 'FESOM-PDAF', 'fid is', fid

! define dimensions
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define dimensions'

call check( NF90_DEF_DIM(fid, 'nod2', mesh_fesom%nod2d, dimID_nod2))
call check( NF90_DEF_DIM(fid, 'nz1',  nlmax,            dimID_nz))
call check( NF90_DEF_DIM(fid, 'time', NF90_UNLIMITED,   dimID_time))

! dimension variables
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define dimension variables'

call check( nf90_def_var(fid, 'nod2', NF90_INT,   dimID_nod2, varid_nod2))
call check( nf90_def_var(fid, 'nz1',  NF90_FLOAT, dimID_nz,   varid_nz))
call check( nf90_def_var(fid, 'time', NF90_INT,   dimID_time, varid_time))

call check( nf90_def_var(fid, 'lon', NF90_FLOAT,   dimID_nod2, varid_lon))
call check( nf90_def_var(fid, 'lat', NF90_FLOAT,   dimID_nod2, varid_lat))

! dimension description
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define dimension description'

call check( nf90_put_att(fid, varid_nod2, 'long_name', 'surface nodes'))
call check( nf90_put_att(fid, varid_nz,   'long_name', 'vertical layers (mid-layer depths)'))
call check( nf90_put_att(fid, varid_nz,   'units', 'm'))
call check( nf90_put_att(fid, varid_time, 'long_name', 'time'))

call check( nf90_put_att(fid, varid_lon,  'long_name', 'longitude'))
call check( nf90_put_att(fid, varid_lat,  'long_name', 'latitude'))
call check( nf90_put_att(fid, varid_lon,  'units', 'degE [-180;180]'))
call check( nf90_put_att(fid, varid_lat,  'units', 'degN [-90;90]'))

! dimension description time
write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearnew, '-', 1, '-', 1, ' 0:0:0'

call check( nf90_put_att(fid, varid_time, 'long_name', 'time'))
call check( nf90_put_att(fid, varid_time, 'standard_name', 'time'))
call check( nf90_put_att(fid, varid_time, 'units', trim(att_text)))
call check( nf90_put_att(fid, varid_time, 'axis', 'T'))
call check( nf90_put_att(fid, varid_time, 'stored_direction', 'increasing'))

! forgetting factor
IF ((memb==0) .and. freq=='day') call check( nf90_def_var(fid, 'forget', NF90_FLOAT, dimID_time, varid_forget))

! fill dimension variables
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'write spatial coordinates'

call check (nf90_enddef(fid))
call check (nf90_put_var(fid, varid_nz,   mesh_fesom% Z   (1:nlmax)))
call check (nf90_put_var(fid, varid_nod2, [(n,n=1,mesh_fesom% nod2D)]))

call check (nf90_put_var(fid, varid_lon,  REAL(180./pi * lon, 4)))
call check (nf90_put_var(fid, varid_lat,  REAL(180./pi * lat, 4)))

if (debug) write (*,'(a, 10x,a,1x,i0)') 'FESOM-PDAF', 'fid is', fid
             
END SUBROUTINE netCDF_deffile
     
     
! ********************************
! ***   netCDF_defvar          ***
! ********************************
! This routine defines variables
! It is called from netCDF_init on writepe
! It closes the netCDF file

SUBROUTINE netCDF_defvar(memb)
! *** Arguments ***
INTEGER, intent (in)  :: memb  ! Zero:    Mean state
                               ! Number:  Ensemble member
! *** Local ***
LOGICAL :: defthis       ! True:  Define this field
LOGICAL :: writedaily    ! True:  Field is written at any day
                         ! False: Field is written at now_to_write_monthly

if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define output variables'

! switch to netCDF variable definition mode

IF ((memb==0) .and. w_dayensm) THEN
   call check (nf90_redef(fidphyday))
   call check (nf90_redef(fidbgcday))
ENDIF

IF ((memb >0) .and. w_daymemb) THEN
   call check (nf90_redef(fidphyday))
   call check (nf90_redef(fidbgcday))
ENDIF

IF ((memb==0) .and. w_monensm) THEN
   call check (nf90_redef(fidphymon))
   call check (nf90_redef(fidbgcmon))
ENDIF

IF ((memb >0) .and. w_monmemb) THEN
   call check (nf90_redef(fidphymon))
   call check (nf90_redef(fidbgcmon))
ENDIF

! LOOP: define state fields
DO j = 1, 8 ! ini / forc / ana / mean / ...
  DO i = 1, nfields ! state fields
     
    ! define this field?
    !         ___activated_____________       ___ens-mean/member_________________________
    defthis = (sfields(i)%output(j,oo)) .and. ((memb >0) .eqv. (sfields(i)%output(j,ee)))
    
    if (debug) write(*,'(a10,1x,a3,1x,i3,1x,l2)') sfields(i)%variable, IFA(j), memb, defthis
        
    ! yes, define this field:
    IF (defthis) THEN
     
         ! daily / monthly field?
         writedaily = sfields(i)%output(j,dd)
     
         ! choose file according to fields
         IF (      (sfields(i)%bgc) .and. writedaily) fileid = fidbgcday
         IF (.not. (sfields(i)%bgc) .and. writedaily) fileid = fidphyday
         IF (      (sfields(i)%bgc) .and. .not. writedaily) fileid = fidbgcmon
         IF (.not. (sfields(i)%bgc) .and. .not. writedaily) fileid = fidphymon
         
         ! set dimensions according to field before defining variables
         dimIDs(1) = dimID_nod2
         IF (sfields(i)% ndims == 1) THEN     ! surface fields
           dimIDs(2) = dimID_time
         ELSEIF (sfields(i)% ndims == 2) THEN ! 3D-fields
           dimIDs(2) = dimID_nz
           dimIDs(3) = dimID_time
         ENDIF
         
         ! number of dimensions
         IF (ANY(IFA(j)==(/'ii','si'/))) THEN
           ndims = sfields(i)% ndims     ! initial field: no iteration
         ELSE
           ndims = sfields(i)% ndims +1  ! plus iteration
         ENDIF
         
         ! define variable
         call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
         ! variable description
         call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
         call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))
         call check( nf90_put_att(fileid, sfields(i)% varid(j), '_FillValue', REAL(fill_value,4)))
         call check( nf90_put_att(fileid, sfields(i)% varid(j), 'CDI_grid_type', "unstructured"))
    
    ENDIF ! defthis
    
    IF (memb==0 .and. (ANY(IFA(j)==(/'ii','aa','ff'/)))) THEN
    ! daily ensemble standard deviation for each field
    IF (      (sfields(i)%bgc)) fileid = fidbgcday
    IF (.not. (sfields(i)%bgc)) fileid = fidphyday
    
    ! standard deviation global ocean volume averaged
    IF (IFA(j)=='ii') THEN
       ndims = 0
    ELSE
       ndims = 1
    ENDIF
    dimIDs(1) = dimID_time
    ! define variable
    call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j)//'_STD_global', NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
    ! variable description
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))//' global average of ensemble standard deviation'))
    call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))
    call check( nf90_put_att(fileid, sfields(i)% varid(j), '_FillValue', REAL(fill_value,4)))
    
    ! standard deviation layerwise
    IF (sfields(i)% ndims==2) THEN
      IF (IFA(j)=='ii') THEN
        ndims     = 1
      ELSE
        ndims     = 2
      ENDIF
      dimIDs(1) = dimID_nz
      dimIDs(2) = dimID_time
      ! define variable
      call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j)//'_STD_layerwise', NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
      ! variable description
      call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))//' layerwise average of ensemble standard deviation'))
      call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))
      call check( nf90_put_att(fileid, sfields(i)% varid(j), '_FillValue', REAL(fill_value,4)))
    ENDIF ! ndims==2
    ENDIF ! memb==0
    
  ENDDO ! state fields
ENDDO ! ini / forc / ana / ...

! close file
IF ((memb==0) .and. w_dayensm) THEN
   call check (nf90_enddef(fidphyday))
   call check (nf90_enddef(fidbgcday))
   call check (nf90_close (fidphyday))
   call check (nf90_close (fidbgcday))
ENDIF

IF ((memb >0) .and. w_daymemb) THEN
   call check (nf90_enddef(fidphyday))
   call check (nf90_enddef(fidbgcday))
   call check (nf90_close (fidphyday))
   call check (nf90_close (fidbgcday))
ENDIF

IF ((memb==0) .and. w_monensm) THEN
   call check (nf90_enddef(fidphymon))
   call check (nf90_enddef(fidbgcmon))
   call check (nf90_close (fidphymon))
   call check (nf90_close (fidbgcmon))
ENDIF

IF ((memb >0) .and. w_monmemb) THEN
   call check (nf90_enddef(fidphymon))
   call check (nf90_enddef(fidbgcmon))
   call check (nf90_close (fidphymon))
   call check (nf90_close (fidbgcmon))
ENDIF

END SUBROUTINE netCDF_defvar



! ********************************
! ***                          ***
! ***   netCDF_out             ***
! ***                          ***
! ********************************

SUBROUTINE netCDF_out(writetype, state_p, memb, now_to_write_monthly, stdev_surf_g, stdev_volo_g, forget, m_state_p)

USE g_clock, ONLY: daynew, month, day_in_month, fleapyear
USE recom_config, ONLY: SecondsPerDay

! ARGUMENTS:
CHARACTER(len=2), intent(in) :: writetype                 ! Write (i) initial, (a) assimilated, (f) forecast, (m) daily-average fields
REAL, INTENT(in)    :: state_p(dim_state_p)               ! State vector; can be ensemble-mean or member
INTEGER, INTENT(in) :: memb                               ! Number of ensemble-member or zero for ensemble-mean
LOGICAL, INTENT(in) :: now_to_write_monthly               ! Whether it's time to write monthly output
REAL, INTENT(in), optional :: stdev_surf_g(nfields,nlmax) ! Ensemble spread layerwise
REAL, INTENT(in), optional :: stdev_volo_g(nfields)       ! Ensemble spread ocean volume
REAL, INTENT(in), optional :: forget                      ! Forgetting factor; passed e.g. during forecast phase for ensemble-mean
REAL, INTENT(in), optional :: m_state_p(dim_state_p)      ! Monthly mean or snapshot; passed at the end of month

! Local variables:
REAL, allocatable :: myData2(:)                    ! Temporary array for pe-local surface fields
REAL, allocatable :: myData3(:,:)                  ! Temporary array for pe-local 3D-fields
REAL, allocatable :: data2_g(:)                    ! Temporary array for global surface fields
REAL, allocatable :: data3_g(:,:)                  ! Temporary array for global 3D-fields
LOGICAL :: writethisnow  ! True:  Write this field
LOGICAL :: writedaily    ! True:  Field is written at any day
                         ! False: Field is written at now_to_write_monthly
INTEGER :: writepos

! Print screen information:
IF (writepe) THEN
  IF (writetype == 'ii') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write initial ocean state to NetCDF at day ', &
          daynew, 'mens/memb', memb
  ELSE IF (writetype== 'ff') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write forecast to NetCDF at day ', &
          daynew, 'mens/memb', memb
  ELSE IF (writetype== 'aa') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write analysis to NetCDF at day ', &
          daynew, 'mens/memb', memb
  ELSE IF (writetype== 'mm') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write m-field to NetCDF at day ', &
          daynew, 'mens/memb', memb
  END IF
END IF ! writepe

IF (writetype == 'ii') j=ii
IF (writetype == 'aa') j=aa
IF (writetype == 'ff') j=ff
IF (writetype == 'mm') j=mm
IF (writetype == 'sf') j=sf
IF (writetype == 'sa') j=sa
IF (writetype == 'si') j=si
IF (writetype == 'sm') j=sm

IF (writepe) THEN
   ! Open netCDF files
   ! daily file
   IF ((memb==0) .and. w_dayensm) THEN
      CALL netCDF_openfile('phy','day',int0,fidphyday)
      CALL netCDF_openfile('bgc','day',int0,fidbgcday)
   ENDIF
   IF ((memb >0) .and. w_daymemb) THEN
      CALL netCDF_openfile('phy','day',memb,fidphyday)
      CALL netCDF_openfile('bgc','day',memb,fidbgcday)
   ENDIF
   ! monthly file
   IF ((memb==0) .and. w_monensm .and. now_to_write_monthly) THEN
      CALL netCDF_openfile('phy','mon',int0,fidphymon)
      CALL netCDF_openfile('bgc','mon',int0,fidbgcmon)
   ENDIF
   IF ((memb >0) .and. w_monmemb .and. now_to_write_monthly) THEN
      CALL netCDF_openfile('phy','mon',memb,fidphymon)
      CALL netCDF_openfile('bgc','mon',memb,fidbgcmon)
   ENDIF
END IF ! writepe

IF (writepe) THEN
    ! time-variable
    ! written during forecast phase
    IF (writetype=='ff') THEN
      
      ! write time coordinate for day
      IF (     ((memb==0) .and. w_dayensm) &
          .or. ((memb >0) .and. w_daymemb)) THEN
      call check( nf90_inq_varid(fidphyday, "time"  , varid_time))
      call check( nf90_put_var  (fidphyday, varid_time, (daynew-1)*SecondsPerDay, &
                                 start=(/ daynew /)))
      call check( nf90_inq_varid(fidbgcday, "time"  , varid_time))
      call check( nf90_put_var  (fidbgcday, varid_time, (daynew-1)*SecondsPerDay, &
                                 start=(/ daynew /)))
      
      ENDIF ! w_day
              
      ! write time coordinate for month
      IF (     ((memb==0) .and. w_monensm .and. now_to_write_monthly) &
          .or. ((memb >0) .and. w_monmemb .and. now_to_write_monthly)) THEN
          call check( nf90_inq_varid(fidphymon, "time",  varid_time))
          call check( nf90_put_var  (fidphymon, varid_time, (daynew+1-num_day_in_month(fleapyear,month))*SecondsPerDay, &
                                     start=(/ month /)))
          call check( nf90_inq_varid(fidbgcmon, "time",  varid_time))
          call check( nf90_put_var  (fidbgcmon, varid_time, (daynew+1-num_day_in_month(fleapyear,month))*SecondsPerDay, &
                                     start=(/ month /)))
          
      ENDIF ! w_mon
    ENDIF ! writetype=='ff'
    
    ! write forgetting factor
    ! recommended to pass daily for memb==0
    IF (present(forget)) THEN
       call check( nf90_inq_varid(fidphyday, "forget", varid_forget))
       call check( nf90_put_var  (fidphyday, varid_forget, REAL(forget,4), &
                                  start=(/ daynew /)))
       call check( nf90_inq_varid(fidbgcday, "forget", varid_forget))
       call check( nf90_put_var  (fidbgcday, varid_forget, REAL(forget,4), &
                                  start=(/ daynew /)))
    ENDIF ! present(forget)
  ENDIF ! writepe

  ! Field-specific:
  DO i=1, nfields
    
    writedaily = sfields(i)%output(j,dd)
    
    ! write this field?
    writethisnow =       (sfields(i)%output(j,oo)) &                    ! - output for field activated
                   .and. ((memb >0) .eqv. (sfields(i)%output(j,ee))) &  ! - memb-output to memb-file; ensm-output to ensm-file
                   .and. (writedaily .or. now_to_write_monthly)         ! - do it at the right time
    
    if (debug .and. writepe) write(*,'(a10,1x,a3,1x,i3,1x,l2)') sfields(i)%variable, IFA(j), memb, writethisnow
    
    ! yes, write this field:
    IF (writethisnow) THEN
    
      ! choose file according to fields
      IF (writepe) THEN
        IF (      (sfields(i)%bgc) .and. writedaily) fileid = fidbgcday
        IF (.not. (sfields(i)%bgc) .and. writedaily) fileid = fidphyday
        IF (      (sfields(i)%bgc) .and. .not. writedaily) fileid = fidbgcmon
        IF (.not. (sfields(i)%bgc) .and. .not. writedaily) fileid = fidphymon
      ENDIF ! writepe
      
      ! --------------
      ! surface fields
      ! --------------
      IF (sfields(i)% ndims == 1) THEN
         ! select pe-local field from state vector
         allocate(myData2(myDim_nod2d))
         DO n = 1, myDim_nod2D
           IF (writedaily) THEN
              myData2(n) = state_p(n + offset(i))
              writepos = daynew
           ELSE
              myData2(n) = m_state_p(n + offset(i))
              writepos = month
           ENDIF ! writedaily
         END DO ! n = 1, myDim_nod2D
         ! gather global field
         allocate(data2_g(mesh_fesom% nod2D))
         CALL gather_nod(myData2, data2_g)
         deallocate(myData2)
         
         IF (writepe) THEN
           ! Inquire variable ID
           call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(1)))
           ! Write variable to netCDF
           IF (ANY(writetype==(/'ii','si'/))) THEN   ! initial field
              call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(data2_g,4)))
           ELSE                                      ! forecast, analysis or m-fields
              call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(data2_g,4), &
                                   start=(/ 1, writepos /), &
                                   count=(/ mesh_fesom% nod2D, 1 /) ))
                                   ! (dims: 1-nod2, 2-time)
           ENDIF ! writetype
         ENDIF ! writepe
         deallocate(data2_g)
       
       ! ---------
       ! 3D fields
       ! ---------
       ELSEIF ((sfields(i)% ndims == 2)) THEN
         ! select pe-local field from state vector
         allocate(myData3(nlmax, myDim_nod2d))
         DO n = 1, myDim_nod2D
         DO l = 1, nlmax
           IF (writedaily) THEN
              myData3(l, n) = state_p((n-1) * nlmax + l + offset(i))
              writepos = daynew
           ELSE
              myData3(l, n) = m_state_p((n-1) * nlmax + l + offset(i))
              writepos = month
           ENDIF ! writedaily
         END DO
         END DO
         ! gather global field
         allocate(data3_g(nlmax,mesh_fesom% nod2D))
         CALL gather_nod(myData3, data3_g)
         WHERE (topography3D_g == 0) data3_g = fill_value
         deallocate(myData3)
         
         IF (writepe) THEN
           ! Inquire variable ID
           call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(1)))
           ! Write variable to netCDF:
           IF (ANY(writetype==(/'ii','si'/))) THEN   ! initial field
             call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(TRANSPOSE(data3_g),4)))
           ELSE                                      ! forecast, analysis and mean fields
             call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(TRANSPOSE(data3_g),4), &
                                      start=(/ 1, 1, writepos /), &
                                      count=(/ mesh_fesom% nod2D, nlmax, 1 /) ))
                                      ! dims: 1-nod2, 2-nz / nz1, 3-iter
           ENDIF ! writetype
         ENDIF ! writepe
         deallocate(data3_g)
         
       ENDIF ! surface / 3D-fields  
    ENDIF ! writethisnow
    
    ! ---------------------------
    !    Daily ensemble spread
    ! ---------------------------
    
    IF (present(stdev_volo_g) .and. writepe) THEN
       ! choose file
       IF (      (sfields(i)%bgc)) fileid = fidbgcday
       IF (.not. (sfields(i)%bgc)) fileid = fidphyday
       ! inquire variable ID
       call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype//'_STD_global', sfields(i)% varid(2)))
       ! write variable to netCDF
       IF (writetype=='ii') THEN
          call check( nf90_put_var(fileid, sfields(i)% varid(2), REAL(stdev_volo_g(i),4)))
       ENDIF ! writetype=='ii'
       IF (ANY(writetype==(/'ff','aa'/))) THEN
          call check( nf90_put_var(fileid, sfields(i)% varid(2), REAL(stdev_volo_g(i),4), &
                                      start=(/ daynew /)))
       ENDIF ! writetype==(/'ff','aa'/)
    ENDIF ! present(stdev_volo_g)
    
    IF (present(stdev_surf_g) .and. writepe .and. (sfields(i)%ndims==2)) THEN
       ! choose file
       IF (      (sfields(i)%bgc)) fileid = fidbgcday
       IF (.not. (sfields(i)%bgc)) fileid = fidphyday
       ! inquire variable ID
       call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype//'_STD_layerwise', sfields(i)% varid(3)))
       ! write variable to netCDF
       IF (writetype=='ii') THEN
          call check( nf90_put_var(fileid, sfields(i)% varid(3), REAL(stdev_surf_g(i,:),4)))
       ENDIF ! writetype=='ii'
       IF (ANY(writetype==(/'ff','aa'/))) THEN
          call check( nf90_put_var(fileid, sfields(i)% varid(3), REAL(stdev_surf_g(i,:),4), &
                                      start=(/     1, daynew /), &
                                      count=(/ nlmax,      1 /) ))
       ENDIF ! writetype==(/'ff','aa'/)
    ENDIF ! present(stdev_surf_g)
  ENDDO ! nfields
  
  ! Close file:
  IF (writepe) THEN
    ! daily file
    IF ((memb==0) .and. w_dayensm) THEN
       call check (nf90_close(fidphyday))
       call check (nf90_close(fidbgcday))
    ENDIF
    IF ((memb >0) .and. w_daymemb) THEN
       call check (nf90_close(fidphyday))
       call check (nf90_close(fidbgcday))
    ENDIF
    ! monthly file
    IF ((memb==0) .and. w_monensm .and. now_to_write_monthly) THEN
       call check (nf90_close(fidphymon))
       call check (nf90_close(fidbgcmon))
    ENDIF
    IF ((memb >0) .and. w_monmemb .and. now_to_write_monthly) THEN
       call check (nf90_close(fidphymon))
       call check (nf90_close(fidbgcmon))
    ENDIF
  ENDIF ! writepe
  

END SUBROUTINE netCDF_out

! ********************************
! ***   netCDF_openfile         ***
! ********************************
! This routine opens the netCDF output file into put-var mode
! It is called from netCDF_out on writepe

SUBROUTINE netCDF_openfile(typ,freq,memb,fid)
! *** Arguments ***
CHARACTER(len=3), intent (in)  :: typ   ! phy or bgc
CHARACTER(len=3), intent (in)  :: freq  ! day or mon
INTEGER,          intent (in)  :: memb    ! Zero:    Mean state
                                        ! Number:  Ensemble member
INTEGER,          intent (out) :: fid   ! file ID
              
IF (memb==0) THEN
  ! mean state
  filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.'//typ//'.'//cyearnew//'.'//freq//'.nc'
ELSE
  ! ensemble member
  WRITE(memberstr,'(i2.2)') memb
  filename = TRIM(DAoutput_path)//'fesom-recom-pdaf.'//typ//'.'//cyearnew//'.'//freq//'.'//memberstr//'.nc'
ENDIF

! open file
call check(NF90_OPEN(trim(filename),NF90_write,fid))
END SUBROUTINE netCDF_openfile

END MODULE mod_nc_out_routines
