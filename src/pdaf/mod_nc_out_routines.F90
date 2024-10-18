MODULE mod_nc_out_routines

! contains:
!   - check()
!   - netCDF_init()
!        -> calls netCDF_init_mm()
!   - netCDF_init_rmsens()
!   - netCDF_out()
!   - netCDF_STD_init()
!   - netCDF_STD_out()


! USES:
USE mod_nc_out_variables
USE mod_assim_pdaf, &
   ONLY: nfields, id, mesh_fesom, nlmax, DAoutput_path, dim_ens, dim_state, &
         dim_state_p, offset, phymin, phymax, bgcmin, bgcmax
USE mod_parallel_pdaf, &
   ONLY: abort_parallel, writepe
USE mod_obs_f_pdaf, &
   ONLY: pi
USE g_config, &
   ONLY: runid
USE g_clock, &
   ONLY: cyearnew, timeold, dayold, yearold, num_day_in_month, fleapyear
USE g_parsup, &
   ONLY: myDim_nod2D
USE g_comm_auto, &
   ONLY: gather_nod
USE netcdf

IMPLICIT NONE

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
CHARACTER, dimension(4) :: IFA              ! Type character ('i','a','f','m)
CHARACTER(len=8), dimension(4) :: IFA_long  ! Type characterG
INTEGER :: dimID_nod2, dimID_time, &
           dimID_nz, dimID_iter
INTEGER :: dimIDs(3)
INTEGER :: varid_nod2, varid_nz, &    ! netCDF IDs for dimensions
           varid_time, varid_iter, &              
           varid_lon, varid_lat, &
           varid_forget               ! netCDF ID for forgetting factor
INTEGER :: ndims
INTEGER :: nlay

REAL, allocatable :: lon(:)
REAL, allocatable :: lat(:)

INTEGER, parameter :: int0=0

LOGICAL :: debug = .true.

! SHOULD BE DELETED AT THE END:
LOGICAL :: DEBUGOUTPUT = .false.

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

INTEGER :: memb ! Zero:    Mean state
                ! Number:  Ensemble member
                
! initialize:
IFA(ff)='f' ! "f" (forecast)
IFA(aa)='a' ! "a" (analysis)
IFA(mm)='m' ! "m" (daily means)
IFA(ii)='i' ! "i" (initial)
IFA_long(ff)='forecast'
IFA_long(aa)='analysis'
IFA_long(mm)='daymean'
IFA_long(ii)='initial'

! gather GEO coordinates (from all PEs)
allocate(lon(mesh_fesom% nod2D),lat(mesh_fesom% nod2D))
call gather_nod(mesh_fesom%geo_coord_nod2D(1, 1:myDim_nod2D), lon)
call gather_nod(mesh_fesom%geo_coord_nod2D(2, 1:myDim_nod2D), lat)

! initialize file (on main PE)
IF (writepe) THEN
     
     DO memb=0,dim_ens
        ! _______________________________________
        ! (1) netCDF_deffile: create and open file, define dimensions and coordinates
        ! (2) netCDF_defvar:  define variables and close file
        IF ((memb==0) .and. w_dayensm) THEN
           CALL netCDF_deffile('phy','day',int0,fidphyday)
           CALL netCDF_deffile('bgc','day',int0,fidbgcday)
           CALL netCDF_defvar(int0)
        ENDIF
        
        IF ((memb >0) .and. w_daymemb) THEN
           CALL netCDF_deffile('phy','day',memb,fidphyday)
           CALL netCDF_deffile('bgc','day',memb,fidbgcday)
           CALL netCDF_defvar(memb)
        ENDIF
        
        IF ((memb==0) .and. w_monensm) THEN
           CALL netCDF_deffile('phy','mon',int0,fidphymon)
           CALL netCDF_deffile('bgc','mon',int0,fidbgcmon)
           CALL netCDF_defvar(int0)
        ENDIF
        
        IF ((memb >0) .and. w_monmemb) THEN
           CALL netCDF_deffile('phy','mon',memb,fidphymon)
           CALL netCDF_deffile('bgc','mon',memb,fidbgcmon)
           CALL netCDF_defvar(memb)
        ENDIF
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
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'Init output to file '//filename

call check(NF90_CREATE(trim(filename),NF90_NETCDF4,fid))

! define dimensions
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define dimensions'

call check( NF90_DEF_DIM(fid, 'nod2', mesh_fesom%nod2d, dimID_nod2))
call check( NF90_DEF_DIM(fid, 'nz',   nlmax,            dimID_nz))
call check( NF90_DEF_DIM(fid, 'time', NF90_UNLIMITED,   dimID_time))

! dimension variables
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define dimension variables'

call check( nf90_def_var(fid, 'nod2', NF90_INT,   dimID_nod2, varid_nod2))
call check( nf90_def_var(fid, 'nz',   NF90_FLOAT, dimID_nz,   varid_nz))
call check( nf90_def_var(fid, 'time', NF90_INT,   dimID_time, varid_time))

call check( nf90_def_var(fid, 'lon', NF90_FLOAT,   dimID_nod2, varid_lon))
call check( nf90_def_var(fid, 'lat', NF90_FLOAT,   dimID_nod2, varid_lat))

! dimension description
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'define dimension description'

call check( nf90_put_att(fid, varid_nod2, 'long_name', 'surface nodes'))
call check( nf90_put_att(fid, varid_nz,  'long_name', 'vertical layers (mid-layer depths)'))
call check( nf90_put_att(fid, varid_nz,  'units', 'm'))
call check( nf90_put_att(fid, varid_time, 'long_name', 'time'))

call check( nf90_put_att(fid, varid_lon,  'long_name', 'longitude'))
call check( nf90_put_att(fid, varid_lat,  'long_name', 'latitude'))
call check( nf90_put_att(fid, varid_lon,  'units', 'degE [-180;180]'))
call check( nf90_put_att(fid, varid_lat,  'units', 'degN [-90;90]'))

! dimension description time
write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'

call check( nf90_put_att(fid, varid_time, 'long_name', 'time'))
call check( nf90_put_att(fid, varid_time, 'standard_name', 'time'))
call check( nf90_put_att(fid, varid_time, 'units', trim(att_text)))
call check( nf90_put_att(fid, varid_time, 'axis', 'T'))
call check( nf90_put_att(fid, varid_time, 'stored_direction', 'increasing'))

! fill dimension variables
if (debug) write (*,'(a, 10x,a)') 'FESOM-PDAF', 'fill dimension variables'

call check (nf90_enddef(fid))
call check (nf90_put_var(fid, varid_nz,   mesh_fesom% Z   (1:nlmax)))
call check (nf90_put_var(fid, varid_nod2, [(n,n=1,mesh_fesom% nod2D)]))

call check (nf90_put_var(fid, varid_lon,  REAL(180./pi * lon, 4)))
call check (nf90_put_var(fid, varid_lat,  REAL(180./pi * lat, 4)))
             
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
   call check (nf90_redef(memb,fidphyday))
   call check (nf90_redef(memb,fidbgcday))
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
DO j = 1, 4 ! ini / forc / ana / mean
  DO i = 1, nfields ! state fields
     
    ! define this field?
    !         ___activated_____________       ___ens-mean/member_________________________
    defthis = (sfields(i)%output(j,oo)) .and. ((memb >0) .eqv. (sfields(i)%output(j,ee)))
    
    if (debug) write(*,'(a10,1x,a2,1x,i3,1x,l2)') sfields(i)%variable, IFA(j), memb, defthis
!~     if (debug) write(*,'(a25,1x,l2)') '(sfields(i)%output(j,oo))', (sfields(i)%output(j,oo))
!~     if (debug) write(*,'(a25,1x,l2)') '(memb >0)', (memb >0)
!~     if (debug) write(*,'(a25,1x,l2)') '(sfields(i)%output(j,ee))', (sfields(i)%output(j,ee))
    
    
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
         IF (IFA(j)=='i') THEN
           ndims = sfields(i)% ndims     ! initial field: no iteration
         ELSE
           ndims = sfields(i)% ndims +1  ! plus iteration
         ENDIF
         
         ! define variable
         call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
         ! variable description
         call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
         call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))
    
    ENDIF ! defthis
  ENDDO ! state fields
ENDDO ! ini / forc / ana

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

SUBROUTINE netCDF_out(writetype, state_p, memb, now_to_write_monthly, rms, forget, m_state_p)

USE g_clock, ONLY: daynew, month, day_in_month, fleapyear
USE recom_config, ONLY: SecondsPerDay

! ARGUMENTS:
CHARACTER(len=1), intent(in) :: writetype              ! Write (i) initial, (a) assimilated, (f) forecast, (m) daily-average fields
REAL, INTENT(in)    :: state_p(dim_state_p)            ! State vector; can be ensemble-mean or member
INTEGER, INTENT(in) :: memb                            ! Number of ensemble-member or zero for ensemble-mean
LOGICAL, INTENT(in) :: now_to_write_monthly            ! Whether it's time to write monthly output
REAL, INTENT(in), optional :: rms(nfields)             ! RMS ensemble spread; passed during analysis phase
REAL, INTENT(in), optional :: forget                   ! Forgetting factor; passed during analysis phase
REAL, INTENT(in), optional :: m_state_p(dim_state_p)   ! Monthly mean or snapshot; passed at the end of month

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
  IF (writetype == 'i') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write initial ocean state to NetCDF at day ', &
          daynew, 'mens/memb', memb
  ELSE IF (writetype== 'f') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write forecast to NetCDF at day ', &
          daynew, 'mens/memb', memb
  ELSE IF (writetype== 'a') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write analysis to NetCDF at day ', &
          daynew, 'mens/memb', memb
  ELSE IF (writetype== 'm') THEN
     WRITE (*, '(a, 8x, a, i4,1x,a,i3)') 'FESOM-PDAF', 'Write m-field to NetCDF at day ', &
          daynew, 'mens/memb', memb
  END IF
END IF ! writepe

IF (writetype == 'i') j=ii
IF (writetype == 'a') j=aa
IF (writetype == 'f') j=ff
IF (writetype == 'm') j=mm

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
    IF (writetype=='f') THEN
      
      ! write day
      IF (w_day) THEN
      call check( nf90_inq_varid(fidphyday, "day"  , varid_time))
      call check( nf90_put_var  (fidphyday, varid_time, (daynew-1)*SecondsPerDay, &
                                 start=(/ daynew /)))
      call check( nf90_inq_varid(fidbgcday, "day"  , varid_time))
      call check( nf90_put_var  (fidbgcday, varid_time, (daynew-1)*SecondsPerDay, &
                                 start=(/ daynew /)))
      ENDIF ! w_day
              
      ! write month
      IF (w_mon .and. now_to_write_monthly) THEN                          
          call check( nf90_inq_varid(fidphymon, "month",  varid_time))
          call check( nf90_put_var  (fidphymon, varid_time, (daynew+1-num_day_in_month(fleapyear,month))*SecondsPerDay, &
                                     start=(/ month /)))
          call check( nf90_inq_varid(fidbgcmon, "month",  varid_time))
          call check( nf90_put_var  (fidbgcmon, varid_time, (daynew+1-num_day_in_month(fleapyear,month))*SecondsPerDay, &
                                     start=(/ month /)))
      ENDIF ! w_mon
        
!~         ! write forgetting factor
!~         if (.not. present(forget)) WRITE (*, '(a, 8x, a)') 'FESOM-PDAF', 'Please pass forget argument!'
!~         call check( nf90_inq_varid(fileid_phy, "forget", varid_forget))
!~         call check( nf90_put_var  (fileid_phy, varid_forget, REAL(forget,4), &
!~                                    start=(/ daynew /)))
    ENDIF ! writetype
  ENDIF ! writepe

  ! Field-specific:
  DO i=1, nfields
    
    writedaily = sfields(i)%output(i,dd)
    
    ! write this field?
    writethisnow =       (sfields(i)%output(j,oo)) &                    ! - output for field activated
                   .and. ((memb >0) .eqv. (sfields(i)%output(j,ee))) &  ! - memb-output to memb-file; ensm-output to ensm-file
                   .and. (writedaily .or. now_to_write_monthly)         ! - do it at the right time
    
    if (debug) write(*,'(a10,1x,a2,1x,i3,1x,l2)') sfields(i)%variable, IFA(j), memb, writethisnow
    
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
         
         IF (writepe) THEN
           ! Inquire variable ID
           call check( nf90_inq_varid(fileid, TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(1)))
           ! Write variable to netCDF
           IF (writetype=='i') THEN ! initial field
              call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(data2_g,4)))
           ELSE                     ! forecast, analysis or m-fields
              call check( nf90_put_var(fileid, sfields(i)% varid(1), REAL(data2_g,4), &
                                   start=(/ 1, writepos /), &
                                   count=(/ mesh_fesom% nod2D, 1 /) ))
                                   ! (dims: 1-nod2, 2-time)
           ENDIF ! writetype
         ENDIF ! writepe
         deallocate(myData2, data2_g)
       
       ! ---------
       ! 3D fields
       ! ---------
       ELSEIF ((sfields(i)% ndims == 2)) THEN
         nlay = nlmax
         ! select pe-local field from state vector
         allocate(myData3(nlay, myDim_nod2d))
         DO n = 1, myDim_nod2D
         DO l = 1, nlay
           IF (writedaily) THEN
              myData3(l, n) = state_p((n-1) * nlay + l + offset(i))
              writepos = daynew
           ELSE
              myData3(l, n) = m_state_p((n-1) * nlay + l + offset(i))
              writepos = month
           ENDIF ! writedaily
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
    ENDIF ! writethisnow
    
!~     ! ---------------------
!~     ! RMS (Ensemble spread)
!~     ! ---------------------
!~     IF (member/=0) CYCLE
!~     IF (writetype=='m') CYCLE
    
!~     IF (writepe) THEN
!~       ! Inquire variable ID
!~       call check( nf90_inq_varid(fileid, 'rms_'//TRIM(sfields(i)% variable)//'_'//writetype, sfields(i)% varid(2)))
!~       ! Write RMS to netCDF
!~       call check( nf90_put_var(fileid, sfields(i)% varid(2), REAL(rms(i),4), &
!~                                start=(/ writepos /)))
!~     ENDIF ! writepe
!~     ENDIF ! don't write analysis of not-updated variables            

  ENDDO ! nfields
  
  ! Close file:
  IF (writepe) THEN
    ! daily file
    IF ((memb==0) .and. w_dayensm) THEN
       call check (nf90_close(fidphyday)
       call check (nf90_close(fidbgcday)
    ENDIF
    IF ((memb >0) .and. w_daymemb) THEN
       call check (nf90_close(fidphyday)
       call check (nf90_close(fidbgcday)
    ENDIF
    ! monthly file
    IF ((memb==0) .and. w_monensm .and. now_to_write_monthly) THEN
       call check (nf90_close(fidphymon)
       call check (nf90_close(fidbgcmon)
    ENDIF
    IF ((memb >0) .and. w_monmemb .and. now_to_write_monthly) THEN
       call check (nf90_close(fidphymon)
       call check (nf90_close(fidbgcmon)
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


!~ ! ********************************
!~ ! ***                          ***
!~ ! ***   netCDF_init_rmsens     ***
!~ ! ***                          ***
!~ ! ********************************

!~ ! Inits netCDF output for ensemble spread (root-mean-square deviation)

!~ SUBROUTINE netCDF_init_rmsens(filename, istart, iend)

!~ character(len=200), intent(in) :: filename ! Full name of output file
!~ integer, intent(in) :: istart              ! First field in state vector
!~ integer, intent(in) :: iend                ! Last field in state vector

!~ ! open file
!~ call check( nf90_open(TRIM(filename), nf90_write, fileid))

!~ ! inquire dimension ID
!~ call check( nf90_inq_dimid( fileid, "time", dimID_iter))
!~ call check( nf90_inq_dimid( fileid, "nod2", dimID_nod2))
!~ call check( nf90_inq_dimid( fileid, "nz",   dimID_nz))
!~ call check( nf90_inq_dimid( fileid, "nz1",  dimID_nz1))

!~ ! ********************
!~ ! define rms variables
!~ ! ********************
!~ DO j = 1, 3 ! ini / forc / ana
!~   DO i = istart, iend ! state fields
    
!~     ! don't write analysis of not-updated variables
!~     IF ((.not. DEBUGOUTPUT) .and. (j==3) .and. (.not. (sfields(i)% updated))) CYCLE
!~     ! set dimension
!~     IF (j==1) THEN
!~       ! initial RMS
!~       ! define dimensionless (1D) variable
!~     call check( NF90_DEF_VAR(fileid, 'rms_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, sfields(i)% varid(j)))
!~     ELSE
!~       ! forecast / analysis RMS
!~       ! define timeseries variable
!~     call check( NF90_DEF_VAR(fileid, 'rms_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimID_iter, sfields(i)% varid(j)))
!~     ENDIF
    
!~     ! variable description
!~     call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', 'Ensemble spread surface (RMS) for '//trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    
!~   ENDDO ! state fields
!~ ENDDO ! ini / ana / forc

!~ ! ***************************
!~ ! define daily mean variables
!~ ! ***************************

!~ j = 4 ! daily mean data
!~ DO i = istart, iend ! state fields

!~     ! set dimensions according to field before defining variables
!~     dimIDs(1) = dimID_nod2
!~     IF (sfields(i)% ndims == 1) THEN ! surface fields
!~       dimIDs(2) = dimID_iter
!~     ELSEIF (sfields(i)% ndims == 2) THEN ! 3D-fields
!~       dimIDs(2) = dimID_nz1
!~       IF (.not.(sfields(i)%nz1)) dimIDs(2) = dimID_nz ! variables on levels (vertical velocity)
!~       dimIDs(3) = dimID_iter
!~     ENDIF
    
!~     ! number of dimensions
!~     ndims = sfields(i)% ndims +1  ! plus iteration
    
!~     ! define variable
!~     call check( NF90_DEF_VAR(fileid, TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j)))
!~     ! variable description
!~     call check( nf90_put_att(fileid, sfields(i)% varid(j), 'long_name', trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
!~     call check( nf90_put_att(fileid, sfields(i)% varid(j), 'units',     sfields(i)% units))

!~   ENDDO ! state fields
  
!~ ! ***************************
!~ ! define forget
!~ ! ***************************

!~ ! define variable
!~ call check( NF90_DEF_VAR(fileid,'forget', NF90_FLOAT, dimID_iter, varID_forget))
!~ ! variable description
!~ call check( nf90_put_att(fileid, varID_forget, 'long_name', 'variable forgetting factor, snapshots'))


!~ ! close file
!~ call check (nf90_close(fileid))
!~ END SUBROUTINE netCDF_init_rmsens




! ********************************
! ***                          ***
! ***   netCDF_STD_init        ***
! ***                          ***
! ********************************
! Initializes a netCDF file for ensemble properties

SUBROUTINE netCDF_STD_init()

! *** Arguments ***
! *** Local variables ***
character(40) :: att_text
integer :: surf
integer :: prof
integer :: prme

surf=0  ! counters for variable IDs; sets of three for ini, forc and ana
prof=3
prme=6

! gather GEO coordinates (from all PEs)
allocate(lon(mesh_fesom% nod2D),lat(mesh_fesom% nod2D))
call gather_nod(mesh_fesom%geo_coord_nod2D(1, 1:myDim_nod2D), lon)
call gather_nod(mesh_fesom%geo_coord_nod2D(2, 1:myDim_nod2D), lat)

! initialize one global file on main PE
IF (writepe) THEN
  
filename_std = TRIM(DAoutput_path)//'fesom-recom-pdaf.STD.'//cyearnew//'.nc'
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initializing netCDF file: '//TRIM(filename_std)

! open file
call check(NF90_CREATE(trim(filename_std),NF90_NETCDF4,fileid))

! define dimensions
call check( NF90_DEF_DIM(fileid, 'time', NF90_UNLIMITED,   dimID_iter))
call check( NF90_DEF_DIM(fileid, 'nod2', mesh_fesom%nod2d, dimID_nod2))
call check( NF90_DEF_DIM(fileid, 'nz',   nlmax,            dimID_nz))

! dimension variables
call check( nf90_def_var(fileid, 'time', NF90_FLOAT, dimID_iter, varid_time))
call check( nf90_def_var(fileid, 'step', NF90_INT,   dimID_iter, varid_iter))
call check( nf90_def_var(fileid, 'nod2', NF90_INT,   dimID_nod2, varid_nod2))
call check( nf90_def_var(fileid, 'lon',  NF90_FLOAT, dimID_nod2, varid_lon))
call check( nf90_def_var(fileid, 'lat',  NF90_FLOAT, dimID_nod2, varid_lat))

call check( nf90_def_var(fileid, 'nz',   NF90_FLOAT, dimID_nz,   varid_nz))

! dimension description
write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearold, '-', 1, '-', 1, ' 0:0:0'
call check( nf90_put_att(fileid, varid_time, 'long_name', 'time'))
call check( nf90_put_att(fileid, varid_time, 'standard_name', 'time'))
call check( nf90_put_att(fileid, varid_time, 'units', trim(att_text)))
call check( nf90_put_att(fileid, varid_time, 'axis', 'T'))
call check( nf90_put_att(fileid, varid_time, 'stored_direction', 'increasing'))
call check( nf90_put_att(fileid, varid_iter, 'long_name', 'iteration'))

call check( nf90_put_att(fileid, varid_nod2, 'long_name', 'surface nodes'))
call check( nf90_put_att(fileid, varid_lon,  'long_name', 'longitude'))
call check( nf90_put_att(fileid, varid_lat,  'long_name', 'latitude'))
call check( nf90_put_att(fileid, varid_lon,  'units', 'degE [-180;180]'))
call check( nf90_put_att(fileid, varid_lat,  'units', 'degN [-90;90]'))

call check( nf90_put_att(fileid, varid_nz,   'long_name', 'vertical layers/levels'))
call check( nf90_put_att(fileid, varid_nz,   'units', 'm'))

! fill dimension variables
call check (nf90_enddef(fileid))
call check (nf90_put_var(fileid, varid_nz,   mesh_fesom% Z   (1:nlmax)))
call check (nf90_put_var(fileid, varid_nod2, [(n,n=1,mesh_fesom% nod2D)]))
call check (nf90_put_var(fileid, varid_lon,  REAL(180./pi * lon, 4)))
call check (nf90_put_var(fileid, varid_lat,  REAL(180./pi * lat, 4)))

! define field variables
call check (nf90_redef(fileid))

DO j = 1, 4 ! ini/forc/ana
  DO i = 1, nfields
 
    ! don't write analysis for not-updated variables
    ! IF ((DEBUGOUTPUT) .or. (j<3) .or. ((sfields(i)% updated))) THEN
    
    ! **********************
    ! *** surface STD fields
    ! **********************
    
    ! set dimensions
    dimIDs(1) = dimID_nod2
    dimIDs(2) = dimID_iter
    
    ! number of dimensions
    IF (IFA(j)=='i') THEN
      ndims = 1 ! no iteration
    ELSE
      ndims = 2 ! plus iteration
    ENDIF
    
    ! define surface variables
    call check( NF90_DEF_VAR(fileid, 'surf_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j+surf)))
    ! variable description
    call check( nf90_put_att(fileid, sfields(i)% varid(j+surf), 'long_name', 'Ensemble standard deviation at surface '//trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    call check( nf90_put_att(fileid, sfields(i)% varid(j+surf), 'units',      sfields(i)% units))
    
    ! ****************************
    ! *** profiles of mean and STD
    ! ****************************
    IF (sfields(i)% ndims == 2) THEN ! 3D-fields
    
    ! set dimensions
    dimIDs(1) = dimID_nz
    dimIDs(2) = dimID_iter
    
    ! number of dimensions
    IF (IFA(j)=='i') THEN
      ndims = 1 ! no iteration
    ELSE
      ndims = 2 ! plus iteration
    ENDIF
    
    ! define profile variables
    call check( NF90_DEF_VAR(fileid, 'prof_std_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j+prof)))
    call check( NF90_DEF_VAR(fileid, 'prof_gmn_'//TRIM(sfields(i)% variable)//'_'//IFA(j), NF90_FLOAT, dimIDs(1:ndims), sfields(i)% varid(j+prme)))
    ! variable description
    call check( nf90_put_att(fileid, sfields(i)% varid(j+prof), 'long_name', 'Global mean of ensemble standard deviation, on profile, '//trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    call check( nf90_put_att(fileid, sfields(i)% varid(j+prme), 'long_name', 'Global mean of ensemble mean of absolute state, on profile, '                   //trim(sfields(i)% long_name)//' '//trim(IFA_long(j))))
    call check( nf90_put_att(fileid, sfields(i)% varid(j+prof), 'units',      sfields(i)% units))
    call check( nf90_put_att(fileid, sfields(i)% varid(j+prme), 'units',      sfields(i)% units))
    
    ENDIF ! 3D-fields
    ! ENDIF ! don't write analysis for not-updated variables
  ENDDO ! i, nfields
ENDDO ! j, ini/forc/ana

call check(NF90_ENDDEF(fileid))
call check(NF90_CLOSE(fileid))

ENDIF ! writepe
deallocate(lat,lon)

END SUBROUTINE netCDF_STD_init

! ********************************
! ***                          ***
! ***   netCDF_STD_out         ***
! ***                          ***
! ********************************
! Write ensemble characteristics to netCDF file
SUBROUTINE netCDF_STD_out(writetype,writepos,iteration,stdev_p,stdevprof_g,meanprof_g)

! *** Arguments ***
CHARACTER(len=1), intent(in) :: writetype             ! Write (i) initial, (a) assimilated, (f) forecast, (m) daily-average fields
INTEGER, INTENT(in) :: writepos                       ! Write position
INTEGER, INTENT(in) :: iteration                      ! Current model time step
REAL, INTENT(in)    :: stdev_p(dim_state_p)           ! Ensemble standard deviation state vector
REAL, INTENT(in)    :: stdevprof_g(nlmax*nfields_3D)  ! Global profile of ensemble standard deviation
REAL, INTENT(in)    :: meanprof_g(nlmax*nfields_3D)   ! Global profile of ensemble mean

! *** Local variables ***
integer :: i3d                ! counter (3D-fields)
REAL    :: ctime              ! current time in seconds from the beginning of the year
integer :: surf, prof, prme   ! counters for variable IDs
REAL, allocatable :: myDataSurf(:)                    ! Temporary array for pe-local surface fields
REAL, allocatable :: dataSurf_g(:)                    ! Temporary array for global surface fields
REAL, allocatable :: dataprof(:)                      ! Temporary array for global profile fields

! counters for variable IDs; sets of three for ini, forc and ana
surf=0  
prof=3
prme=6

! current time in seconds from the beginning of the year
ctime=timeold+(dayold-1.)*86400

! Print screen information:
IF (writepe) THEN
  IF (writetype == 'i') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write initial ensemble properties to NetCDF at step ', &
          iteration, ' position ', writepos
  ELSE IF (writetype== 'f') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ensemble properties forecast to NetCDF at step ', &
          iteration, ' position ', writepos
  ELSE IF (writetype== 'a') THEN
     WRITE (*, '(a, 8x, a, i9, a, i5)') 'FESOM-PDAF', 'Write ensemble properties analysis to NetCDF at step ', &
          iteration, ' position ', writepos
  END IF
END IF ! writepe

! Open netCDF file:
filename_std = TRIM(DAoutput_path)//'fesom-recom-pdaf.STD.'//cyearnew//'.nc'
IF (writepe) THEN
  call check( nf90_open(TRIM(filename_std), nf90_write, fileid))
  
  ! Write time:
  IF (writetype=='a') THEN
      call check( nf90_inq_varid(fileid, "step", varid_iter))
      call check( nf90_inq_varid(fileid, "time", varid_time))
      
      call check( nf90_put_var  (fileid, varid_iter, iteration, &
                                 start=(/ writepos /)))
      call check( nf90_put_var  (fileid, varid_time, ctime, &
                                 start=(/ writepos /)))
  ENDIF ! writetype 'a'
ENDIF ! writepe
  
  i3D = 1
  DO i = 1, nfields
 
    ! don't write analysis for not-updated variables
    ! IF ((DEBUGOUTPUT) .or. (.not. (writetype=='a')) .or. ((sfields(i)% updated))) THEN
    
    ! **********************
    ! *** surface STD fields
    ! **********************
    ! get pe-local fields
    allocate(myDataSurf(myDim_nod2D))
    allocate(dataSurf_g(mesh_fesom%nod2D))
    
    IF (sfields(i)% ndims == 2) THEN ! 3D-fields
      DO n = 1, myDim_nod2D
             myDataSurf(n) = stdev_p((n-1) * (nlmax) + 1 + offset(i))
      ENDDO ! n, myDim_nod2D
    ELSEIF (sfields(i)% ndims == 1) THEN ! surface fields
      myDataSurf = stdev_p(offset(i)+1:offset(i)+myDim_nod2D)
    ENDIF ! 2D/3D-fields
    
    ! gather global surface fields
    CALL gather_nod(myDataSurf, dataSurf_g)
    
    ! write global field to netCDF
    IF (writepe) THEN
      ! Inquire variable ID
      call check( nf90_inq_varid(fileid, &
                  'surf_'//TRIM(sfields(i)% variable)//'_'//writetype, &
                  sfields(i)% varid(surf)))
      ! Write variable to netCDF:
      IF (writetype=='i') THEN ! initial field
        call check( nf90_put_var(fileid, sfields(i)% varid(surf), REAL(dataSurf_g,4)))
        
      ELSE                     ! forecast, analysis and mean fields
        call check( nf90_put_var(fileid, sfields(i)% varid(surf), REAL(dataSurf_g,4), &
                                 start=(/ 1, writepos /), &
                                 count=(/ mesh_fesom% nod2D, 1 /) ))
                                 ! (dims: 1-nod2, 2-time)
      ENDIF ! writetype (i,f,a,m)
    ENDIF ! writepe
    deallocate(myDataSurf,dataSurf_g)
    
    ! ****************************
    ! *** profiles of mean and STD
    ! ****************************
    IF (sfields(i)% ndims == 2) THEN ! 3D-fields
    IF (writepe) THEN
      ! *** STD ***
      allocate(dataprof(nlmax))
      dataprof = stdevprof_g( (i3D-1)*nlmax+1 : (i3D-1)*nlmax+nlmax )
      ! inquire variable ID
      call check( nf90_inq_varid(fileid, &
              'prof_std_'//TRIM(sfields(i)% variable)//'_'//writetype, &
              sfields(i)% varid(prof)))
      ! write variable to netCDF
      IF (writetype=='i') THEN ! initial field
        call check( nf90_put_var(fileid, sfields(i)% varid(prof), REAL(dataprof,4)))
      ELSE                     ! forecast, analysis and mean fields
        call check( nf90_put_var(fileid, sfields(i)% varid(prof), REAL(dataprof,4), &
                                 start=(/ 1, writepos /), &
                                 count=(/ nlmax, 1 /) ))
                                 ! (dims: 1-nz, 2-time)
      ENDIF ! writetype (i,f,a)
      deallocate(dataprof)
      ! *** Mean ***
      allocate(dataprof(nlmax))
      dataprof = meanprof_g( (i3D-1)*nlmax+1 : (i3D-1)*nlmax+nlmax )
      ! inquire variable ID
      call check( nf90_inq_varid(fileid, &
              'prof_gmn_'//TRIM(sfields(i)% variable)//'_'//writetype, &
              sfields(i)% varid(prme)))
      ! write variable to netCDF
      IF (writetype=='i') THEN ! initial field
        call check( nf90_put_var(fileid, sfields(i)% varid(prme), REAL(dataprof,4)))
      ELSE                     ! forecast, analysis and mean fields
        call check( nf90_put_var(fileid, sfields(i)% varid(prme), REAL(dataprof,4), &
                                 start=(/ 1, writepos /), &
                                 count=(/ nlmax, 1 /) ))
                                 ! (dims: 1-nz, 2-time)
      ENDIF ! writetype (i,f,a)
      deallocate(dataprof)
    ENDIF ! writepe
    i3D = i3D+1
    ENDIF ! 3D-fields
    ! ENDIF ! don't write analysis for not-updated variables
    ENDDO ! i, nfields

IF (writepe) call check (nf90_close(fileid))
END SUBROUTINE netCDF_STD_out


END MODULE mod_nc_out_routines
