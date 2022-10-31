! This program interpolates FESOM output to EN4 observation locations
PROGRAM interpolateEN4

! Use:
USE mpi

IMPLICIT NONE

INCLUDE 'netcdf.inc'

! Variables:
! ----------------------------------------------------------------------
!   G at the end of name:      global
!   no G:                      distributed on grid

INTEGER :: i,j,k,n, d                   ! counters
INTEGER :: dmax0, dmax                  ! number of days
INTEGER :: Nnodes                       ! global number of nodes
INTEGER :: nz0, nz                      ! (maximum) number of layers
INTEGER :: s                            ! auxiliary: status counter
INTEGER :: stat(500)                    ! auxiliary: status array

INTEGER :: year
CHARACTER(200) :: path                  ! Path to netCDF file
CHARACTER(200) :: runname
CHARACTER(200) :: filename              ! Name of netCDF file
CHARACTER(4)   :: yearstr               ! Year as str

INTEGER :: fileid                       ! ID of netCDF file

! To read FESOM-PDAF output:
INTEGER :: DimId_iter                    ! dimension: iteration
INTEGER :: DimId_n2D, DimId_nz           ! dimension: nodes, layers
INTEGER :: VarId_tempa, VarId_tempf      ! variable: temperature (analysis, forecast)
INTEGER :: VarId_salta, VarId_saltf      ! variable: salinity (analysis, forecast)
                                         ! Analysis/Forecast values for Temp/Salt
REAL(4), ALLOCATABLE :: t_anaG(:,:,:), t_ana(:,:,:)
REAL(4), ALLOCATABLE :: t_forG(:,:,:), t_for(:,:,:)
REAL(4), ALLOCATABLE :: s_anaG(:,:,:), s_ana(:,:,:)
REAL(4), ALLOCATABLE :: s_forG(:,:,:), s_for(:,:,:)

! To read nodal coordinates:
INTEGER :: varid_coordinate_n2d
REAL(4), ALLOCATABLE :: coordinate_n2dG(:,:) ! geo. coordinates [-Pi:Pi]
REAL(4), ALLOCATABLE :: coordinate_n2d(:,:)

! Node partitioning:
INTEGER :: npes                          ! number of model processes
INTEGER :: myDim_nod2D                   ! number of pe-local nodes
INTEGER :: eDim_nod2D                    ! number of pe-local halo nodes
INTEGER, ALLOCATABLE :: myList_nod2D(:)  ! pe-local node lists
INTEGER, ALLOCATABLE :: Nnodes_part(:)   ! rank partioning vector (global list of number of nodes per process)
TYPE :: nlistparttype                    ! holds gathered node lists on process 0
	INTEGER, ALLOCATABLE :: nodes(:)
END TYPE
TYPE(nlistparttype), DIMENSION(72) :: nlistpart

! To read EN4 observation info:
INTEGER, ALLOCATABLE :: n2d_temp(:,:), nl1_temp(:), &
                        n2d_salt(:,:), nl1_salt(:)                  ! nodes and layer of observations
INTEGER, ALLOCATABLE :: offset_temp(:), nobs_temp(:), &
                        offset_salt(:), nobs_salt(:)                ! daily offset and number of daily observations
REAL(4), ALLOCATABLE :: coordinate_temp(:,:), temperature(:)        ! geo. coordinates [-Pi:Pi]
REAL(4), ALLOCATABLE :: coordinate_salt(:,:), salinity(:)
REAL(4), ALLOCATABLE :: temperatureG(:), salinityG(:)
INTEGER :: varid_n2d_temp, varid_nl1_temp, varid_offset_temp
INTEGER :: varid_n2d_salt, varid_nl1_salt, varid_offset_salt
INTEGER :: varid_coordinate_temp, varid_temperature, varid_nobs_temp
INTEGER :: varid_coordinate_salt, varid_salinity,    varid_nobs_salt
INTEGER :: dimobs_temp, dimobs_salt                                 ! number of observations (in one year)
INTEGER :: dimobs_tempG, dimobs_saltG
INTEGER :: DimId_n_obs_temp, DimID_n_obs_sal
INTEGER :: ndays                                                    ! number of days in one year
INTEGER, ALLOCATABLE :: dimobs_temp_part(:), dimobs_salt_part(:)    ! rank partioning vector (global list of number of observations per process)
INTEGER, ALLOCATABLE :: dimobs_temp_dspl(:), dimobs_salt_dspl(:)    ! displacement (other word for offset)

! To calculate the interpolation:
REAL(4), ALLOCATABLE :: lon1(:), lon2(:), lon3(:)          ! geo. coordinates of 3 nodes that are associated with an observation
REAL(4), ALLOCATABLE :: lat1(:), lat2(:), lat3(:)
REAL(4), ALLOCATABLE :: t_ana1(:), t_ana2(:), t_ana3(:)    ! values at these 3 nodes
REAL(4), ALLOCATABLE :: t_for1(:), t_for2(:), t_for3(:)
REAL(4), ALLOCATABLE :: s_ana1(:), s_ana2(:), s_ana3(:)
REAL(4), ALLOCATABLE :: s_for1(:), s_for2(:), s_for3(:)
REAL(4), ALLOCATABLE :: t_ana_n(:), t_for_n(:), &
                        s_ana_n(:), s_for_n(:)             ! new interpolated values at observation location
REAL(4), ALLOCATABLE :: latT(:), lonT(:), latS(:), lonS(:) ! observation geo. coordinates for the day
INTEGER, ALLOCATABLE :: day_temp(:), day_salt(:)
REAL(4), ALLOCATABLE :: t_ana_nG(:), t_for_nG(:), &
                        s_ana_nG(:), s_for_nG(:)
REAL(4), ALLOCATABLE :: w1(:), w2(:), w3(:)                ! weights for the interpolation

! To gather interpolated values on process 0L
INTEGER, ALLOCATABLE :: day_tempG(:), day_saltG(:)
INTEGER, ALLOCATABLE :: nl1_tempG(:)
REAL(4), ALLOCATABLE :: lon_tempG(:)
REAL(4), ALLOCATABLE :: lat_tempG(:)

INTEGER :: offset1, offset2

! Parallel:
INTEGER        :: myPE
CHARACTER(4)   :: myPEstr4               ! myPE as str
CHARACTER(5)   :: myPEstr5               ! myPE as str
INTEGER        :: ierror, sizecomm, &
                  mpistatus(MPI_STATUS_SIZE), &
                  tag
                  
! Write to NetCDF:
INTEGER :: dimid_temp, dimid_salt, &
           varid_daytemp, varid_lontemp, varid_lattemp, varid_levtemp, &
           varid_daysalt, varid_lonsalt, varid_latsalt, varid_levsalt

! **********************************************************************
! *** Init:
! **********************************************************************
CALL MPI_INIT(ierror)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, myPE, IERROR)
CALL MPI_COMM_SIZE(mpi_comm_world, npes, ierror)

year  = 2016
ndays = 366

! **********************************************************************
! *** Read FESOM-PDAF output:
! **********************************************************************
! process 0 reads global file

IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Read FESOM-PDAF output:'
WRITE(*,*) '**********************************************************************'
END IF

IF (myPE==0) THEN ! read FESOM-PDAF output
s = 1

WRITE(yearstr,'(i4.4)') year
runname  ='FESOM_PDAF_jobscript' !'FESOM_PDAF_TEST72'
path     ='/work/ollie/frbunsen/model_runs/fesom2/'//TRIM(runname)//'/work/01/'
filename ='fesom.'//yearstr//'.oce.'//TRIM(runname)//'.nc'

WRITE(*,*) 'Reading FESOM-PDAF output from file: ', TRIM(path)//TRIM(filename)

! Open FESOM-PDAF NetCDF file:
stat(s) = NF_OPEN(TRIM(TRIM(path)//TRIM(filename)), NF_NOWRITE, fileid)

IF (stat(s) /= NF_NOERR) WRITE(*,*) 'NetCDF error opening FESOM-PDAF NetCDF'

! Inquire dimensions:
s = 1
stat(s) = NF_INQ_DIMID(fileid, 'iteration', DimId_iter)
s = s+1
stat(s) = NF_INQ_DIMLEN(fileid, DimId_iter,dmax0)
s = s+1
WRITE(*,*) 'Days of FESOM-PDAF output: ', dmax0
dmax=dmax0

stat(s) = NF_INQ_DIMID(fileid, 'nod2', DimId_n2D)
s = s+1
stat(s) = NF_INQ_DIMLEN(fileid, DimId_n2D, Nnodes)
s = s+1
WRITE(*,*) 'Global number of nodes: ', Nnodes

stat(s) = NF_INQ_DIMID(fileid, 'nz1', DimId_nz)
s = s+1
stat(s) = NF_INQ_DIMLEN(fileid, DimId_nz, nz0)
s = s+1
WRITE(*,*) 'Maximum depth levels: ', nz0
nz=nz0

! Inquire variables:
stat(s) = NF_INQ_VARID(fileid, 'temp_a', VarId_tempa)
s = s+1
stat(s) = NF_INQ_VARID(fileid, 'temp_f', VarId_tempf)
s = s+1
stat(s) = NF_INQ_VARID(fileid, 'salt_a', VarId_salta)
s = s+1
stat(s) = NF_INQ_VARID(fileid, 'salt_f', VarId_saltf)
s = s+1

ALLOCATE(t_anaG (Nnodes,nz0,dmax0))
ALLOCATE(t_forG (Nnodes,nz0,dmax0))
ALLOCATE(s_anaG (Nnodes,nz0,dmax0))
ALLOCATE(s_forG (Nnodes,nz0,dmax0))

! Read data:
stat(s) = NF_GET_VAR_REAL(fileid, VarId_tempa, t_anaG)
s = s+1
stat(s) = NF_GET_VAR_REAL(fileid, VarId_tempf, t_forG)
s = s+1
stat(s) = NF_GET_VAR_REAL(fileid, VarId_salta, s_anaG)
s = s+1
stat(s) = NF_GET_VAR_REAL(fileid, VarId_saltf, s_forG)
s = s+1
WRITE(*,*) 'debug t_anaG(200,:,1): ', t_anaG(200,:,1)

! Close FESOM-PDAF NetCDF file:  
stat(s) = NF_CLOSE(fileid)
s = s+1

DO i = 1,  s-1
   IF (stat(i) /= NF_NOERR) &
		WRITE(*,*) 'NetCDF error in reading from FESOM-PDAF NetCDF, no.', i
END DO
END IF ! read FESOM-PDAF output

CALL MPI_SCATTER(dmax0,1,mpi_integer,&
                 dmax, 1,mpi_integer,0,MPI_COMM_WORLD,ierror)
CALL MPI_SCATTER(nz0,  1,mpi_integer,&
                 nz,   1,mpi_integer,0,MPI_COMM_WORLD,ierror)
WRITE(*,*) 'debug ierror: ', ierror
WRITE(*,*) 'debug dmax: ', dmax
WRITE(*,*) 'debug nz: ', nz


! **********************************************************************
! *** Read global geographic nodal coordinates:
! **********************************************************************
! process 0 reads global file

IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Read nodal geo coordinates:'
WRITE(*,*) '**********************************************************************'
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

IF (myPE==0) THEN ! read nodal coordinates
s = 1

path     ='/work/ollie/frbunsen/model_runs/fesom2/helpfiles/'
filename ='fesom.mesh.diag.nc'

WRITE(*,*) 'Reading nodal coordinates from file: ', TRIM(path)//TRIM(filename)

! Open mesh diagnostics file:
stat(s) = NF_OPEN(TRIM(TRIM(path)//TRIM(filename)), NF_NOWRITE, fileid)

IF (stat(s) /= NF_NOERR) WRITE(*,*) 'NetCDF error opening mesh diagnostics'
s = 1

ALLOCATE(coordinate_n2dG(Nnodes,2))

stat(s) = NF_INQ_VARID(fileid, 'nodes', varid_coordinate_n2d)
s = s+1
stat(s) = NF_GET_VAR_REAL(fileid, varid_coordinate_n2d, coordinate_n2dG)

! Check whether values look ok:
WRITE(*,*) 'debug coordinate_n2dG(5:10,2): ', coordinate_n2dG(5:10,2)

stat(s) = NF_CLOSE(fileid)

END IF ! read nodal coordinates


! **********************************************************************
! *** Read mesh partitioning:
! **********************************************************************
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Read mesh partitioning:'
WRITE(*,*) '**********************************************************************'
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! File:
path     = '/work/ollie/projects/clidyn/FESOM2/meshes/core2/dist_72/'
write(myPEstr5,'(i5.5)') myPE
filename = trim(path)//'my_list'//trim(myPEstr5)//'.out'
fileID   = myPE+9                     ! any number 9-99, but unique for each file
IF (myPE==0) WRITE(*,*) 'Reading mesh partitioning from file: ', TRIM(filename)

! Reading mesh partitioning
! Every process reads its file
open(fileID, file=trim(filename))
read(fileID,*) n                      ! myPE
 
read(fileID,*) myDim_nod2D            ! number of vertices that belong to myPE
read(fileID,*) eDim_nod2D             ! number of halo vertices (share a common
                                      ! triangle with vertices that belong to myPE,
                                      ! yet do not belong to myPE themselves.
allocate(myList_nod2D(myDim_nod2D))
read(fileID,*) myList_nod2D           ! list of vertices
close(fileID)

! Check whether data is correct:
IF (myPE==0) WRITE(*,*) 'debug myList_nod2D(1:5): ',    myList_nod2D(1:5)
IF (myPE==0) WRITE(*,*) 'debug myList_nod2D(-5:end): ', myList_nod2D(myDim_nod2D-5:myDim_nod2D)

!~ ! Gather number of nodes on process 0:
IF (myPE==0) ALLOCATE(Nnodes_part(npes))
CALL MPI_GATHER(myDim_nod2D,1,MPI_INT,&
                Nnodes_part,1,MPI_INT,&
                0,MPI_COMM_WORLD,ierror)
! Check whether data is correct:
IF (myPE==0) WRITE(*,*) 'debug Nnodes_part: ', Nnodes_part

! Gather node lists on main:
! 1. Allocate space for all node lists in nlistpart on process 0.
IF (myPE == 0) THEN; DO i = 1, npes
ALLOCATE(nlistpart(i)%nodes(Nnodes_part(i)))
END DO; END IF
! 2. Send myList_nod2D from all other processes to process 0.
IF (myPE/=0) THEN
WRITE(*,*) '1 debug mype: ', mype
WRITE(*,*) '1 debug size(myList_nod2D): ', size(myList_nod2D)
WRITE(*,*) '1 debug myList_nod2D(1:5): ', myList_nod2D(1:5)
WRITE(*,*) '1 debug myDim_nod2D: ', myDim_nod2D
CALL MPI_SEND(myList_nod2D,myDim_nod2D,MPI_INT,&
              0, myPE, MPI_COMM_WORLD,ierror)
! 3. Gather myList_nod2D in nlistpart on 0.
ELSE
nlistpart(1)%nodes=myList_nod2D
DO i = 1, npes-1
WRITE(*,*) '2 debug mype: ', mype
WRITE(*,*) '2 debug size(nlistpart(i+1)%nodes): ', size(nlistpart(i+1)%nodes)
WRITE(*,*) '2 debug Nnodes_part(i+1): ', Nnodes_part(i+1)
CALL MPI_RECV(nlistpart(i+1)%nodes,Nnodes_part(i+1),MPI_INT,&
              i, i, MPI_COMM_WORLD, mpistatus, ierror)
END DO; END IF

! **********************************************************************
! *** Distribute analysis, forecast and model grid coordinates
! **********************************************************************
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Distribute analysis, forecast and model grid coordinates:'
WRITE(*,*) '**********************************************************************'
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! allocate space on every process
ALLOCATE(t_ana(myDim_nod2D,nz,dmax), &
         t_for(myDim_nod2D,nz,dmax)) !, &
!~          s_ana(myDim_nod2D,nz,dmax), &
!~          s_for(myDim_nod2D,nz,dmax))

ALLOCATE(coordinate_n2d(myDim_nod2D,2))

! 1. Send from process 0 to all other.
IF (myPE==0) THEN
t_ana = t_anaG(nlistpart(1)%nodes,:,:)
t_for = t_forG(nlistpart(1)%nodes,:,:)
coordinate_n2d = coordinate_n2dG(nlistpart(1)%nodes,:)
!!!!!!!!!!! SHOULD BE (i = 1, npes-1) IN THE NEXT LINE BUT THIS GIVES MEMORY ERROR. !!!!!!!!!!!!
DO i = 1, 1 ! npes-1
WRITE(*,*) '1 debug npes-1', npes-1
WRITE(*,*) '1 debug i: ', i
WRITE(*,*) '1 debug size(t_anaG(nlistpart(i+1)%nodes,:,:)): ', size(t_anaG(nlistpart(i+1)%nodes,:,:))
WRITE(*,*) '1 debug Nnodes_part(i+1)*nz*dmax: ', Nnodes_part(i+1)*nz*dmax
WRITE(*,*) '1 debug send dest', i
WRITE(*,*) '1 debug send tag', i
CALL MPI_SEND(t_anaG(nlistpart(i+1)%nodes,:,:), Nnodes_part(i+1)*nz*dmax, MPI_REAL4, &
              i, i, MPI_COMM_WORLD, ierror)
WRITE(*,*) '1 debug send ierror', ierror
CALL MPI_SEND(t_forG(nlistpart(i+1)%nodes,:,:), Nnodes_part(i+1)*nz*dmax, MPI_REAL4, &
              i, i+100, MPI_COMM_WORLD, ierror)
CALL MPI_SEND(coordinate_n2dG(nlistpart(i+1)%nodes,:), Nnodes_part(i+1)*2, MPI_REAL4, &
              i, i+200, MPI_COMM_WORLD, ierror)
END DO
DEALLOCATE(t_anaG, t_forG, s_anaG, s_forG, coordinate_n2dG)
! 2. Receive on every other process.
ELSE
WRITE(*,*) '2 debug myPE: ', myPE
WRITE(*,*) '2 debug size(t_ana): ', size(t_ana)
WRITE(*,*) '2 debug myDim_nod2D*nz*dmax: ', myDim_nod2D*nz*dmax
CALL MPI_RECV(t_ana, myDim_nod2D*nz*dmax, MPI_REAL4, &
              0, myPE, MPI_COMM_WORLD, mpistatus, ierror)
WRITE(*,*) '2 debug recv ierror', ierror
CALL MPI_RECV(t_for, myDim_nod2D*nz*dmax, MPI_REAL4, &
              0, myPE+100, MPI_COMM_WORLD, mpistatus, ierror)
CALL MPI_RECV(coordinate_n2d, myDim_nod2D*2, MPI_REAL4, &
              0, myPE+200, MPI_COMM_WORLD, mpistatus, ierror)
END IF

WRITE(*,*) 'debug myPE: ', myPE, 't_ana(200,:,1): ', t_ana(200,:,1)

! **********************************************************************
! *** Read EN4 observation info from distributed files:
! **********************************************************************
! Every process reads its file
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Read EN4 from distributed files:'
WRITE(*,*) '**********************************************************************'
END IF

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
s = 1

WRITE(myPEstr4,'(i4.4)') myPE

path     = '/work/ollie/frbunsen/data/physics/Temp-Salt-3D/EN4/fesom-pdaf-processed_dist_72/2016/'
filename = 'obs_profile_'//myPEstr4//'.nc'

IF (myPE==0) WRITE(*,*) 'Reading EN4 observation info from file: ', TRIM(path)//TRIM(filename)

! Open EN4 NetCDF file:
stat(s) = NF_OPEN(TRIM(TRIM(path)//TRIM(filename)), NF_NOWRITE, fileid)
IF (stat(s) /= NF_NOERR) WRITE(*,*) 'NetCDF error opening EN4 NetCDF'

! Inquire dimensions:
s = 1
stat(s) = NF_INQ_DIMID(fileid, 'n_obs_temp', DimId_n_obs_temp)
s = s+1
stat(s) = NF_INQ_DIMLEN(fileid, DimId_n_obs_temp,dimobs_temp)
s = s+1
IF (myPE==0) WRITE(*,*) 'Number of temperature observations: ', dimobs_temp

stat(s) = NF_INQ_DIMID(fileid, 'n_obs_sal', DimId_n_obs_sal)
s = s+1
stat(s) = NF_INQ_DIMLEN(fileid, DimId_n_obs_sal,dimobs_salt)
s = s+1
IF (myPE==0) WRITE(*,*) 'Number of salinity observations: ', dimobs_salt

DO i = 1,  s-1
   IF (stat(i) /= NF_NOERR) &
		WRITE(*,*) 'NetCDF error in reading from EN4 NetCDF, no.', i
END DO

IF (dimobs_temp /= 0) THEN ! have observations (temp)
	
	stat(s) = NF_INQ_VARID(fileid, 'n2d_temp',        varid_n2d_temp)
	s = s+1
	stat(s) = NF_INQ_VARID(fileid, 'offset_temp',     varid_offset_temp)
	s = s+1
	stat(s) = NF_INQ_VARID(fileid, 'coordinate_temp', varid_coordinate_temp)
	s = s+1
	stat(s) = NF_INQ_VARID(fileid, 'temperature',     varid_temperature)
	s = s+1
	stat(s) = NF_INQ_VARID(fileid, 'nl1_temp',        varid_nl1_temp)
	s = s+1
	stat(s) = NF_INQ_VARID(fileid, 'nobs_temp',       varid_nobs_temp)
	s = s+1
	
	ALLOCATE(n2d_temp(3,dimobs_temp))
	ALLOCATE(offset_temp(ndays))
	ALLOCATE(coordinate_temp(2,dimobs_temp))
	ALLOCATE(temperature(dimobs_temp))
	ALLOCATE(nl1_temp(dimobs_temp))
	ALLOCATE(nobs_temp(ndays))
	         
	stat(s) = NF_GET_VAR_INT(fileid, varid_n2d_temp, n2d_temp)
	s = s+1
	stat(s) = NF_GET_VAR_INT(fileid, varid_offset_temp, offset_temp)
	s = s+1
	stat(s) = NF_GET_VAR_REAL(fileid, varid_coordinate_temp, coordinate_temp)
	s = s+1
	stat(s) = NF_GET_VAR_REAL(fileid, varid_temperature, temperature)
	s = s+1
	stat(s) = NF_GET_VAR_INT(fileid, varid_nl1_temp, nl1_temp)
	s = s+1
	stat(s) = NF_GET_VAR_INT(fileid, varid_nobs_temp, nobs_temp)
	s = s+1
	
	DO i = 1,  s-1
    IF (stat(i) /= NF_NOERR) &
		WRITE(*,*) 'NetCDF error in reading EN4 temperature, no.', i
	END DO
	
	! Check whether data looks good:
	IF (myPe==0) THEN
	WRITE(*,*) 'debug n2d_temp(3,61:71): ',        n2d_temp(3,61:71)
	WRITE(*,*) 'debug offset_temp(6:11): ',        offset_temp(6:11)
	WRITE(*,*) 'debug coordinate_temp(2,61:71): ', coordinate_temp(2,61:71)
	WRITE(*,*) 'debug temperature(6:11): ',        temperature(6:11)
	WRITE(*,*) 'debug nl1_temp(6:11): ',           nl1_temp(6:11)
	WRITE(*,*) 'debug nobs_temp(6:11): ',          nobs_temp(6:11)
	END IF
	
END IF ! have observations (temp)

IF (dimobs_salt /= 0) THEN ! have observations (salt)
! TO-DO: The same for salt.
! n2d_salt, offset_salt, coordinate_salt, salinity, nl1_salt
END IF ! have observations (salt)

! Close EN4 NetCDF file:  
stat(s) = NF_CLOSE(fileid)
s = s+1

! **********************************************************************
! *** Interpolate the model output to observation coordinates:
! **********************************************************************
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Interpolate the model output to observation coordinates:'
WRITE(*,*) '**********************************************************************'
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

IF (dimobs_temp/=0) THEN ! have observations

! allocate array of days that will serve as time coordinate later
ALLOCATE(day_temp(dimobs_temp))

ALLOCATE(t_ana_n(dimobs_temp))
ALLOCATE(t_for_n(dimobs_temp))

DO d=1,dmax ! do for each day (d)
	
	! Select daily data by observation day offset:
	offset1 = offset_temp(d)
	offset2 = offset_temp(d)+nobs_temp(d)
	
	ALLOCATE(lonT(nobs_temp(d)),latT(nobs_temp(d)))
	
	lonT = coordinate_temp(1,offset1:offset2)
	latT = coordinate_temp(2,offset1:offset2)
	
	! Select corresponding model output by nodes and levels (given with the observations):
	ALLOCATE(t_ana1(nobs_temp(d)),&
	         t_ana2(nobs_temp(d)),&
	         t_ana3(nobs_temp(d)),&
	         t_for1(nobs_temp(d)),&
	         t_for2(nobs_temp(d)),&
	         t_for3(nobs_temp(d)),&
	         lon1(nobs_temp(d)),&
	         lon2(nobs_temp(d)),&
	         lon3(nobs_temp(d)),&
	         lat1(nobs_temp(d)),&
	         lat2(nobs_temp(d)),&
	         lat3(nobs_temp(d)))
	
	t_ana1 = pack(  t_ana( n2d_temp(1,offset1:offset2), nl1_temp(offset1:offset2), d ), .true.)
	t_ana2 = pack(  t_ana( n2d_temp(2,offset1:offset2), nl1_temp(offset1:offset2), d ), .true.)
	t_ana3 = pack(  t_ana( n2d_temp(3,offset1:offset2), nl1_temp(offset1:offset2), d ), .true.)
	
	t_for1 = pack(  t_ana( n2d_temp(1,offset1:offset2), nl1_temp(offset1:offset2), d ), .true.)
	t_for2 = pack(  t_ana( n2d_temp(2,offset1:offset2), nl1_temp(offset1:offset2), d ), .true.)
	t_for3 = pack(  t_ana( n2d_temp(3,offset1:offset2), nl1_temp(offset1:offset2), d ), .true.)
	
	lon1 = pack(  coordinate_n2d( n2d_temp(1,offset1:offset2) , 1 ), .true.)
	lon2 = pack(  coordinate_n2d( n2d_temp(2,offset1:offset2) , 1 ), .true.)
	lon3 = pack(  coordinate_n2d( n2d_temp(3,offset1:offset2) , 1), .true.)
	
	lat1 = pack(  coordinate_n2d( n2d_temp(1,offset1:offset2) , 2 ), .true.)
	lat2 = pack(  coordinate_n2d( n2d_temp(2,offset1:offset2) , 2 ), .true.)
	lat3 = pack(  coordinate_n2d( n2d_temp(3,offset1:offset2) , 2 ), .true.)
	
	! calculate weights for barycentric interpolation:
	ALLOCATE(w1(nobs_temp(d)), w2(nobs_temp(d)), w3(nobs_temp(d)))
	
	w1 =   ((lat2-lat3) * (lonT-lon3) + (lon3-lon2) * (latT-lat3)) &
         / ((lat2-lat3) * (lon1-lon3) + (lon3-lon2) * (lat1-lat3))
    
    w2 =   ((lat3-lat1) * (lonT-lon3) + (lon1-lon3) * (latT-lat3)) &
         / ((lat2-lat3) * (lon1-lon3) + (lon3-lon2) * (lat1-lat3))
    
    w3 = 1 - w1 - w2
    
    ! interpolate / weighted average:
    t_ana_n (offset1:offset2) = t_ana1*w1 + t_ana2*w2 + t_ana3*w3
    t_for_n (offset1:offset2) = t_for1*w1 + t_for2*w2 + t_for3*w3
    
    IF (d==3) THEN
    WRITE(*,*) 'debug myPE: ', myPE, 't_ana1(1:5): ', t_ana1(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 't_ana_n(1:5): ', t_ana_n(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 't_for1(1:5): ', t_for1(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 't_for_n(1:5): ', t_for_n(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 'lon1(1:5): ', lon1(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 'lat1(1:5): ', lat1(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 'w1(1:5): ', w1(1:5)
    WRITE(*,*) 'debug myPE: ', myPE, 'latT(1:5): ', latT(1:5)
    END IF

	
	DEALLOCATE(t_ana1, t_ana2, t_ana3)
	DEALLOCATE(t_for1, t_for2, t_for3)
	DEALLOCATE(lon1, lon2, lon3)
	DEALLOCATE(lat1, lat2, lat3)
	DEALLOCATE(w1, w2, w3)
	DEALLOCATE(latT,lonT)
	
	! fill array of days that will serve as time coordinate later
	IF ( myPE==0) WRITE(*,*) 'debug offset1: ', offset1
	IF ( myPE==0) WRITE(*,*) 'debug offset2: ', offset2
	
	day_temp(offset1:offset2) = d
	
END DO ! do for each day (d)
WRITE(*,*) 'debug myPE: ', myPE, 'day_temp(1:5): ', day_temp(1:5)
END IF ! have observations


! **********************************************************************
! *** Gather the interpolated values on process 0:
! **********************************************************************
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Gather interpolated values on process 0:'
WRITE(*,*) '**********************************************************************'
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

! Global number of observations (dimobs_tempG):
CALL MPI_REDUCE(dimobs_temp, dimobs_tempG, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
IF (myPE==0) WRITE(*,*) 'debug dimobs_tempG: ', dimobs_tempG

! Gather number of observations per process:
ALLOCATE(dimobs_temp_part(npes))
CALL MPI_GATHER(dimobs_temp,         1, MPI_INT,       &
                dimobs_temp_part,    1, MPI_INT, 0,    &
                MPI_COMM_WORLD, ierror)
IF (myPE==0) WRITE(*,*) 'debug dimobs_temp_part: ', dimobs_temp_part

IF (myPE==0) THEN
ALLOCATE(dimobs_temp_dspl(npes))
dimobs_temp_dspl(1) = 0
IF (npes>1) THEN
	DO i=2, npes
	dimobs_temp_dspl(i) = SUM(dimobs_temp_part(1:i-1))
	END DO
END IF
WRITE(*,*) 'debug dimobs_temp_dspl: ', dimobs_temp_dspl
END IF

! Gather analysis, forecast, observations, day, depth level, lon and lat:
! (Gather analysis temperature)
IF (myPE==0) ALLOCATE(t_ana_nG(dimobs_tempG))
WRITE(*,*) 'debug myPE: ', myPE, 't_ana_n(1:4) ', t_ana_n(1:4)
CALL MPI_GATHERV(t_ana_n,  dimobs_temp,                        MPI_REAL,       &
                 t_ana_nG, dimobs_temp_part, dimobs_temp_dspl, MPI_REAL, 0,    & 
                 MPI_COMM_WORLD, ierror)
DEALLOCATE(t_ana_n)
IF (myPE==0) WRITE(*,*) 'debug t_ana_nG(1:4) ', t_ana_nG(1:4)
IF (myPE==0) WRITE(*,*) 'debug t_ana_nG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', t_ana_nG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! (Gather forecast temperature)
IF (myPE==0)  ALLOCATE(t_for_nG(dimobs_tempG))
CALL MPI_GATHERV(t_for_n,  dimobs_temp,                        MPI_REAL,       &
                 t_for_nG, dimobs_temp_part, dimobs_temp_dspl, MPI_REAL, 0,    & 
                 MPI_COMM_WORLD, ierror)
DEALLOCATE(t_for_n)
IF (myPE==0) WRITE(*,*) 'debug t_for_nG(1:4) ', t_for_nG(1:4)
IF (myPE==0) WRITE(*,*) 'debug t_for_nG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', t_for_nG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! (Gather observations temperature)
IF (myPE==0)  ALLOCATE(temperatureG(dimobs_tempG))
CALL MPI_GATHERV(temperature,  dimobs_temp,                        MPI_REAL,       &
                 temperatureG, dimobs_temp_part, dimobs_temp_dspl, MPI_REAL, 0,    & 
                 MPI_COMM_WORLD, ierror)
DEALLOCATE(temperature)
IF (myPE==0) WRITE(*,*) 'debug temperatureG(1:4) ', temperatureG(1:4)
IF (myPE==0) WRITE(*,*) 'debug temperatureG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', temperatureG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! (Gather day)
IF (myPE==0)  ALLOCATE(day_tempG(dimobs_tempG))
CALL MPI_GATHERV(day_temp,  dimobs_temp,                        MPI_INT,       &
                 day_tempG, dimobs_temp_part, dimobs_temp_dspl, MPI_INT, 0,    & 
                 MPI_COMM_WORLD, ierror)
DEALLOCATE(day_temp)
IF (myPE==0) WRITE(*,*) 'debug day_tempG(1:4) ', day_tempG(1:4)
IF (myPE==0) WRITE(*,*) 'debug day_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', day_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! (Gather layer)
IF (myPE==0)  ALLOCATE(nl1_tempG(dimobs_tempG))
CALL MPI_GATHERV(nl1_temp,  dimobs_temp,                        MPI_REAL,       &
                 nl1_tempG, dimobs_temp_part, dimobs_temp_dspl, MPI_REAL, 0,    & 
                 MPI_COMM_WORLD, ierror)
DEALLOCATE(nl1_temp)
IF (myPE==0) WRITE(*,*) 'debug nl1_tempG(1:4) ', nl1_tempG(1:4)
IF (myPE==0) WRITE(*,*) 'debug nl1_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', nl1_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! (Gather longitude)
IF (myPE==0)  ALLOCATE(lon_tempG(dimobs_tempG))
CALL MPI_GATHERV(coordinate_temp(1,:), dimobs_temp,                        MPI_REAL,       &
                 lon_tempG,            dimobs_temp_part, dimobs_temp_dspl, MPI_REAL, 0,    & 
                 MPI_COMM_WORLD, ierror)
IF (myPE==0) WRITE(*,*) 'debug lon_tempG(1:4) ', lon_tempG(1:4)
IF (myPE==0) WRITE(*,*) 'debug lon_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', lon_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! (Gather latitude)
IF (myPE==0)  ALLOCATE(lat_tempG(dimobs_tempG))
CALL MPI_GATHERV(coordinate_temp(2,:), dimobs_temp,                        MPI_REAL,       &
                 lat_tempG,            dimobs_temp_part, dimobs_temp_dspl, MPI_REAL, 0,    & 
                 MPI_COMM_WORLD, ierror)
DEALLOCATE(coordinate_temp)
IF (myPE==0) WRITE(*,*) 'debug lat_tempG(1:4) ', lat_tempG(1:4)
IF (myPE==0) WRITE(*,*) 'debug lat_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4) ', lat_tempG(dimobs_temp_part(1)+1:dimobs_temp_part(1)+4)

! **********************************************************************
! *** Write global observations to NetCDF:
! **********************************************************************
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Write global observations to NetCDF:'
WRITE(*,*) '**********************************************************************'
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

IF (myPE==0) THEN ! write global obs. to netCDF (myPE 0)
s=1

filename = trim('/work/ollie/frbunsen/model_runs/fesom2/postprocessing/observations/global_observations.nc')

stat(s) = NF_create(trim(filename), 0, fileID)
s=s+1
stat(s) = NF_def_dim(fileID,'dim_temp',dimobs_tempG,dimid_temp)
s=s+1
stat(s) = NF_def_var(fileID,'temperature',NF_FLOAT,1,dimid_temp,varid_temperature)
s=s+1
stat(s) = NF_def_var(fileID,'day_temp',NF_INT,1,dimid_temp,varid_daytemp)
s=s+1
stat(s) = NF_def_var(fileID,'lon_temp',NF_FLOAT,1,dimid_temp,varid_lontemp)
s=s+1
stat(s) = NF_def_var(fileID,'lat_temp',NF_FLOAT,1,dimid_temp,varid_lattemp)
s=s+1
stat(s) = NF_def_var(fileID,'lev_temp',NF_FLOAT,1,dimid_temp,varid_levtemp)
s=s+1
stat(s) = NF_enddef(fileID)
s=s+1

stat(s) = NF_put_var_real(fileID,varid_temperature,temperatureG)
stat(s) = NF_put_var_int (fileID,varid_daytemp,day_tempG)
stat(s) = NF_put_var_real(fileID,varid_lontemp,lon_tempG)
stat(s) = NF_put_var_real(fileID,varid_lattemp,lat_tempG)
stat(s) = NF_put_var_real(fileID,varid_levtemp,nl1_tempG)

stat(s) = NF_close(fileID)
s=s+1

DO i = 1, s-1
   IF (stat(i) /= NF_NOERR) &
		WRITE(*,*) 'NetCDF error in writing observations to NetCDF, no.', i
END DO

END IF! write global obs. to netCDF (myPE 0)


! **********************************************************************
! *** Write global interpolated analysis and forecast to NetCDF:
! **********************************************************************
IF (myPE==0) THEN
WRITE(*,*) ''
WRITE(*,*) '**********************************************************************'
WRITE(*,*) '*** Write global interpolated analysis and forecast to NetCDF:'
WRITE(*,*) '**********************************************************************'
END IF
CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

IF (myPE==0) THEN ! write global ana/forecast to netCDF (myPE 0)
s=1

filename = trim('/work/ollie/frbunsen/model_runs/fesom2/postprocessing/')//trim(runname)//trim('/FESOM-PDAF-out-interpolated-to-EN4.nc')

WRITE(*,*) 'Writing netCDF to ', trim(filename)

stat(s) = NF_create(trim(filename), 0, fileID)
s=s+1
stat(s) = NF_def_dim(fileID,'dim_temp',dimobs_tempG,dimid_temp)
s=s+1
stat(s) = NF_def_var(fileID,'t_ana',NF_FLOAT,1,dimid_temp,varid_tempa)
s=s+1
stat(s) = NF_def_var(fileID,'t_for',NF_FLOAT,1,dimid_temp,varid_tempf)
s=s+1
stat(s) = NF_def_var(fileID,'day_temp',NF_INT,1,dimid_temp,varid_daytemp)
s=s+1
stat(s) = NF_def_var(fileID,'lon_temp',NF_FLOAT,1,dimid_temp,varid_lontemp)
s=s+1
stat(s) = NF_def_var(fileID,'lat_temp',NF_FLOAT,1,dimid_temp,varid_lattemp)
s=s+1
stat(s) = NF_def_var(fileID,'lev_temp',NF_FLOAT,1,dimid_temp,varid_levtemp)
s=s+1
stat(s) = NF_enddef(fileID)
s=s+1

stat(s) = NF_put_var_real(fileID,varid_tempa,t_ana_nG)
stat(s) = NF_put_var_real(fileID,varid_tempf,t_for_nG)
stat(s) = NF_put_var_int (fileID,varid_daytemp,day_tempG)
stat(s) = NF_put_var_real(fileID,varid_lontemp,lon_tempG)
stat(s) = NF_put_var_real(fileID,varid_lattemp,lat_tempG)
stat(s) = NF_put_var_real(fileID,varid_levtemp,nl1_tempG)

stat(s) = NF_close(fileID)
s=s+1

DO i = 1, s-1
   IF (stat(i) /= NF_NOERR) &
		WRITE(*,*) 'NetCDF error in writing analysis/forecast to NetCDF, no.', i
END DO

END IF! write global ana/forecast to netCDF (myPE 0)

! **********************************************************************
! *** Finalize:
! **********************************************************************
CALL MPI_FINALIZE(ierror)
END PROGRAM interpolateEN4

