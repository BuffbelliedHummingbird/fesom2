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

INTEGER :: i,j,k,n                      ! counters
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
INTEGER, ALLOCATABLE :: n2d_temp(:,:), nl1_temp(:), offset_temp(:)  ! nodes and layer of observation, offset in file
INTEGER, ALLOCATABLE :: n2d_salt(:,:), nl1_salt(:), offset_salt(:)
REAL(4), ALLOCATABLE :: coordinate_temp(:,:), temperature(:)        ! geo. coordinates [-Pi:Pi]
REAL(4), ALLOCATABLE :: coordinate_salt(:,:), salinity(:)
INTEGER :: varid_n2d_temp, varid_nl1_temp, varid_offset_temp
INTEGER :: varid_n2d_salt, varid_nl1_salt, varid_offset_salt
INTEGER :: varid_coordinate_temp, varid_temperature
INTEGER :: varid_coordinate_salt, varid_salinity
INTEGER :: dimobs_temp, dimobs_salt                                 ! number of observations (in one year)
INTEGER :: DimId_n_obs_temp, DimID_n_obs_sal
INTEGER :: day, ndays                                               ! day counter, number of days in one year


! Parallel:
INTEGER        :: myPE
CHARACTER(4)   :: myPEstr4               ! myPE as str
CHARACTER(5)   :: myPEstr5               ! myPE as str
INTEGER        :: ierror, sizecomm, &
                  mpistatus(MPI_STATUS_SIZE), &
                  tag
                  
REAL(4), ALLOCATABLE :: TEST(:,:,:)

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
!!!!!!!!!!! SHOULD BE npes-1 IN THE NEXT LINE BUT THIS GIVES MEMORY ERROR. !!!!!!!!!!!!
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
END DO
DEALLOCATE(t_anaG)
DEALLOCATE(t_forG)
DEALLOCATE(s_anaG)
DEALLOCATE(s_forG)
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

END IF
IF (myPE==1) WRITE(*,*) 'debug t_ana(200,:,1): ', t_ana(200,:,1)

CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)



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
	
	ALLOCATE(n2d_temp(3,dimobs_temp))
	ALLOCATE(offset_temp(ndays))
	ALLOCATE(coordinate_temp(2,dimobs_temp))
	ALLOCATE(temperature(dimobs_temp))
	ALLOCATE(nl1_temp(dimobs_temp))
	         
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
! *** Finalize:
! **********************************************************************
CALL MPI_FINALIZE(ierror)
END PROGRAM interpolateEN4

