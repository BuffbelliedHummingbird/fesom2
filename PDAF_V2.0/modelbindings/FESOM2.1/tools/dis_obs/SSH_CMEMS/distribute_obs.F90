! $Id: distribute_obs.F90 644 2012-06-20 09:27:59Z lnerger $
!
! Program: distribute_obs --- Generate distributed observation files
!
PROGRAM distribute_obs

! DESCRIPTION:
! This program generates files holding distributed observation
! information according to the distribution information of the
! mesh. The observation infromation is read from a NetCDF file 
! holding the global observation information.

! USES:
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

! Local variables
  INTEGER :: i, s, iter, pe
  LOGICAL :: twin_data
  CHARACTER(len=120) :: inpath, outpath, distpath
  CHARACTER(len=120) :: infile, outfile, partfile, mylistfile
  CHARACTER(len=150) :: ncfile_in, ncfile_out
  CHARACTER(len=150) :: attstr
  CHARACTER(len=120) :: file_name                                       ! myListXXXXX numbered filenames
  CHARACTER(len=10)  :: mype_string
  INTEGER :: fileID                                                     ! ID for myListXXXXX files
  INTEGER :: ncid_in, ncid_out
  INTEGER :: id_dim, id_stderr, id_time, id_list
  INTEGER :: id_obs, id_my_obs, id_std, id_my_std
  INTEGER :: dimid_iter, dimid_one, dimid_n2d
  INTEGER :: varid_time
  INTEGER :: n2d, steps
  INTEGER :: idummy, n
  INTEGER :: sum_n2d
  INTEGER :: stat(100)
  INTEGER :: countv(2), startv(2)
  INTEGER :: dimids(2)

  REAL, ALLOCATABLE :: times(:)
  REAL, ALLOCATABLE :: obs(:), my_obs(:), std(:), my_std(:)
  REAL :: stderr_obs
  INTEGER, ALLOCATABLE :: cntobs(:), cntobsday(:)
  
  INTEGER :: npes                                                       ! number of model processes
  INTEGER :: myDim_nod2D, myDim_elem2D                                  ! number of pe-local nodes / elements
  INTEGER :: eDim_nod2D, eDim_elem2D                                    ! number of pe-local halo nodes / elements
  INTEGER :: eXDim_elem2D
  INTEGER, ALLOCATABLE :: nod2D_part(:), elem2D_part(:)                 ! rank partioning vectors (list of nodes, elems per process)
  INTEGER, ALLOCATABLE :: List_nod2D(:), List_elem2D(:)                 ! (dummy) arrays to read from myList files
  INTEGER, ALLOCATABLE :: myList_nod2D(:), myList_elem2D(:)             ! (dummy) arrays to read from myList files
  
  INTEGER :: year
  CHARACTER(len=4)  :: year_string

! ************************************************
! *** Configuration                            ***
! ************************************************

  year = 2016
  write(year_string,'(i4.4)') year

  ! Path to and name of file holding global observations on model grid
  ! inpath = '/work/ollie/frbunsen/data/physics/SSH/CMEMS/CORE2/'
  inpath = '/albedo/work/projects/p_recompdaf/frbunsen/data/physics/SSH/CMEMS/CORE2/'
  infile = 'SSH_track_'//year_string//'.nc'
  
  ! Path to and name stub of output files
  ! outpath = '/work/ollie/frbunsen/data/physics/SSH/CMEMS/CORE2/dist_144/'
  outpath = '/albedo/work/projects/p_recompdaf/frbunsen/data/physics/SSH/CMEMS/CORE2/dist_128/'
  outfile = 'SSH_track_'//year_string//'_dist128'

  ! Path to mesh partioning
  ! distpath   = '/work/ollie/projects/clidyn/FESOM2/meshes/core2/dist_128/'
  distpath = '/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/dist_128/'
  partfile   = 'rpart.out'
  mylistfile = 'my_list'

  ! Whether the observation file hold twin data 
  ! only in this case the estimated stderr_obs and time array are handled
  twin_data = .FALSE.


! ************************************************
! *** Init                                     ***
! ************************************************

  WRITE (*,'(3x,a)')  '********************************************************'
  WRITE (*,'(3x,a)')  '*** Generate distributed observation files for FESOM ***'
  WRITE (*,'(3x,a/)') '********************************************************'

  ncfile_in = TRIM(inpath)//TRIM(infile)
  WRITE (*,*) 'Read fields from file: ',TRIM(ncfile_in)

  ncfile_out = TRIM(outpath)//TRIM(outfile)
  WRITE (*,*) 'Write distributed fields to files: ',TRIM(ncfile_out),'_XXXX.nc'

  WRITE (*,*) 'Read partitioning information from directory: ',TRIM(distpath)


! *******************************************
! *** Open obs file and read dimensions   ***
! *******************************************

  s = 1
  stat(s) = NF_OPEN(TRIM(ncfile_in), NF_NOWRITE, ncid_in)
  s = s + 1

  ! Get dimensions
  stat(s) = NF_INQ_DIMID(ncid_in, 'nodes_2D', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, n2d)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'time', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, steps)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
  END DO

  WRITE (*,'(/1x,a)') 'Global dimensions of experiment:'
  WRITE (*,'(5x,1x,a3,i10)') 'n2d', n2d
  WRITE (*,'(5x,a4,i10)') 'iter', steps


! ***********************************
! *** Read mesh partitioning data ***
! ***********************************

  OPEN(unit=1,file=TRIM(distpath)//TRIM(partfile), status='old', form='formatted')

  READ(1,*) npes
  CLOSE(1)

  ALLOCATE( nod2D_part(npes))
  ALLOCATE(elem2D_part(npes))
    
! *******************************************************
! *** Lists of nodes and elements in global indexing. ***
! *** every proc reads its file                       ***
! *******************************************************
  
  DO i = 1, npes
  
	write(mype_string,'(i5.5)') i-1
	file_name = trim(distpath)//'my_list'//trim(mype_string)//'.out'  
	fileID = i+1
		
	open(fileID, file=trim(file_name))
	read(fileID,*) n                      ! myPE
 
	read(fileID,*) myDim_nod2D            ! number of vertices that belong to myPE
	read(fileID,*) eDim_nod2D             ! number of halo vertices (share a common
	                                      ! triangle with vertices that belong to myPE,
	                                      ! yet do not belong to myPE themselves.
	allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
	read(fileID,*) myList_nod2D           ! list of vertices
	
	nod2D_part(i) = myDim_nod2D

	read(fileID,*) myDim_elem2D           ! triangles that belong to myPE
	read(fileID,*) eDim_elem2D            ! triangles sharing an edge with my triangles
	read(fileID,*) eXDim_elem2D           ! triangles sharing a vertix with my triangles
	allocate(myList_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
	read(fileID,*) myList_elem2D
	
	elem2D_part(i) = myDim_elem2D

	close(fileID)
	
	deallocate(myList_nod2D)
	deallocate(myList_elem2D)
  
  END DO

  ! *** Consistency check ***
  IF (SUM(nod2D_part) /= n2D) THEN
	WRITE (*,*) 'Partitioned mesh not consistent with global mesh (nodes)'
	STOP
  END IF
  
  
  ! *** Write ***
  WRITE (*,'(/1x,a)') 'Mesh distribution:'
  WRITE (*,'(5x,a,i6)') 'Number of PEs:', npes
  WRITE (*,'(5x,a)') 'Local nodes:   nod2d'
  DO i = 1, npes
     WRITE (*,'(18x,i7)') nod2d_part(i)
  END DO


! ************************************
! *** Initialize distributed files ***
! ************************************

  ! Read stderr and time information from input file

  ALLOCATE(times(steps))

  IF (twin_data) THEN
    s = 1
    stat(s) = NF_INQ_VARID(ncid_in, 'time', id_time)
    s = s + 1
    stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_time, times)
    s = s + 1
    stat(s) = NF_INQ_VARID(ncid_in, 'std', id_stderr)
    s = s + 1
    stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_stderr, stderr_obs)

    DO i = 1,  s
       IF (stat(i) /= NF_NOERR) &
            WRITE(*, *) 'NetCDF error in reading from input file, no.', i
    END DO
  ENDIF

  ! Initialize one file for each PE
  initoutfiles: DO pe = 0, npes - 1

     WRITE(mype_string,'(i5.5)') pe

     s = 1
     stat(s) = NF_CREATE(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', 0, ncid_out) 

     attstr  = 'daily observations from year '//year_string
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(attstr), &
          TRIM(attstr)) 

     ! Define dimensions
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'iter', steps, dimid_iter)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'one', 1, dimid_one)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'nodes_2D', nod2d_part(pe + 1), dimid_n2d)

     ! Define variables
     IF (twin_data) THEN
        s = s + 1
        stat(s) = NF_DEF_VAR(ncid_out, 'stderr_obs', NF_DOUBLE, 1, dimid_one, Id_stderr) 
        s = s + 1
        stat(s) = NF_DEF_VAR(ncid_out, 'time', NF_DOUBLE, 1, dimid_iter, Id_time) 
     ENDIF
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'nodlist', NF_INT, 1, dimid_n2d, id_list)

     dimids(1) = DimId_n2d
     dimids(2) = dimid_iter

     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'obs', NF_DOUBLE, 2, dimids(1:2), id_my_obs)
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'std', NF_DOUBLE, 2, dimids(1:2), id_my_std)
     s = s + 1
     stat(s) = NF_ENDDEF(ncid_out) 

     IF (twin_data) THEN
        ! Write std error
        s = s + 1
        stat(s) = NF_PUT_VAR_DOUBLE(ncid_out, Id_stderr, stderr_obs)

        ! Write times
        s = s + 1
        stat(s) = NF_PUT_VAR_DOUBLE(ncid_out, Id_time, times)
     ENDIF

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in init of output file, no.', i
     END DO

  END DO initoutfiles


! ***********************************
! *** Generate distributed fields ***
! ***********************************

  WRITE (*,'(/1x,a/)') '------- Generate distributed files -------------'

  ! Allocate global obs field
  ALLOCATE(obs(n2d))
  ALLOCATE(std(n2d))
  ALLOCATE(cntobs(steps))
  ALLOCATE(cntobsday(npes))
  cntobs = 0
  cntobsday = 0

  peloop: DO pe = 0, npes - 1

! Read my_listXXXXX files
    WRITE (*,'(1x,a)') '*** Read my_listXXXXX files ***'
    
	write(mype_string,'(i5.5)') pe
	file_name = trim(distpath)//'my_list'//trim(mype_string)//'.out'  
	fileID = pe+1
		
	open(fileID, file=trim(file_name))
	read(fileID,*) n
 
	read(fileID,*) myDim_nod2D
	read(fileID,*) eDim_nod2D
	allocate(List_nod2D(myDim_nod2D+eDim_nod2D))
    read(fileID,*) List_nod2D
	
	allocate(myList_nod2D(myDim_nod2D))
	myList_nod2D = List_nod2D(1:myDim_nod2D)
	
	read(fileID,*) myDim_elem2D
	read(fileID,*) eDim_elem2D
	read(fileID,*) eXDim_elem2D
	allocate(List_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
	read(fileID,*) List_elem2D

	allocate(myList_elem2D(myDim_elem2D))
	myList_elem2D = List_elem2D(1:myDim_elem2D)
	
	close(fileID)
	
	deallocate(List_nod2D)
	deallocate(List_elem2D)
	
    ALLOCATE(my_obs     (myDim_nod2D))
    ALLOCATE(my_std     (myDim_nod2D))

     ! Open output file
     s = 1
     stat(s) = NF_OPEN(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', NF_WRITE, ncid_out)

     ! Write MyList
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, myList_nod2D)

     ! Get ID for global obs field
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'ssh', id_obs)

     ! Get ID for local obs field
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'obs', id_my_obs)


     ! Get ID for global STD
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'std', id_std)

     ! Get ID for local STD
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'std', id_my_std)

     ! Loop through time slices and generate and write local fields
     stepsloop: DO iter = 1, steps

        ! Read global obs field
        startv(2) = iter
        countv(2) = 1
        startv(1) = 1
        countv(1) = n2d  
        stat(1) = NF_GET_VARA_DOUBLE(ncid_in, id_obs, startv, countv, obs)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading obs'

        ! Read global std field
        startv(2) = iter
        countv(2) = 1
        startv(1) = 1
        countv(1) = n2d  
        stat(1) = NF_GET_VARA_DOUBLE(ncid_in, id_std, startv, countv, std)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading STD'

        ! Select local nodes
        DO i = 1, nod2d_part(pe+1)
           my_obs(i) = obs(myList_nod2D(i))
           my_std(i) = std(myList_nod2D(i))
           IF (my_obs(i)>-100.0 .AND. my_obs(i)<100.0) cntobs(iter) = cntobs(iter)+1
           IF (iter==1) THEN
              IF (my_obs(i)>-100.0 .AND. my_obs(i)<100.0) cntobsday(pe+1) = cntobsday(pe+1)+1
           END IF
        END DO

        ! Write local obs field
        startv(2) = iter
        countv(2) = 1
        startv(1) = 1
        countv(1) = nod2d_part(pe+1)
        
        stat(1) = NF_PUT_VARA_DOUBLE(ncid_out, id_my_obs, startv, countv, my_obs)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in writing local obs'

        ! Write local std field
        startv(2) = iter
        countv(2) = 1
        startv(1) = 1
        countv(1) = nod2d_part(pe+1)
        stat(1) = NF_PUT_VARA_DOUBLE(ncid_out, id_my_std, startv, countv, my_std)
        IF (stat(1) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in writing local STD'

     END DO stepsloop

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in writing to output file, no.', i
     END DO

     DEALLOCATE(myList_nod2D, myList_elem2D, my_obs, my_std)

  END DO peloop

  stat(1) = nf_close(ncid_in)

  DO iter=1,steps
     WRITE (*,*) 'step, cntobs', iter, cntobs(iter)
  END DO
  DO pe=1,npes
     write (*,*) 'pe, cntobsday1', pe-1, cntobsday(pe)
  END DO

  WRITE (*,'(1x,a/)') '------- END -------------'

  DEALLOCATE(cntobs, cntobsday, obs, std)

END PROGRAM distribute_obs
