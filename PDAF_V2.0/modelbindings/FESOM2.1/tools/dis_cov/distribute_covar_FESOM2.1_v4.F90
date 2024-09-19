! PDAF_V2.0

! $Id: distribute_covar.F90 1604 2016-05-30 06:42:16Z lnerger $

! Modified by QT     2018-02-27 ---- for AWI-CM
! Modified by Frauke 2022-02-28 ---- for FESOM2.1 with core2 mesh

! !Program: distribute_covar --- Generate distributed covariance matrix files
!
! !INTERFACE:
PROGRAM distribute_covar

! !DESCRIPTION:
! This program generates files holding distributed covariance matrix
! information according to the distribution information of the
! mesh. The global covariance matrix is read from a NetCDF file.
! The input file contains values in single precision. For compatibility
! with FESOM the distributed output files will be in double precision.
!
! Version 2: Velocities are interpolated from elements on nodes.

  USE mpi
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'


  ! Local variables
  INTEGER :: mpierror,errorcode                                         ! Error status for MPI
  INTEGER :: i, s, iter, pe, irank, n, b                                ! counters
  INTEGER :: j, k                                                       ! counters
  CHARACTER(len=120) :: inpath, outpath, dist_mesh_dir
  CHARACTER(len=120) :: infile, outfile, partfile, mylistfile
  CHARACTER(len=150) :: ncfile_in, ncfile_out
  CHARACTER(len=10)  :: mype_string                                     ! for numbered filenames for each process
  CHARACTER(len=120) :: file_name                                       ! myListXXXXX numbered filenames
  CHARACTER(len=150) :: title                                           ! reading global cov.nc
  CHARACTER(len=150) :: fieldsstr                                       ! reading global cov.nc
  INTEGER :: fileID                                                     ! ID for myListXXXXX files
  INTEGER :: ncid_in, ncid_out
  INTEGER :: id_dim
  INTEGER :: id_sigma, id_list
  INTEGER :: id_mean, id_svec
  INTEGER :: dimid_rank, dimid_one, dimid_nfields, dimid_state          ! NC IDs for pe-distributed outfiles
  INTEGER :: dimid_2D, dimid_tracer3D                                   ! NC IDs for pe-distributed outfiles
  INTEGER :: id_out
  INTEGER :: startpos
  INTEGER :: rank, dim_state, nfields
  INTEGER :: nlevels                                                    ! Number of vertical levels
  INTEGER :: nod2D, elem2D
  INTEGER :: stat(100)
  INTEGER :: dim_state_l                                                ! length of pe-local state vector
  INTEGER :: countv(2), startv(2)
  INTEGER :: dimids(2)
  INTEGER :: npes                                                       ! number of model processes
  CHARACTER(len=2) :: npes_char                                         ! number of model processes in file name
  INTEGER :: myDim_nod2D, myDim_elem2D                                  ! number of pe-local nodes / elements
  INTEGER :: eDim_nod2D, eDim_elem2D                                    ! number of pe-local halo nodes / elements
  INTEGER :: eXDim_elem2D
  INTEGER, ALLOCATABLE :: nod2D_part(:), elem2D_part(:)                 ! rank partioning vectors (list of nodes, elems per process)
  REAL(kind=4), ALLOCATABLE :: field_surf(:), field_3D(:)               ! global field vectors
  REAL(kind=4), ALLOCATABLE :: svals(:)
  REAL(kind=4), ALLOCATABLE :: myfield_surf(:), myfield_3D(:)           ! pe-local field vectors
  REAL(kind=8), ALLOCATABLE :: state_l(:)                               ! pe-local state vector
  INTEGER, ALLOCATABLE :: List_nod2D(:), List_elem2D(:)                 ! (dummy) arrays to read from myList files
  INTEGER, ALLOCATABLE :: myList_nod2D(:), myList_elem2D(:)             ! (dummy) arrays to read from myList files
  INTEGER, ALLOCATABLE :: my_nodlist_surf(:), my_nodlist_3D(:)          ! pe-local node lists (level-wise repeated lists for 3D variables)
  
  INTEGER :: dim_fields3D, dim_fieldssurf
  INTEGER, ALLOCATABLE :: dim_fields(:)                                              ! Field dimensions
  INTEGER, ALLOCATABLE :: offsets(:)                                                 ! Field offsets in state vector
  INTEGER, ALLOCATABLE :: dim_fields_l(:)                                            ! dimensions of pe-local fields  
  
integer :: fmin, fmax                     ! Number of fields
  
  TYPE field_ids                          ! Type declaration, holding the indices of model fields in the state vector
     INTEGER :: ssh
     INTEGER :: u
     INTEGER :: v
     INTEGER :: temp
     INTEGER :: salt
     INTEGER :: DIC
     INTEGER :: Alk
     INTEGER :: DIN
     INTEGER :: O2
END TYPE field_ids
  
  TYPE(field_ids) :: idx                  ! Type variable holding field IDs in state vector
                          
  type state_field
   integer :: ndims = 0                   ! Number of field dimensions (1 or 2)
   character(len= 10) :: variable = ''    ! Name of field
   character(len=250) :: filename = ''    ! File name input incl. path
   integer :: ncid = 0                    ! ID nc-file input
   integer :: fid = 0                     ! Field ID (in/output file)
   integer :: mid = 0                     ! Mean field ID (output file)
   integer :: sid = 0                     ! Singular vector ID
  end type state_field

  type(state_field), allocatable :: fproperties(:) ! Type variable holding the
                                                   ! definitions of model fields

! ************************************************
! *** Configuration                            ***
! ************************************************

  CALL MPI_INIT(mpierror)
  
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, npes, MPIERROR)
  write(npes_char,'(i2.2)') npes
  
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, pe,   MPIERROR)
  write(mype_string,'(i5.5)') pe
  
  
  ! Path to and name of file holding global covariance matrix
  inpath = '/albedo/work/projects/p_recompdaf/frbunsen/modelruns/BGC9_covar/'
  infile = 'cov.nc'

  ! Path to mesh
  dist_mesh_dir = '/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/dist_'//trim(npes_char)//'/'
  partfile      = 'rpart.out'
  mylistfile    = 'my_list'

  ! Path to and name stub of output files
  outpath = '/albedo/work/projects/p_recompdaf/frbunsen/modelruns/BGC9_covar/dist'//trim(npes_char)//'/'
  outfile = 'cov'
  
  ! Composition of state vector:
  nfields = 9
  
  fmin = 1
  fmax = nfields
  
  idx% ssh  = 1
  idx% u    = 2
  idx% v    = 3
  idx% temp = 4
  idx% salt = 5
  idx% DIC  = 6
  idx% Alk  = 7
  idx% DIN  = 8
  idx% O2   = 9
  
  ! Field-specific variables:
  allocate(fproperties(nfields))
  
  fproperties(idx% ssh) % ndims = 1
  fproperties(idx% ssh) % variable = 'ssh'
  
  fproperties(idx% u) % ndims = 2
  fproperties(idx% u) % variable = 'u'

  fproperties(idx% v) % ndims = 2
  fproperties(idx% v) % variable = 'v'

  fproperties(idx% temp) % ndims = 2
  fproperties(idx% temp) % variable = 'temp'

  fproperties(idx% salt) % ndims = 2
  fproperties(idx% salt) % variable = 'salt'

  fproperties(idx% DIC) % ndims = 2
  fproperties(idx% DIC) % variable = 'DIC' 

  fproperties(idx% Alk) % ndims = 2
  fproperties(idx% Alk) % variable = 'Alk'

  fproperties(idx% DIN) % ndims = 2
  fproperties(idx% DIN) % variable = 'DIN'

  fproperties(idx% O2) % ndims = 2
  fproperties(idx% O2) % variable = 'O2'
  

! ************************************************
! *** Init                                     ***
! ************************************************

  if (pe==0) WRITE (*,'(3x,a)')  '*****************************************************************'
  if (pe==0) WRITE (*,'(3x,a)')  '*** Generate distributed covariance matrix files for FESOM2.1 ***'
  if (pe==0) WRITE (*,'(3x,a/)') '*****************************************************************'

  ncfile_in = TRIM(inpath)//TRIM(infile)
  if (pe==0) WRITE (*,*) 'Read fields from file: ',TRIM(ncfile_in)

  ncfile_out = TRIM(outpath)//TRIM(outfile)
  if (pe==0) WRITE (*,*) 'Write distributed fields to files: ',TRIM(ncfile_out),'_XXXX.nc'
  if (pe==0) WRITE (*,*) 'Read partitioning information from directory: ',TRIM(dist_mesh_dir)


! ********************************************
! *** Open global file and read dimensions ***
! ********************************************

  if (pe==0) then ! netCDF pe
  s = 1
  stat(s) = NF_OPEN(TRIM(ncfile_in), NF_NOWRITE, ncid_in)

  ! Get number of fields
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nfields', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, nfields)
  endif ! netCDF pe
  
  if (pe==0) WRITE(*,*) 'Broadcasting nfields'
  CALL MPI_BCAST(nfields,1,MPI_INT,0,MPI_COMM_WORLD)
  
  allocate(dim_fields(nfields))
  allocate(offsets(nfields))
  allocate(dim_fields_l(nfields))
  
  if (pe==0) then ! netCDF pe
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, rank)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nod2', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_fieldssurf) ! for surface fields, e. SSH
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nz_x_nod2', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_fields3D) ! for 3D fields, eg. temp, salt, u, v
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_state)
  
  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
  END DO
  
  endif ! netCDF pe
  
  if (pe==0) WRITE(*,*) 'Broadcasting SVD rank'
  CALL MPI_BCAST(rank,        1,MPI_INT,0,MPI_COMM_WORLD, MPIERROR)
  if (pe==0) WRITE(*,*) 'Broadcasting dim_fields'
  CALL MPI_BCAST(dim_fieldssurf,1,MPI_INT,0,MPI_COMM_WORLD, MPIERROR)
  CALL MPI_BCAST(dim_fields3D,1,MPI_INT,0,MPI_COMM_WORLD, MPIERROR)
    
  ! dimensions of each field
  do b=fmin, fmax
    if (fproperties(b)% ndims == 1) dim_fields(b) = dim_fieldssurf ! surface fields
    if (fproperties(b)% ndims == 2) dim_fields(b) = dim_fields3D   ! 3D fields
  enddo
  
  nlevels = 46         ! number of vertical levels
  nod2D   = 126858     ! number of nodes
  elem2D  = 244659     ! number of elements (not used)
  
  if (pe==0) then ! write pe
  ! Write dimensions
  WRITE (*,'(/1x,a)') 'Global dimensions of experiment:'
  do b=fmin, fmax
    WRITE (*,'(10x,1x,a30,i12)') 'dim_fields '//trim(fproperties(b)% variable), dim_fields(b)
  enddo
  
  WRITE (*,'(10x,1x,a30,i12)') 'nlevels                       ', nlevels
  WRITE (*,'(10x,1x,a30,i12)') 'rank                          ', rank
  WRITE (*,'(10x,1x,a30,i12)') 'dim_state                     ', dim_state
  WRITE (*,'(10x,1x,a30,i12)') 'nfields                       ', nfields
  
  ! Bug check
  IF (dim_fields(idx% ssh)*nlevels /= dim_fields(idx% temp)) THEN
	WRITE (*,*) 'Number of vertical levels and horizontal nodes not consistent with 3D nodes!'
	STOP
  END IF
  endif ! write pe

! ***********************************
! *** Read mesh partitioning data ***
! ***********************************
    
! *******************************************************
! *** Lists of nodes and elements in global indexing. ***
! *** every proc reads its file                       ***
! *******************************************************
  
	write(mype_string,'(i5.5)') pe
	file_name = trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
	fileID = pe+2
		
	open(fileID, file=trim(file_name))
	read(fileID,*) n                      ! myPE
 
	read(fileID,*) myDim_nod2D            ! number of vertices that belong to myPE
	read(fileID,*) eDim_nod2D             ! number of halo vertices (share a common
	                                      ! triangle with vertices that belong to myPE,
	                                      ! yet do not belong to myPE themselves.
	allocate(myList_nod2D(myDim_nod2D+eDim_nod2D)) 	 
	read(fileID,*) myList_nod2D           ! list of vertices
	
	read(fileID,*) myDim_elem2D           ! triangles that belong to myPE
	read(fileID,*) eDim_elem2D            ! triangles sharing an edge with my triangles
	read(fileID,*) eXDim_elem2D           ! triangles sharing a vertix with my triangles
	allocate(myList_elem2D(myDim_elem2D+eDim_elem2D+eXDim_elem2D))
	read(fileID,*) myList_elem2D
	
	close(fileID)
	
        if (pe==0) WRITE (*,'(5x,a,4x,a)') 'Local nodes:   nod2d','elem2d'
        CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERROR)
        WRITE(*,'(18x,i5,i7,i10)') pe, myDim_nod2D, myDim_elem2D
	
	ALLOCATE( nod2D_part(npes))
    ALLOCATE(elem2D_part(npes))
	
	CALL MPI_ALLGATHER(myDim_nod2D,  1, MPI_INT, nod2D_part,  1, MPI_INT, MPI_COMM_WORLD, MPIERROR)
	CALL MPI_ALLGATHER(myDim_elem2D, 1, MPI_INT, elem2D_part, 1, MPI_INT, MPI_COMM_WORLD, MPIERROR)

  ! elem2D_part(i) is the number of triangles that belong to each PE,
  ! which are those that contain at least one vertex that belongs to the PE.
  ! Thus, triangles with vertices that belong to several PE are counted
  ! as my_triangle by each PE.
  
  if (pe==0) then ! write pe
  
    WRITE(*,*) 'Sum of elements in partitioned mesh: ', SUM(elem2D_part)
    WRITE(*,*) 'Sum of elements in global mesh:      ', elem2D
    
    WRITE (*,'(/1x,a)') 'Mesh distribution:'
    WRITE (*,'(5x,a,i6)') 'Number of PEs:', npes
    WRITE (*,*) 'Local nodes gathered:'
    DO i = 1, npes
       WRITE (*,*) 'mype:', i, 'nod2D_part:', nod2D_part(i), 'elem2D_part', elem2D_part(i)
    END DO
    
    ! *** Consistency check ***
    IF (SUM(nod2D_part) /= nod2D) THEN
          WRITE (*,*) 'Partitioned mesh not consistent with global mesh (nodes)'
          CALL MPI_ABORT(MPI_COMM_WORLD,errorcode,MPIERROR)
    END IF

    WRITE (*,'(/1x,a/)') '------- Initialize distributed files -------------'

  endif ! write pe

! ************************************
! *** Initialize distributed files ***
! ************************************

  ! Read singular values from input file

  ALLOCATE(svals(rank))

  if (pe==0) then ! netCDF pe
  s = 1
  stat(s) = NF_INQ_VARID(ncid_in, 'sigma', id_sigma)
  s = s + 1
  stat(s) = NF_GET_VAR_REAL(ncid_in, id_sigma, svals)
  s = s + 1
  title = ''
  stat(s) = NF_GET_ATT_TEXT(ncid_in, NF_GLOBAL, 'title', title)
  s = s + 1
  fieldsstr = ''
  stat(s) = NF_GET_ATT_TEXT(ncid_in, NF_GLOBAL, 'state_fields', fieldsstr)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading from input file, no.', i
  END DO
  endif ! netCDF pe
  
  CALL MPI_BCAST(svals,     rank,       MPI_DOUBLE,    0, MPI_COMM_WORLD, MPIERROR)
  CALL MPI_BCAST(title,     len(title), MPI_CHARACTER, 0, MPI_COMM_WORLD, MPIERROR)
  CALL MPI_BCAST(fieldsstr, len(title), MPI_CHARACTER, 0, MPI_COMM_WORLD, MPIERROR)

     ! PE-local dimensions of fields in state vector
     do b=fmin,fmax
       if (fproperties(b)% ndims == 1) dim_fields_l(b) = nod2D_part(pe+1)
       if (fproperties(b)% ndims == 2) dim_fields_l(b) = nod2D_part(pe+1) * nlevels
     enddo

     ! PE-local dimension of state vector
     dim_state_l = SUM(dim_fields_l)
     
     write(mype_string,'(i4.4)') pe
     WRITE (*,'(a, i4, 4x, i7)') 'pe and dim_state_l: ', pe, dim_state_l

     ! Initialize one file on each PE
     s = 1
     stat(s) = NF_CREATE(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', NF_64BIT_OFFSET, ncid_out)
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(title), &
          TRIM(title)) 
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'state_fields', LEN_TRIM(fieldsstr), &
          TRIM(fieldsstr)) 

	  ! Define dimensions
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'rank',  rank, dimid_rank)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'nod2',  nod2D_part(pe+1), dimid_2D)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'nz_x_nod2', nod2D_part(pe+1)*nlevels, dimid_tracer3D)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'dim_state', dim_state_l, dimid_state)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'one',  1, dimid_one)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'nfields',  nfields, dimid_nfields)

     ! Define variables

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'sigma', NF_DOUBLE, 1, dimid_rank, id_sigma)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_surf', NF_INT, 1, dimid_2D, id_list)              ! surface nodes
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_tracer3D',  NF_INT, 1, dimid_tracer3D,  id_list)  ! nodes repeated on levels

     ! running mean state (for the last snap shot)
     dimids(1) = dimid_state
     dimids(2) = dimid_one

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'running_meanstate', NF_DOUBLE, 2, dimids, id_mean)

     ! singular state vectors
     dimids(1) = dimId_state
     dimids(2) = dimid_rank

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'V', NF_DOUBLE, 2, dimids, Id_svec)

     ! End define mode
     s = s + 1
     stat(s)=NF_ENDDEF(ncid_out)

     ! Write singular values
     s = s + 1
     stat(s) = NF_PUT_VAR_DOUBLE(ncid_out, id_sigma, REAL(svals(1:rank),8))

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in init of output file, no.', i
     END DO


! ****************************************************
! *** Generate and write distributed state vectors ***
! ****************************************************

  if (pe==0) WRITE (*,'(1x,a/)') '------- Generate and write distributed files -------------'

    ! Allocate global fields
     ALLOCATE(field_surf (dim_fields(idx% ssh)))
     ALLOCATE(field_3D   (dim_fields(idx% temp)))

     ! Allocate local fields
	 ALLOCATE(myfield_surf (dim_fields_l(idx% ssh)))
	 ALLOCATE(myfield_3D   (dim_fields_l(idx% temp)))
	 
     ALLOCATE(state_l(dim_state_l))

     ! Define offsets
     offsets (1)          = 0
     do b=2,fmax
       offsets(b) = offsets(b-1) + dim_fields_l(b-1)
     enddo
     
     if (pe==0) THEN
             WRITE (*,'(1x,a)') '*** Create list of 3D nodes / elements ***'
             do b=fmin,fmax
                WRITE (*,'(1X,A30,1X,A5,1X,I7,1X,I7)') 'PE0-dim and offset of field ', &
                                                        trim(fproperties(b)%variable), &
                                                        dim_fields_l(b), &
                                                        offsets(b)
             enddo
     endif
	
	allocate(my_nodlist_surf (dim_fields_l(idx% ssh)))
	allocate(my_nodlist_3D   (dim_fields_l(idx% temp)))
	
	my_nodlist_surf = myList_nod2D

	! indeces of pe-local nodes / elements in global variable vector for 3D-fields
	DO k = 1, nlevels
      DO j = 1, myDim_nod2D
         my_nodlist_3D((j-1)*nlevels+k) = (myList_nod2D(j)-1) * nlevels + k
      END DO
    END DO
	
	deallocate(myList_nod2D)
	deallocate(myList_elem2D)

     ! Open output file
     write(mype_string,'(i4.4)') pe
     
     s = 1
     stat(s) = NF_OPEN(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', NF_WRITE, ncid_out)

     ! Write MyLists
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_surf', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, my_nodlist_surf)

     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_tracer3D', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, my_nodlist_3D)
     
     ! Inquire IDs for mean state and singular vectors
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'running_meanstate', id_mean)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'V', id_svec)

     if (pe==0) then ! netCDF pe
     
     ! Get IDs for global fields from input file
     do b=fmin,fmax
       s = s + 1
       stat(s) = NF_INQ_VARID(ncid_in, trim(fproperties(b)%variable)//'_mean', fproperties(b)% mid)
       s = s + 1
       stat(s) = NF_INQ_VARID(ncid_in, trim(fproperties(b)%variable)//'_svd', fproperties(b)% sid)
     enddo

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in inquire from output file, no.', i, ' pe ',pe
     END DO
     endif ! netCDF pe

     loopfields: DO irank = 0, rank

        IF (irank == 0) THEN
           ! Treat mean state

           if (pe==0) then ! netCDF pe
           ! Mean field Ids for input
           do b=fmin,fmax
             fproperties(b)% fid = fproperties(b)% mid
           enddo
           endif ! netCDF pe

           ! Id for output
           id_out = id_mean
           
           ! Column of matrix in file
           startpos = 1
           
        ELSE
           ! Treat singular vectors

           if (pe==0) then ! netCDF pe
           ! Singular vector Ids for input
           do b=fmin,fmax
             fproperties(b)% fid = fproperties(b)% sid
           enddo
           endif ! netCDF pe
           
           ! Id for output
           id_out = id_svec
           
           ! Column of matrix in file
           startpos = irank
           
        END IF

        if (pe==0) WRITE (*,*) '*** Generate local state vector (rank ',irank,') ***'
        ! *** Generate local state vector
        ! Read global state vector on PE 0 and broadcast
        s=1
        
        startv(2) = startpos
        countv(2) = 1
        startv(1) = 1
        
        do b=fmin,fmax
        
          countv(1) = dim_fields(b)
          
          if (fproperties(b)% ndims==1) then
            ! surface fields
            if (pe==0) then ! netCDF pe
              ! read global field
              s = s + 1
              stat(s) = NF_GET_VARA_REAL(ncid_in, fproperties(b)% fid, startv, countv, field_surf)
            endif ! netCDF pe
            CALL MPI_BCAST(field_surf, dim_fieldssurf, mpi_float, 0, mpi_comm_world, mpierror)
            ! fill local state vector:
             do i = 1, dim_fields_l(b)
              state_l(i + offsets(b)) = REAL(field_surf(my_nodlist_surf(i)),8)
             enddo
             
          elseif (fproperties(b)% ndims==2) then
            ! 3D-fields
            ! read global field
            if (pe==0) then ! netCDF pe
              s = s + 1
              stat(s) = NF_GET_VARA_REAL(ncid_in, fproperties(b)% fid, startv, countv, field_3D)
            endif ! netCDF pe
            CALL MPI_BCAST(field_3D, dim_fields3D, mpi_float, 0, mpi_comm_world, mpierror)
            ! fill local state vector:
            do i = 1, dim_fields_l(b)
              state_l(i + offsets(b)) = REAL(field_3D(my_nodlist_3D(i)),8)
            enddo
          endif
        enddo

        ! Write local state vector
        startv(2) = startpos
        countv(2) = 1
        startv(1) = 1
        countv(1) = dim_state_l
        s = s + 1
        stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_out, startv, countv, state_l)

        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in writing to output file, no.', i, &
                ' pe, irank',pe, irank
        END DO

     END DO loopfields

     s = s + 1
     stat(s) = nf_close(ncid_out)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in closing output file, no.', i, ' pe ', pe
     END DO

     DEALLOCATE(myfield_surf, myfield_3D)
     DEALLOCATE(state_l, my_nodlist_surf, my_nodlist_3D)

  if (pe==0) stat(1) = nf_close(ncid_in)

  if (pe==0) WRITE (*,'(/1x,a/)') '------- END -------------'
  CALL MPI_FINALIZE(mpierror)

END PROGRAM distribute_covar
