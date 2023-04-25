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

! !USES:
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'


  ! Local variables
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
  INTEGER :: dimid_2D, dimid_tracer3D, dimid_w                          ! NC IDs for pe-distributed outfiles
  INTEGER :: id_mu, id_mv, id_mw, id_mz, id_ms, id_mt, id_mi            ! NC IDs global running mean state variables
  INTEGER :: id_svdu, id_svdv, id_svdw, id_svdz, id_svds, id_svdt, &
             id_svdi                                                    ! NC IDs global scaled singular vectors variables
  INTEGER :: id_u, id_v, id_w, id_z, id_t, id_s, id_i, id_out
  INTEGER :: startpos
  INTEGER :: rank, dim_state, nfields
  INTEGER :: nlevels                                                    ! Number of vertical levels
  INTEGER :: nod2D, elem2D
  INTEGER :: stat(100)
  INTEGER :: dim_state_l                                                ! length of pe-local state vector
  INTEGER :: countv(2), startv(2)
  INTEGER :: dimids(2)
  INTEGER :: npes                                                       ! number of model processes
  INTEGER :: myDim_nod2D, myDim_elem2D                                  ! number of pe-local nodes / elements
  INTEGER :: eDim_nod2D, eDim_elem2D                                    ! number of pe-local halo nodes / elements
  INTEGER :: eXDim_elem2D
  INTEGER, ALLOCATABLE :: nod2D_part(:), elem2D_part(:)                 ! rank partioning vectors (list of nodes, elems per process)
  REAL(kind=4), ALLOCATABLE :: field_2D(:), field_uv(:)
  REAL(kind=4), ALLOCATABLE :: field_w(:)  , field_ts(:)
  REAL(kind=4), ALLOCATABLE :: svals(:)
  REAL(kind=4), ALLOCATABLE :: myfield_2D(:), myfield_uv(:)             ! pe-local field vectors
  REAL(kind=4), ALLOCATABLE :: myfield_w(:)  , myfield_ts(:)            ! pe-local field vectors
  REAL(kind=8), ALLOCATABLE :: state_l(:)                               ! pe-local state vector
  INTEGER, ALLOCATABLE :: List_nod2D(:), List_elem2D(:)                 ! (dummy) arrays to read from myList files
  INTEGER, ALLOCATABLE :: myList_nod2D(:), myList_elem2D(:)             ! (dummy) arrays to read from myList files
  INTEGER, ALLOCATABLE :: my_nodlist_2D(:), my_nodlist_uv(:)            ! pe-local node lists (repeated lists for 3D variables)
  INTEGER, ALLOCATABLE :: my_nodlist_w(:)  , my_nodlist_ts(:)           ! pe-local node lists (repeated lists for 3D variables)
  
  INTEGER, ALLOCATABLE :: dim_fields(:)                                              ! Field dimensions for SSH, u, v, w, temp, salt
  INTEGER, ALLOCATABLE :: offsets(:)                                                 ! Field offsets in state vector
  INTEGER, ALLOCATABLE :: dim_fields_l(:)                                            ! dimensions of pe-local fields
  
  INTEGER :: biomin, biomax ! Start and end indices of biogeochemistry fields in state vector
  
  ! Declare Fortran type holding the indices of model fields in the state vector
  TYPE field_ids
     INTEGER :: ssh
     INTEGER :: u 
     INTEGER :: v 
     INTEGER :: w 
     INTEGER :: temp 
     INTEGER :: salt
     INTEGER :: a_ice
     INTEGER :: MLD1
     INTEGER :: PhyChl
     INTEGER :: DiaChl
     INTEGER :: DIC
     INTEGER :: DOC
     INTEGER :: Alk
     INTEGER :: DIN
     INTEGER :: DON
     INTEGER :: O2
     INTEGER :: pCO2s
     INTEGER :: CO2f
  END TYPE field_ids
  ! Type variable holding field IDs in state vector
  TYPE(field_ids) :: ids
  
  ! Declare Fortran type holding model field-specific variabels
  type state_field
   integer :: ndims = 0                   ! Number of field dimensions (1 or 2)
   character(len= 10) :: variable = ''    ! Name of field
   integer :: fid = 0                     ! Field ID (in/output file)
   integer :: mid = 0                     ! Mean field ID (output file)
   integer :: sid = 0                     ! Singular vector ID
  end type state_field

  type(state_field), allocatable :: biofields(:) ! Type variable holding the
                                                 ! definitions of model fields

! ************************************************
! *** Configuration                            ***
! ************************************************

  ! Path to and name of file holding global covariance matrix
  ! inpath = '/work/ollie/frbunsen/model_runs/fesom2/test_control/cov/'
  inpath = '/albedo/work/projects/p_recompdaf/frbunsen/modelruns/cov/'
  infile = 'cov.nc'

  ! Path to mesh
  ! dist_mesh_dir = '/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/dist_72/'
  dist_mesh_dir = '/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/dist_72/'
  partfile      = 'rpart.out'
  mylistfile    = 'my_list'

  ! Path to and name stub of output files
  outpath = '/albedo/work/projects/p_recompdaf/frbunsen/modelruns/cov/dist72/'
  outfile = 'cov'
  
  ! Composition of state vector:
  nfields = 18
  
  ids% ssh    = 1
  ids% u      = 2
  ids% v      = 3
  ids% w      = 4
  ids% temp   = 5
  ids% salt   = 6
  ids% a_ice  = 7
  ids% MLD1   = 8
  ids% PhyChl = 9
  ids% DiaChl = 10
  ids% DIC    = 11
  ids% DOC    = 12
  ids% Alk    = 13
  ids% DIN    = 14
  ids% DON    = 15
  ids% O2     = 16
  ids% pCO2s  = 17
  ids% CO2f   = 18
  
  biomin = 8
  biomax = 18
  
  ! Field-specific variables:
  allocate(biofields(nfields))
  
  biofields (ids% MLD1) % ndims = 1
  biofields (ids% MLD1) % variable = 'MLD1'

  biofields (ids% PhyChl) % ndims = 2 
  biofields (ids% PhyChl) % variable = 'PhyChl'

  biofields (ids% DiaChl) % ndims = 2
  biofields (ids% DiaChl) % variable = 'DiaChl'

  biofields(ids% DIC) % ndims = 2
  biofields(ids% DIC) % variable = 'DIC'

  biofields(ids% DOC) % ndims = 2
  biofields(ids% DOC) % variable = 'DOC'

  biofields(ids% Alk) % ndims = 2
  biofields(ids% Alk) % variable = 'Alk'

  biofields(ids% DIN) % ndims = 2
  biofields(ids% DIN) % variable = 'DIN'

  biofields(ids% DON) % ndims = 2
  biofields(ids% DON) % variable = 'DON'

  biofields(ids% O2) % ndims = 2
  biofields(ids% O2) % variable = 'O2'

  biofields(ids% pCO2s) % ndims = 1
  biofields(ids% pCO2s) % variable = 'pCO2s'

  biofields(ids% CO2f) % ndims = 1
  biofields(ids% CO2f) % variable = 'CO2f'

! ************************************************
! *** Init                                     ***
! ************************************************

  WRITE (*,'(3x,a)')  '*****************************************************************'
  WRITE (*,'(3x,a)')  '*** Generate distributed covariance matrix files for FESOM2.1 ***'
  WRITE (*,'(3x,a/)') '*****************************************************************'

  ncfile_in = TRIM(inpath)//TRIM(infile)
  WRITE (*,*) 'Read fields from file: ',TRIM(ncfile_in)

  ncfile_out = TRIM(outpath)//TRIM(outfile)
  WRITE (*,*) 'Write distributed fields to files: ',TRIM(ncfile_out),'_XXXX.nc'
  WRITE (*,*) 'Read partitioning information from directory: ',TRIM(dist_mesh_dir)


! ********************************************
! *** Open global file and read dimensions ***
! ********************************************

  s = 1
  stat(s) = NF_OPEN(TRIM(ncfile_in), NF_NOWRITE, ncid_in)

  ! Get dimensions
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nfields', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, nfields)
  
  allocate(dim_fields(nfields))
  allocate(offsets(nfields))
  allocate(dim_fields_l(nfields))
  
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, rank)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nod2', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_fields(1)) ! for SSH (SIC)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nz1_x_nod2', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_fields(4)) ! for w
    s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nz_x_nod2', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_fields(5)) ! for temp (salt, u, v)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_state)
  
  ! dimensions of u, v and salt are like temp (5)
  dim_fields(2) = dim_fields(5) ! u
  dim_fields(3) = dim_fields(5) ! v
  dim_fields(6) = dim_fields(5) ! salt
  ! dimensions of SIC are like SSH
  dim_fields(ids% a_ice) = dim_fields(ids% SSH)
  
  ! dimensions of biogeochemistry fields
  do b=biomin, biomax
    if (biofields(b)% ndims == 1) dim_fields(b) = dim_fields(ids% SSH)  ! surface fields
    if (biofields(b)% ndims == 2) dim_fields(b) = dim_fields(ids% temp) ! 3D fields
  enddo

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
  END DO
  
  nlevels = 47         ! number of vertical levels
  nod2D   = 126858     ! number of nodes
  elem2D  = 244659     ! number of elements (not used)
  
  ! Write dimensions
  WRITE (*,'(/1x,a)') 'Global dimensions of experiment:'
  WRITE (*,'(10x,1x,a30,i12)') 'dim_fields 1,7  (SSH, SIC)    ', dim_fields(1)
  WRITE (*,'(10x,1x,a30,i12)') 'dim_fields 2,3  (u, v)        ', dim_fields(2)
  WRITE (*,'(10x,1x,a30,i12)') 'dim_fields 4    (w)           ', dim_fields(4)
  WRITE (*,'(10x,1x,a30,i12)') 'dim_fields 5,6  (temp, salt)  ', dim_fields(5)
  
  do b=biomin, biomax
    WRITE (*,'(10x,1x,a30,i12)') 'dim_fields '//trim(biofields(b)% variable), dim_fields(b)
  enddo
  
  WRITE (*,'(10x,1x,a30,i12)') 'nlevels                       ', nlevels
  WRITE (*,'(10x,1x,a30,i12)') 'rank                          ', rank
  WRITE (*,'(10x,1x,a30,i12)') 'dim_state                     ', dim_state
  WRITE (*,'(10x,1x,a30,i12)') 'nfields                       ', nfields
  
  IF (dim_fields(1)*nlevels /= dim_fields(5)) THEN
	WRITE (*,*) 'Number of vertical levels and horizontal nodes not consistent with 3D nodes!'
	STOP
  END IF

! ***********************************
! *** Read mesh partitioning data ***
! ***********************************
 
  OPEN(unit=1,file=TRIM(dist_mesh_dir)//TRIM(partfile), status='old', form='formatted')
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
	file_name = trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
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
  IF (SUM(nod2D_part) /= nod2D) THEN
	WRITE (*,*) 'Partitioned mesh not consistent with global mesh (nodes)'
	STOP
  END IF
  
  ! elem2D_part(i) is the number of triangles that belong to each PE,
  ! which are those that contain at least one vertex that belongs to the PE.
  ! Thus, triangles with vertices that belong to several PE are counted
  ! as my_triangle by each PE.
  
  WRITE(*,*) 'Sum of elements in partitioned mesh: ', SUM(elem2D_part)
  WRITE(*,*) 'Sum of elements in global mesh:      ', elem2D

  WRITE (*,'(/1x,a)') 'Mesh distribution:'
  WRITE (*,'(5x,a,i6)') 'Number of PEs:', npes
  WRITE (*,'(5x,a,4x,a)') 'Local nodes:   nod2d','elem2d'
  WRITE (*,*) 'Local nodes:'
  DO i = 1, npes
     WRITE (*,'(18x,i7,i10)') nod2D_part(i) , elem2D_part(i)
  END DO  



! ************************************
! *** Initialize distributed files ***
! ************************************

  WRITE (*,'(/1x,a/)') '------- Initialize distributed files -------------'

  ! Read singular values from input file

  ALLOCATE(svals(rank))

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

  ! Initialize one file for each PE
  initoutfiles: DO pe = 0, npes - 1
  
     dim_fields_l(1)          = nod2D_part(pe+1)                ! SSH
     dim_fields_l(2)          = nod2D_part(pe+1) * nlevels      ! u
     dim_fields_l(3)          = nod2D_part(pe+1) * nlevels      ! v
     dim_fields_l(4)          = nod2D_part(pe+1) * (nlevels+1)  ! w
     dim_fields_l(5)          = nod2D_part(pe+1) * nlevels      ! temp
     dim_fields_l(6)          = nod2D_part(pe+1) * nlevels      ! salt
     dim_fields_l(ids% a_ice) = nod2D_part(pe+1)
     
     do b=biomin,biomax
       if (biofields(b)% ndims == 1) dim_fields_l(b) = nod2D_part(pe+1)
       if (biofields(b)% ndims == 2) dim_fields_l(b) = nod2D_part(pe+1) * nlevels
     enddo

     dim_state_l = SUM(dim_fields_l)
     
     write(mype_string,'(i4.4)') pe
     
     WRITE (*,'(a, i4, i7)') 'pe and dim_state_l: ', pe, dim_state_l

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
	  stat(s) = NF_DEF_DIM(ncid_out, 'nod2', dim_fields_l(1), dimid_2D)     ! dim for SSH, SIC
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'nz1_x_nod2', dim_fields_l(4), dimid_w) ! dim for w
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'nz_x_nod2', dim_fields_l(5), dimid_tracer3D) ! dim for temp, salt, u, v
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
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_2D', NF_INT, 1, dimid_2D, id_list)  ! nodes for SSH, SIC
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_w',   NF_INT, 1, dimid_w,   id_list)  ! nodes for w
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_tracer3D',  NF_INT, 1, dimid_tracer3D,  id_list)  ! nodes for temp, salt, u, v

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

  END DO initoutfiles


! ****************************************************
! *** Generate and write distributed state vectors ***
! ****************************************************

  WRITE (*,'(1x,a/)') '------- Generate and write distributed files -------------'

    ! Allocate global fields
     ALLOCATE(field_2D (dim_fields(1)))
     ALLOCATE(field_uv (dim_fields(2)))
     ALLOCATE(field_w  (dim_fields(4)))
     ALLOCATE(field_ts (dim_fields(5)))

  ! *** Read, distribute and write mean state

  peloop: DO pe = 0, npes - 1

     WRITE (*,'(5x, a, i5)') '-- Process ', pe
     
     dim_fields_l(1)          = nod2D_part(pe+1)                ! SSH
     dim_fields_l(2)          = nod2D_part(pe+1) * nlevels      ! u
     dim_fields_l(3)          = nod2D_part(pe+1) * nlevels      ! v
     dim_fields_l(4)          = nod2D_part(pe+1) * (nlevels+1)  ! w
     dim_fields_l(5)          = nod2D_part(pe+1) * nlevels      ! temp
     dim_fields_l(6)          = nod2D_part(pe+1) * nlevels      ! salt
     dim_fields_l(ids% a_ice) = nod2D_part(pe+1)
     
     do b=biomin,biomax
       if (biofields(b)% ndims == 1) dim_fields_l(b) = nod2D_part(pe+1)
       if (biofields(b)% ndims == 2) dim_fields_l(b) = nod2D_part(pe+1) * nlevels
     enddo

	 ! Size of local state vector
     dim_state_l = SUM(dim_fields_l)

     ! Allocate fields
	 ALLOCATE(myfield_2D (dim_fields_l(1)))
	 ALLOCATE(myfield_uv (dim_fields_l(2)))
	 ALLOCATE(myfield_w  (dim_fields_l(4)))
	 ALLOCATE(myfield_ts (dim_fields_l(5)))
	 
     ALLOCATE(state_l(dim_state_l))

     ! Define offsets
     offsets (1)          = 0                            ! offset of field SSH
     offsets (2)          = dim_fields_l(1)              ! offset of field u
     offsets (3)          = offsets(2) + dim_fields_l(2) ! offset of field v
     offsets (4)          = offsets(3) + dim_fields_l(3) ! offset of field w
     offsets (5)          = offsets(4) + dim_fields_l(4) ! offset of field t
     offsets (6)          = offsets(5) + dim_fields_l(5) ! offset of field s
     offsets (ids% a_ice) = offsets(6) + dim_fields_l(6)
     
     do b=biomin,biomax
       offsets(b) = offsets(b-1) + dim_fields_l(b-1)
     enddo

    ! Read my_listXXXXX files
    WRITE (*,'(1x,a)') '*** Read my_listXXXXX files ***'
    
	write(mype_string,'(i5.5)') pe
	file_name = trim(dist_mesh_dir)//'my_list'//trim(mype_string)//'.out'  
	fileID = pe
		
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
	
	WRITE (*,'(1x,a)') '*** Create list of 3D nodes / elements ***'
	
	allocate(my_nodlist_2D(dim_fields_l(1)))
	allocate(my_nodlist_uv(dim_fields_l(2)))
	allocate(my_nodlist_w  (dim_fields_l(4)))
	allocate(my_nodlist_ts (dim_fields_l(5)))
	
	my_nodlist_2D = myList_nod2D

	! indeces of pe-local nodes / elements in global variable vector for 3D-fields
    ! velocity w
    DO k = 1, nlevels+1
      DO j = 1, myDim_nod2D
         my_nodlist_w((j-1)*(nlevels+1)+k) = (myList_nod2D(j)-1) * (nlevels+1) + k
      END DO
    END DO
    
	! temp and salinity
	! velocities u and v
	DO k = 1, nlevels
      DO j = 1, myDim_nod2D
         my_nodlist_ts((j-1)*nlevels+k) = (myList_nod2D(j)-1) * nlevels + k
         my_nodlist_uv((j-1)*nlevels+k) = (myList_nod2D(j)-1) * nlevels + k
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
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_2D', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, my_nodlist_2D)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_w', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, my_nodlist_w)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_tracer3D', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, my_nodlist_ts)
     
     ! Inquire IDs for mean state and singular vectors
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'running_meanstate', id_mean)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'V', id_svec)


     ! Get IDs for global fields from input file
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'ssh_mean', id_mz)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'u_mean', id_mu)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'v_mean', id_mv)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'w_mean', id_mw)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'temp_mean', id_mt)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'salt_mean', id_ms)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'a_ice_mean', id_mi)
     
     do b=biomin,biomax
       s = s + 1
       stat(s) = NF_INQ_VARID(ncid_in, trim(biofields(b)%variable)//'_mean', biofields(b)% mid)
     enddo
     
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'ssh_svd', id_svdz)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'u_svd', id_svdu)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'v_svd', id_svdv)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'w_svd', id_svdw)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'temp_svd', id_svdt)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'salt_svd', id_svds)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'a_ice_svd', id_svdi)
     
     do b=biomin,biomax
       s = s + 1
       stat(s) = NF_INQ_VARID(ncid_in, trim(biofields(b)%variable)//'_svd', biofields(b)% sid)
     enddo

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in inquire from output file, no.', i, ' pe ',pe
     END DO

     loopfields: DO irank = 0, rank

        IF (irank == 0) THEN
           ! Treat mean state

           ! Ids for input
           id_z = id_mz
           id_u = id_mu
           id_v = id_mv
           id_w = id_mw
           id_t = id_mt
           id_s = id_ms
           id_i = id_mi
           
           do b=biomin,biomax
             biofields(b)% fid = biofields(b)% mid
           enddo

           ! Id for output
           id_out = id_mean
           
           ! Column of matrix in file
           startpos = 1
           
        ELSE
           ! Treat singular vectors

           ! Ids for input
           id_z = id_svdz
           id_u = id_svdu
           id_v = id_svdv
           id_w = id_svdw
           id_t = id_svdt
           id_s = id_svds
           id_i = id_svdi
           
           do b=biomin,biomax
             biofields(b)% fid = biofields(b)% sid
           enddo

           ! Id for output
           id_out = id_svec
           
           ! Column of matrix in file
           startpos = irank
           
        END IF

        WRITE (*,'(1x,a/)') '*** Generate local state vector ***'
        ! *** Generate local state vector

        ! Read global SSH (mean state or singular vectors, respectively)
        startv(2) = startpos
        countv(2) = 1
        startv(1) = 1
        countv(1) = dim_fields(1)
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_z, startv, countv, field_2D)

        ! Initialize local nodes of pe-local state vector
        DO i = 1, dim_fields_l(1)
              state_l(i + offsets(1)) = field_2D(my_nodlist_2D(i))
        END DO

        ! Read global U
        countv(1) = dim_fields(2) 
        s = 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_u, startv, countv, field_uv)

        ! Initialize local nodes of local state vector
        DO i = 1, dim_fields_l(2)
           state_l(i + offsets(2)) = REAL(field_uv(my_nodlist_uv(i)),8)
        END DO

        ! Read global V
        countv(1) = dim_fields(3)
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_v, startv, countv, field_uv)

        ! Initialize local nodes of local state vector
        DO i = 1, dim_fields_l(3)
           state_l(i + offsets(3)) = REAL(field_uv(my_nodlist_uv(i)),8)
        END DO

        ! Read global W
        countv(1) = dim_fields(4)
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_w, startv, countv, field_w)

        ! Initialize local nodes of local state vector
        DO i = 1, dim_fields_l(4)
           state_l(i + offsets(4)) = REAL(field_w(my_nodlist_w(i)),8)
        END DO

        ! Read global Temperature
        countv(1) = dim_fields(5)
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_t, startv, countv, field_ts)

        ! Initialize local nodes of local state vector
        DO i = 1, dim_fields_l(5)
           state_l(i + offsets(5)) = REAL(field_ts(my_nodlist_ts(i)),8)
        END DO

        ! Read global Salinity
        countv(1) = dim_fields(6)
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_s, startv, countv, field_ts)

        ! Initialize local nodes of local state vector
        DO i = 1, dim_fields_l(6)
           state_l(i + offsets(6)) = REAL(field_ts(my_nodlist_ts(i)),8)
        END DO
        
        ! Read global SIC
        countv(1) = dim_fields(ids% a_ice)
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_i, startv, countv, field_2D)

        ! Initialize local nodes of pe-local state vector
        DO i = 1, dim_fields_l(ids% a_ice)
              state_l(i + offsets(ids% a_ice)) = REAL(field_2D(my_nodlist_2D(i)),8)
        END DO
        
        ! Read global biogeochemistry
        do b=biomin,biomax
        
          countv(1) = dim_fields(b)
          
          if (biofields(b)% ndims==1) then
            ! surface fields
            ! read global field
            s = s + 1
            stat(s) = NF_GET_VARA_REAL(ncid_in, biofields(b)% fid, startv, countv, field_2D)
            ! fill local state vector:
             do i = 1, dim_fields_l(b)
              state_l(i + offsets(b)) = REAL(field_2D(my_nodlist_2D(i)),8)
             enddo
             
          elseif (biofields(b)% ndims==2) then
            ! 3D-fields
            ! read global field
            s = s + 1
            stat(s) = NF_GET_VARA_REAL(ncid_in, biofields(b)% fid, startv, countv, field_ts)
            ! fill local state vector:
            do i = 1, dim_fields_l(b)
              state_l(i + offsets(b)) = REAL(field_ts(my_nodlist_ts(i)),8)
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
             WRITE(*, *) 'NetCDF error in closing output file, no.', i, ' pe ',pe
     END DO

     DEALLOCATE(myfield_2D, myfield_uv, myfield_w, myfield_ts)
     DEALLOCATE(state_l, my_nodlist_2D, my_nodlist_uv, my_nodlist_w, my_nodlist_ts)

  END DO peloop

  stat(1) = nf_close(ncid_in)

  WRITE (*,'(/1x,a/)') '------- END -------------'

END PROGRAM distribute_covar
