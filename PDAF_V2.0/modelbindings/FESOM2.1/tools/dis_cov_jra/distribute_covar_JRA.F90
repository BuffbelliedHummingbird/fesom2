! PDAF_V2.0

! $Id: distribute_covar.F90 1604 2016-05-30 06:42:16Z lnerger $

! Modified by QT     2018-02-27 ---- for AWI-CM
! Modified by Frauke 2022-02-28 ---- for FESOM2.1 with core2 mesh
! Modified by Frauke 2022-11-17 ---- for JRA atmospheric forcing

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

! !USES:
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'


  ! Local variables
  INTEGER :: i, s, iter, pe, irank, n                                   ! counters
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
  INTEGER :: dimid_rank, dimid_nfields, dimid_state                     ! NC IDs for pe-distributed outfiles
  INTEGER :: dimid_nod2                                                 ! NC IDs for pe-distributed outfiles
  INTEGER :: id_mhuss, id_mprra, id_mprsn, id_mpsl, id_mrlds, id_mrsds  ! NC IDs global running mean state variables
  INTEGER :: id_mtas, id_muas, id_mvas
  INTEGER :: id_svdhuss, id_svdprra, id_svdprsn, id_svdpsl, id_svdrlds, &
             id_svdrsds, id_svdtas, id_svduas, id_svdvas                ! NC IDs global scaled singular vectors variables
  INTEGER :: id_huss, id_prra, id_prsn, id_psl, id_rlds, id_rsds, &
             id_tas, id_uas, id_vas
  INTEGER :: startpos, id_out
  INTEGER :: rank, dim_state, nfields
  INTEGER :: dim_fields                                                 ! Field dimensions (all fields equal)
  INTEGER :: offsets(9)                                                 ! Field offsets in state vector
  INTEGER :: nlevels                                                    ! Number of vertical levels
  INTEGER :: nod2D, elem2D
  INTEGER :: stat(100)
  INTEGER :: dim_state_l                                                ! length of pe-local state vector
  INTEGER :: dim_fields_l                                               ! dimensions of pe-local fields
  INTEGER :: countv(2), startv(2)
  INTEGER :: dimids(2)
  INTEGER :: npes                                                       ! number of model processes
  INTEGER :: myDim_nod2D, myDim_elem2D                                  ! number of pe-local nodes / elements
  INTEGER :: eDim_nod2D, eDim_elem2D                                    ! number of pe-local halo nodes / elements
  INTEGER :: eXDim_elem2D
  INTEGER, ALLOCATABLE :: nod2D_part(:), elem2D_part(:)                 ! rank partioning vectors (list of nodes, elems per process)
  REAL(kind=4), ALLOCATABLE :: field_2D(:)
  REAL(kind=8), ALLOCATABLE :: svals(:)
  REAL(kind=4), ALLOCATABLE :: myfield_2D(:)                            ! pe-local field vectors
  REAL(kind=8), ALLOCATABLE :: state_l(:)                               ! pe-local state vector
  INTEGER, ALLOCATABLE :: List_nod2D(:), List_elem2D(:)                 ! (dummy) arrays to read from myList files
  INTEGER, ALLOCATABLE :: myList_nod2D(:), myList_elem2D(:)             ! (dummy) arrays to read from myList files
  ! Declare Fortran type holding the indices of model fields in the state vector
  TYPE field_ids
     INTEGER :: huss
     INTEGER :: prra
     INTEGER :: prsn
     INTEGER :: psl
     INTEGER :: rlds 
     INTEGER :: rsds
     INTEGER :: tas
     INTEGER :: uas
     INTEGER :: vas
  END TYPE field_ids
  ! Type variable holding field IDs in state vector
  TYPE(field_ids) :: ids

! ************************************************
! *** Configuration                            ***
! ************************************************

  ! Path to and name of file holding global covariance matrix
  ! inpath = '/work/ollie/frbunsen/jra_core2/cov/'
  inpath = '/albedo/work/projects/p_recompdaf/frbunsen/data/jra_core2/cov/'
  infile = 'cov.nc'

  ! Path to mesh
  ! dist_mesh_dir = '/work/ollie/projects/clidyn/FESOM2/meshes/core2/dist_72/'
  dist_mesh_dir = '/albedo/work/projects/p_recompdaf/frbunsen/FESOM2/meshes/core2/dist_126/'
  partfile      = 'rpart.out'
  mylistfile    = 'my_list'

  ! Path to and name stub of output files
  ! outpath = '/work/ollie/frbunsen/jra_core2/dist72_cov/'
  outpath = '/albedo/work/projects/p_recompdaf/frbunsen/data/jra_core2/dist126_cov/'
  outfile = 'cov'
  
  ! Composition of state vector:
  ids% huss = 1
  ids% prra = 2
  ids% prsn = 3
  ids% psl  = 4
  ids% rlds = 5
  ids% rsds = 6
  ids% tas  = 7
  ids% uas  = 8
  ids% vas  = 9

! ************************************************
! *** Init                                     ***
! ************************************************

  WRITE (*,'(3x,a)')  '************************************************************'
  WRITE (*,'(3x,a)')  '*** Generate distributed covariance matrix files for JRA ***'
  WRITE (*,'(3x,a/)') '************************************************************'

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
  stat(s) = NF_INQ_DIMID(ncid_in, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, rank)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nod2', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_fields)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, dim_state)
  s = s + 1
  stat(s) = NF_INQ_DIMID(ncid_in, 'nfields', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(ncid_in, id_dim, nfields)
  
  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
  END DO
  
  nlevels = 47         ! number of vertical levels
  nod2D   = 126858     ! number of nodes
  elem2D  = 244659     ! number of elements (not used)
  
  ! Write dimensions
  WRITE (*,'(/1x,a)') 'Global dimensions of experiment:'
  WRITE (*,'(10x,1x,a25,i12)') 'dim_fields                    ', dim_fields
  WRITE (*,'(10x,1x,a25,i12)') 'nlevels                       ', nlevels
  WRITE (*,'(10x,1x,a25,i12)') 'rank                          ', rank
  WRITE (*,'(10x,1x,a25,i12)') 'dim_state                     ', dim_state
  WRITE (*,'(10x,1x,a25,i12)') 'nfields                       ', nfields

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
  stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_sigma, svals)
  s = s + 1
  title = ''
  stat(s) = NF_GET_ATT_TEXT(ncid_in, NF_GLOBAL, 'title', title)
  s = s + 1
  fieldsstr = ''
  stat(s) = NF_GET_ATT_TEXT(ncid_in, NF_GLOBAL, 'state_fields', fieldsstr)
  
  WRITE(*,*) 'Sigma: '
  WRITE(*,*) svals(1:5)
  WRITE(*,*) '...'
  WRITE(*,*) ''

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading from input file, no.', i
  END DO

  ! Initialize one file for each PE
  initoutfiles: DO pe = 0, npes - 1
  
     dim_fields_l = nod2D_part(pe+1)
     dim_state_l  = dim_fields_l * nfields
     
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
	  stat(s) = NF_DEF_DIM(ncid_out, 'nod2', dim_fields_l, dimid_nod2)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'dim_state', dim_state_l, dimid_state)
	  s = s + 1
	  stat(s) = NF_DEF_DIM(ncid_out, 'nfields',  nfields, dimid_nfields)

     ! Define variables

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'sigma', NF_DOUBLE, 1, dimid_rank, id_sigma)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nodlist_2D', NF_INT, 1, dimid_nod2, id_list)

     ! running mean state (for the last snap shot)

     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'running_meanstate', NF_DOUBLE, 1, dimid_state, id_mean)

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
  ALLOCATE(field_2D (dim_fields))

  ! *** Read, distribute and write mean state

  peloop: DO pe = 0, npes - 1

     WRITE (*,'(5x, a, i5)') '-- Process ', pe
     
     dim_fields_l = nod2D_part(pe+1)

	 ! Size of local state vector
     dim_state_l = dim_fields_l * nfields

     ! Allocate fields
	 ALLOCATE(myfield_2D (dim_fields_l))
	 
     ALLOCATE(state_l(dim_state_l))

  ! Define offsets
  offsets (1) = 0                        
  offsets (2) = offsets(1) + dim_fields_l        
  offsets (3) = offsets(2) + dim_fields_l        
  offsets (4) = offsets(3) + dim_fields_l       
  offsets (5) = offsets(4) + dim_fields_l       
  offsets (6) = offsets(5) + dim_fields_l        
  offsets (7) = offsets(6) + dim_fields_l        
  offsets (8) = offsets(7) + dim_fields_l       
  offsets (9) = offsets(8) + dim_fields_l        

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

	deallocate(myList_elem2D)

     ! Open output file
     write(mype_string,'(i4.4)') pe
     
     s = 1
     stat(s) = NF_OPEN(TRIM(ncfile_out)//'_'//TRIM(mype_string)//'.nc', NF_WRITE, ncid_out)

     ! Write MyLists
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'nodlist_2D', id_list)
     s = s + 1
     stat(s) = NF_PUT_VAR_INT(ncid_out, id_list, myList_nod2D)
     
     ! Inquire IDs for mean state and singular vectors
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'running_meanstate', id_mean)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_out, 'V', id_svec)


     ! Get IDs for global fields from input file
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'huss_mean', id_mhuss)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'prra_mean', id_mprra)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'prsn_mean', id_mprsn)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'psl_mean',  id_mpsl)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'rlds_mean', id_mrlds)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'rsds_mean', id_mrsds)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'tas_mean',  id_mtas)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'uas_mean',  id_muas)
          s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'vas_mean',  id_mvas)

     
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'huss_svd', id_svdhuss)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'prra_svd', id_svdprra)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'prsn_svd', id_svdprsn)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'psl_svd', id_svdpsl)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'rlds_svd', id_svdrlds)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'rsds_svd', id_svdrsds)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'tas_svd', id_svdtas)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'uas_svd', id_svduas)
     s = s + 1
     stat(s) = NF_INQ_VARID(ncid_in, 'vas_svd', id_svdvas)


     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in inquire from output file, no.', i, ' pe ',pe
     END DO

     loopfields: DO irank = 0, rank

        IF (irank == 0) THEN
           ! Treat mean state

           ! Ids for input
           id_huss = id_mhuss
           id_prra = id_mprra
           id_prsn = id_mprsn
           id_psl  = id_mpsl
           id_rlds = id_mrlds
           id_rsds = id_mrsds
           id_tas  = id_mtas
           id_uas  = id_muas
           id_vas  = id_mvas

           ! Id for output
           id_out = id_mean
           
           ! Column of matrix in file
           startpos = 1
           
        ELSE
           ! Treat singular vectors

           ! Ids for input
           id_huss = id_svdhuss
           id_prra = id_svdprra
           id_prsn = id_svdprsn
           id_psl  = id_svdpsl
           id_rlds = id_svdrlds
           id_rsds = id_svdrsds
           id_tas  = id_svdtas
           id_uas  = id_svduas
           id_vas  = id_svdvas

           ! Id for output
           id_out = id_svec
           
           ! Column of matrix in file
           startpos = irank
           
        END IF

        WRITE (*,'(1x,a/)') '*** Generate local state vector ***'
        ! *** Generate local state vector

        startv(1) = 1
        countv(1) = dim_fields
        startv(2) = startpos
        countv(2) = 1

        ! Read global field (mean state or singular vectors, respectively)
        s = 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_huss, startv, countv, field_2D)
        ! Fill state vector with local field:
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% huss)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_prra, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% prra)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_prsn, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% prsn)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_psl, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% psl )) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_rlds, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% rlds)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_rsds, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% rsds)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_tas, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% tas)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_uas, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% uas)) = field_2D(myList_nod2D(i))
        END DO
        
        s = s + 1
        stat(s) = NF_GET_VARA_REAL(ncid_in, id_vas, startv, countv, field_2D)
        DO i = 1, dim_fields_l
              state_l(i + offsets(ids% vas)) = field_2D(myList_nod2D(i))
        END DO

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

     DEALLOCATE(myfield_2D)
     DEALLOCATE(state_l, myList_nod2D)

  END DO peloop

  stat(1) = nf_close(ncid_in)

  WRITE (*,'(/1x,a/)') '------- END -------------'

END PROGRAM distribute_covar
