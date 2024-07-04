! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022-02 - Frauke B     - Adapted for FESOM2.1

! !USES:
  USE mod_assim_pdaf, &
       ONLY: file_init, path_init, read_inistate, file_inistate, varscale, &
       offset, ASIM_START_USE_CLIM_STATE, this_is_pdaf_restart, mesh_fesom, &
       id, dim_fields, offset, start_from_ENS_spinup
  USE mod_nc_out_variables, &
       ONLY: sfields
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, COMM_filter, abort_parallel
  USE g_PARSUP, &
       ONLY: MPIerr, edim_nod2d, mydim_nod2d
  USE o_arrays, &
       ONLY: eta_n, tr_arr, uv, wvel
  USE i_arrays, &
       ONLY: a_ice
  USE g_ic3d
  USE recom_config, ONLY: tiny

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'
  INCLUDE 'mpif.h'

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

! !CALLING SEQUENCE:
! Called by: PDAF_init       (as U_init_ens)

! *** local variables ***
  INTEGER :: i, j, member, s, row, col, k ! Counters
  INTEGER :: rank                      ! Rank stored in init file
  INTEGER :: dim_p_cov                 ! Local state dimension in covariance file
  INTEGER,ALLOCATABLE :: offsets_cov(:)! Local field offsets in covariance file
  INTEGER :: fileid                    ! ID for NetCDF file
  INTEGER :: id_dim                    ! ID for dimension
  INTEGER :: id_state,id_svals,id_eof  ! IDs for fields
  INTEGER :: startv(2),countv(2)       ! Vectors for reading fields
  INTEGER :: stat(7)                   ! Status flag for NetCDF commands
  REAL :: fac                          ! Square-root of dim_ens or dim_ens-1
  REAL,ALLOCATABLE :: eof_p(:,:)       ! Matrix of eigenvectors of covariance matrix
  REAL,ALLOCATABLE :: svals(:)         ! Singular values
  REAL,ALLOCATABLE :: omega(:,:)       ! Transformation matrix Omega
  CHARACTER(len=5)   :: mype_string    ! String for process rank
  CHARACTER(len=5)   :: col_string     ! String for ensemble member
  CHARACTER(len=150) :: infile         ! File holding initial state estimate
  INTEGER :: dim_p_read
  LOGICAL :: runningmean               ! True: Initialize state vector from
                                       ! nc-file running mean
  REAL, ALLOCATABLE :: ens_p_per(:,:)  ! Ensemble states of to-be-perturbed fields
  INTEGER :: o1,o2,f1,f2             ! Indeces
  INTEGER, ALLOCATABLE :: idxs(:)      ! Field indeces of to-be-perturbed model fields
  INTEGER :: nfields_per               ! Number of to-be-perturbed fields
  
  TYPE field_ids                       ! Field IDs for covariance matrix
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
  
  TYPE(field_ids) :: id_cov           ! ! Type variable holding field IDs in covariance matrix
  
  INTEGER :: n_treshold_ssh_p,  &
             n_treshold_temp_p, &
             n_treshold_salt_p, &
             n_treshold_sic_p          ! Local counters for treshold-based corrections
             
  INTEGER :: n_treshold_ssh_g,  &
             n_treshold_temp_g, &
             n_treshold_salt_g, &
             n_treshold_sic_g          ! Global counters for treshold-based corrections
             
             
             
! no need to initialize the ensemble in case of restart, just skip this routine:
  IF (this_is_pdaf_restart) THEN
      IF (mype_filter == 0)  WRITE(*,*) 'FESOM-PDAF This is a restart, skipping init_ens_pdaf'
      CALL collect_state_PDAF(dim_p, state_p)
  ELSEIF (start_from_ENS_spinup) THEN
      IF (mype_filter == 0)  WRITE(*,*) 'FESOM-PDAF Starting from a perturbed ensemble, skipping init_ens_pdaf'
      CALL collect_state_PDAF(dim_p, state_p)
  ELSE

! **********************
! *** INITIALIZATION ***
! **********************

  ! Fields in covariance input file:
  id_cov% ssh  = 1
  id_cov% u    = 2
  id_cov% v    = 3
  id_cov% temp = 4
  id_cov% salt = 5
  id_cov% DIC  = 6
  id_cov% Alk  = 7
  id_cov% DIN  = 8
  id_cov% O2   = 9
  
  ! Positions of to-be-perturbed state fields in the full state vector:
  nfields_per = 9
  ALLOCATE(idxs(nfields_per))
  idxs = (/ id%ssh,  &
            id%u,    &
            id%v,    &
            id%temp, &
            id%salt, &
            id%DIC,  &
            id%Alk,  &
            id%DIN,  &
            id%O2      /)
  
  ! Dimension of to-be-perturbed vector is
  ! sum of to to-be-perturbed field dimensions
  dim_p_cov =  0
  DO j=1, nfields_per
     dim_p_cov = dim_p_cov + dim_fields(idxs(j))
  END DO
  
  ! Offsets in covariance matrix:
  ALLOCATE(offsets_cov(nfields_per))
  offsets_cov(1) = 0
  do j=2, nfields_per
     offsets_cov(j) = offsets_cov(j-1) + dim_fields(idxs(j-1))
  enddo
  
  if (mype_filter==0) then
  do j=1,nfields_per
      WRITE (*,'(1X,A40,1X,A5,1X,I7,1X,I7)') 'PE0-dim and offset of covariance field ', &
                                              trim(sfields(idxs(j))%variable), &
                                              dim_fields(idxs(j)), &
                                              offsets_cov(j)
  enddo
  endif
  
  ! Print initialization message:
  write(mype_string,'(i4.4)') mype_filter
  mype0: IF (mype_filter == 0) THEN
    WRITE (*, '(/a, 8x,a)') 'FESOM-PDAF', 'Generate state ensemble from covariance matrix'
    WRITE (*, '(a, 8x,a)') &
         'FESOM-PDAF', '--- use 2nd order exact sampling (SEIK type)'
    WRITE (*, '(a, 8x,a,i5)') 'FESOM-PDAF', '--- number of EOFs:',dim_ens-1
  END IF mype0

  ! Allocate memory for temporary fields
  ALLOCATE(eof_p(dim_p_cov, dim_ens-1))
  ALLOCATE(svals(dim_ens-1))
  ALLOCATE(omega(dim_ens, dim_ens-1))

! *************************************************
! *** Initialize covar matrix                   ***
! *************************************************
 
  infile=Trim(path_init)//Trim(file_init)//TRIM(mype_string)//'.nc'
  
  IF (mype_filter == 0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- Read initial state from file ', infile
  s = 1
  stat(s) = NF_OPEN(infile, NF_NOWRITE, fileid)
  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) THEN
        WRITE(*, *) 'NetCDF error in opening initialization file, no.', i
        STOP
     END IF
  END DO

  ! Read size of state vector
  s = 1
  stat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_p_read)

  ! Read rank stored in file
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, rank)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from file, no.', i
  END DO

  checkdim: IF (dim_p_read == dim_p_cov .AND. rank >= dim_ens-1) THEN

     IF (mype_filter == 0) WRITE (*,'(8x,a)') '--- Read covariance matrix'

     ! Inquire IDs for mean state, singular vectors and values
     s = 1
     stat(s) = NF_INQ_VARID(fileid, 'running_meanstate', id_state)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'V', id_eof)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

     ! EOF and singular values
     startv(2) = 1
     countv(2) = dim_ens-1
     startv(1) = 1
     countv(1) = dim_p_read
     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_eof, startv, countv, eof_p)

     s = s + 1
     stat(s) = NF_GET_VARA_DOUBLE(fileid, id_svals, 1, dim_ens-1, svals)

     s = s + 1
     stat(s) = nf_close(fileid)

     DO i = 1,  s
        IF (stat(i) /= NF_NOERR) &
             WRITE(*, *) 'NetCDF error in reading initialization file, no.', i
     END DO

     IF (mype_filter==0) THEN
        WRITE(*,*) 'svals', svals
     END IF

! ************************************************
! *** Generate ensemble of perturbations       ***
! ************************************************

     IF (dim_ens>1) THEN
        ! Only initialize Omega if ensemble size > 0

        IF (mype_filter==0) THEN

           WRITE (*,'(a,8x,a)') 'FESOM-PDAF','--- generate state ensemble'

           ! *** Generate uniform orthogonal matrix OMEGA ***
           CALL PDAF_seik_omega(dim_ens-1, Omega, 1, 1)

            ! ***      Generate ensemble of states         ***
            ! *** x_i = x + sqrt(FAC) eofV (Omega C^(-1))t ***

            ! A = Omega C^(-1)
           DO col = 1, dim_ens-1
              DO row = 1, dim_ens
                 Omega(row, col) = Omega(row,col) * svals(col)
              END DO
           END DO
        END IF
        CALL MPI_Bcast(Omega, dim_ens*(dim_ens-1), MPI_DOUBLE_PRECISION, 0, &
             COMM_filter, MPIerr)
     END IF

     ! *** state_ens = state + sqrt(dim_ens-1) eofV A^T ***
     
     allocate(ens_p_per(dim_p_cov,dim_ens))
     ens_p_per   = 0.0

     IF (dim_ens>1) THEN
        ! Only add perturbations if ensemble size > 0

        fac = varscale * SQRT(REAL(dim_ens-1))
        
        ! =========          =====             =====        =========
        ! ens_p_per := fac * eof_p * transpose(Omega) + 0 * ens_p_per

        CALL DGEMM('n', 't', dim_p_cov, dim_ens, dim_ens-1, &
             fac, eof_p, dim_p_cov, Omega, dim_ens, 0.0, ens_p_per, dim_p_cov) ! matrix operation

     END IF
   
   ! Ensemble mean state from all initial model fields (state_p)
   CALL collect_state_PDAF(dim_p, state_p)
   DO col = 1,dim_ens
      ens_p(1:dim_p,col) = state_p(1:dim_p)
   END DO
   
   ! Add perturbation to to-be-perturbed part of full state vector:
   DO j = 1,nfields_per
      ! field offset in covariance matrix
      o1 = offsets_cov(j)
      o2 = offsets_cov(j) + dim_fields(idxs(j))
      ! field offset in full state vector
      f1 = offset(idxs(j))
      f2 = offset(idxs(j)) + dim_fields(idxs(j))
      ens_p(f1:f2,:) = ens_p(f1:f2,:) + ens_p_per(o1:o2,:)
      
!~       if (mype_filter==0) then
!~       WRITE (*,'(1X,A15,1X,A5,1X,I7,1X,I7,1X,A15,1X,I7,1X,I7)') 'Index in covar: ', &
!~                                                                  trim(sfields(idxs(j))%variable), &
!~                                                                  o1, o2, &
!~                                                                  'Index in state: ', &
!~                                                                  f1, f2
!~       endif
      
   ENDDO
   
   DEALLOCATE(ens_p_per)

   ! *** Treshold values ***
   treshold: DO col= 1,dim_ens
      
      ! surface fields
      n_treshold_ssh_p = 0
      fields_surf: Do i= 1,myDim_nod2D
          ! SSH: set to +/- 1.7m where larger than that
          ! (note: control simulation has minimum of -1.34 in Jan 2016; and -1.67 to 1.56 full-year max.)
          IF ( ens_p(i+offset(id% ssh),col) < -1.7 ) n_treshold_ssh_p = n_treshold_ssh_p + 1 ! counter
          IF ( ens_p(i+offset(id% ssh),col) > +1.7 ) n_treshold_ssh_p = n_treshold_ssh_p + 1 ! counter
          ens_p(i+offset(id% ssh),col)=min(max(ens_p(i+offset(id% ssh),col),-1.7),1.7)       ! apply correction
      END DO fields_surf
      
      ! write out correction counters for SSH and sea ice:
      CALL MPI_Allreduce(n_treshold_ssh_p, n_treshold_ssh_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      IF (mype_filter == 0) WRITE(*,*) 'FESOM-PDAF init_ens_pdaf - SSH  - Number of field corrections: ', n_treshold_ssh_g, ' on member ', col
      
      ! 3D fields
      n_treshold_temp_p = 0
      n_treshold_salt_p = 0
      fields_3D: DO i = 1, (mesh_fesom%nl-1) * myDim_nod2D 
          ! temp: set to -1.895 Celsius where smaller than that
          IF ( ens_p(i+offset(id% temp),col) < -1.895D0) n_treshold_temp_p = n_treshold_temp_p + 1 ! counter
          ens_p(i+offset(id% temp),col)= max(ens_p(i+offset(id% temp),col),-1.895D0)               ! apply correction
          ! salt: set to null where negative
          IF ( ens_p(i+offset(id% salt),col) < 0.0) n_treshold_salt_p = n_treshold_salt_p + 1      ! counter
          ens_p(i+offset(id% salt),col)= max(ens_p(i+offset(id% salt),col),0.)                     ! apply correction
          ! BGC: set to "tiny" where negative
          ens_p(i+ offset(id% O2     ),col) = max(tiny     ,ens_p(i+ offset(id% O2     ),col))
          ens_p(i+ offset(id% DIN    ),col) = max(tiny*1e-3,ens_p(i+ offset(id% DIN    ),col))
          ens_p(i+ offset(id% DIC    ),col) = max(tiny*1e-3,ens_p(i+ offset(id% DIC    ),col))
          ens_p(i+ offset(id% Alk    ),col) = max(tiny*1e-3,ens_p(i+ offset(id% Alk    ),col))
          
      END DO fields_3D
      
      ! write out correction counters for temperature and salinity:
      CALL MPI_Allreduce(n_treshold_temp_p, n_treshold_temp_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      CALL MPI_Allreduce(n_treshold_salt_p, n_treshold_salt_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      IF (mype_filter == 0) WRITE(*,*) 'FESOM-PDAF init_ens_pdaf - Temp - Number of field corrections: ', n_treshold_temp_g, ' on member ', col
      IF (mype_filter == 0) WRITE(*,*) 'FESOM-PDAF init_ens_pdaf - Salt - Number of field corrections: ', n_treshold_salt_g, ' on member ', col

   END DO treshold

  ELSE checkdim

      ! *** Rank stored in file is smaller than requested EOF rank ***
     WRITE(*,*) 'FESOM-PDAF: ','ERROR: Rank stored in file is smaller than requested EOF rank ...'
     WRITE(*,*) 'FESOM-PDAF: ','... or dim_p not equal--------> init_ens_pdaf'
     CALL abort_parallel()

  END IF checkdim

! ****************
! *** clean up ***
! ****************

  DEALLOCATE(svals, eof_p, omega, idxs, offsets_cov)
  ENDIF ! this_is_pdaf_restart

END SUBROUTINE init_ens_pdaf
  
