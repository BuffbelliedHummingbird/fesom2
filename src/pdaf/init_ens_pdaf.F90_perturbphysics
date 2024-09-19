! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022-02 - Frauke B     - Adapted for FESOM2.1

! !USES:
  USE mod_assim_pdaf, &
       ONLY: file_init, path_init, read_inistate, file_inistate, varscale, &
       offset, ASIM_START_USE_CLIM_STATE, this_is_pdaf_restart, mesh_fesom, &
       id, dim_fields, offset
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, COMM_filter, abort_parallel
  USE g_PARSUP, &
       ONLY: MPIerr, edim_nod2d, mydim_nod2d
  USE o_arrays, &
       ONLY: eta_n, tr_arr, uv, wvel
  USE i_arrays, &
       ONLY: a_ice
  USE g_ic3d

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
  INTEGER :: dim_p_file                ! Local state dimension in init file
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
  CHARACTER(len=150) :: file           ! File holding initial state estimate
  INTEGER :: dim_p_read
  LOGICAL :: runningmean               ! True: Initialize state vector from
                                       ! nc-file running mean
  REAL, ALLOCATABLE :: ens_p_phy(:,:)  ! Ensemble states (physics part)
  INTEGER :: dim_p_phy                 ! Local state dimension (physics part)
  INTEGER :: idx1,idx2                 ! Indeces
  INTEGER, ALLOCATABLE :: idxs(:)
  
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
  ELSE

! **********************
! *** INITIALIZATION ***
! **********************

  dim_p_phy =   dim_fields(id% ssh   ) + &
                dim_fields(id% u     ) + &
                dim_fields(id% v     ) + &
                dim_fields(id% w     ) + &
                dim_fields(id% temp  ) + &
                dim_fields(id% salt  ) + &
                dim_fields(id% a_ice )
  
  dim_p_read = dim_p_phy
  
  write(mype_string,'(i4.4)') mype_filter 

  mype0: IF (mype_filter == 0) THEN
    WRITE (*, '(/a, 8x,a)') 'FESOM-PDAF', 'Generate state ensemble from covariance matrix'
    WRITE (*, '(a, 8x,a)') &
         'FESOM-PDAF', '--- use 2nd order exact sampling (SEIK type)'
    WRITE (*, '(a, 8x,a,i5)') 'FESOM-PDAF', '--- number of EOFs:',dim_ens-1
  END IF mype0

  ! allocate memory for temporary fields
  ALLOCATE(eof_p(dim_p_phy, dim_ens-1))
  ALLOCATE(svals(dim_ens-1))
  ALLOCATE(omega(dim_ens, dim_ens-1))

! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************
 
  file=Trim(path_init)//Trim(file_init)//TRIM(mype_string)//'.nc'
  
  IF (mype_filter == 0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- Read initial state from file ', file

  s = 1
  stat(s) = NF_OPEN(file, NF_NOWRITE, fileid)

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
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_p_file)

  ! Read rank stored in file
  s = s + 1
  stat(s) = NF_INQ_DIMID(fileid, 'rank', id_dim)
  s = s + 1
  stat(s) = NF_INQ_DIMLEN(fileid, id_dim, rank)

  DO i = 1,  s
     IF (stat(i) /= NF_NOERR) &
          WRITE(*, *) 'NetCDF error in reading dimensions from file, no.', i
  END DO

  checkdim: IF (dim_p_read == dim_p_file .AND. rank >= dim_ens-1) THEN

     IF (mype_filter == 0) WRITE (*,'(8x,a)') '--- Read covariance matrix'

     ! Inquire IDs for mean state, singular vectors and values
     s = 1
     stat(s) = NF_INQ_VARID(fileid, 'running_meanstate', id_state)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'V', id_eof)
     s = s + 1
     stat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

     ! *** Read initialization information ***
     
     ! Ensemble mean state vector (state_p)
     runningmean = .False.
     IF (runningmean) THEN
        s = s + 1
        stat(s) = NF_GET_VAR_DOUBLE(fileid, id_state, state_p)
     ELSE
        CALL collect_state_PDAF(dim_p, state_p)
     END IF

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

! **********************************************************
! *** Set variance of diagnostic fields to zero          ***
! **********************************************************

! Is already zero (see gen_cov tool).
! eof_p(offset-first-biogeo:end,1:dim_ens-1)=0

! *****************************************
! *** Generate ensemble of model states ***
! *****************************************

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
     
     DO col = 1,dim_ens
        ens_p(1:dim_p,col) = state_p(1:dim_p)
     END DO
     
     allocate(ens_p_phy(dim_p_phy,dim_ens))
     ens_p_phy   = 0.0

     IF (dim_ens>1) THEN
        ! Only add perturbations if ensemble size > 0

        fac = varscale * SQRT(REAL(dim_ens-1))
        
        ! =========          =====             =====        =========
        ! ens_p_phy := fac * eof_p * transpose(Omega) + 0 * ens_p_phy

        CALL DGEMM('n', 't', dim_p_phy, dim_ens, dim_ens-1, &
             fac, eof_p, dim_p_phy, Omega, dim_ens, 0.0, ens_p_phy, dim_p_phy) ! matrix operation
     END IF
     
   ! Add perturbation to physics part:
   ALLOCATE(idxs(6))
   idxs = (/ id%SSH, id%u, id%v, id%temp, id%salt, id%a_ice /)
   
   DO j = 1,6
      idx1 = offset(idxs(j))
      idx2 = offset(idxs(j)) + dim_fields(idxs(j))
      ens_p(idx1:idx2,:) = ens_p(idx1:idx2,:) + ens_p_phy(idx1:idx2,:)
   ENDDO

   ! *** Treshold values ***
   treshold: DO col= 1,dim_ens
      
      ! surface fields
      n_treshold_ssh_p = 0
      n_treshold_sic_p = 0
      Do i= 1,myDim_nod2D
          ! SSH: set to +/- 1.7m where larger than that
          ! (note: control simulation has minimum of -1.34 in Jan 2016; and -1.67 to 1.56 full-year max.)
          IF ( ens_p(i+offset(id% ssh),col) < -1.7 ) n_treshold_ssh_p = n_treshold_ssh_p + 1 ! counter
          IF ( ens_p(i+offset(id% ssh),col) > +1.7 ) n_treshold_ssh_p = n_treshold_ssh_p + 1 ! counter
          ens_p(i+offset(id% ssh),col)=min(max(ens_p(i+offset(id% ssh),col),-1.7),1.7)       ! apply correction
          ! SIC: set to null where negative
          !      set to 1 where larger than that
          IF ( ens_p(i+offset(id% a_ice),col) < 0.0 ) n_treshold_sic_p = n_treshold_sic_p + 1 ! counter
          IF ( ens_p(i+offset(id% a_ice),col) > 1.0 ) n_treshold_sic_p = n_treshold_sic_p + 1 ! counter
          ens_p(i+offset(id% a_ice),col)= max(ens_p(i+offset(id% a_ice),col),0.)              ! apply correction
          ens_p(i+offset(id% a_ice),col)= min(ens_p(i+offset(id% a_ice),col),1.)              ! apply correction
      END DO
      
      ! write out correction counters for SSH and sea ice:
      CALL MPI_Allreduce(n_treshold_ssh_p, n_treshold_ssh_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      CALL MPI_Allreduce(n_treshold_sic_p, n_treshold_sic_g, 1, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
      IF (mype_filter == 0) WRITE(*,*) 'FESOM-PDAF init_ens_pdaf - SSH  - Number of field corrections: ', n_treshold_ssh_g, ' on member ', col
      IF (mype_filter == 0) WRITE(*,*) 'FESOM-PDAF init_ens_pdaf - SIC  - Number of field corrections: ', n_treshold_sic_g, ' on member ', col
      
      ! 3D fields
      n_treshold_temp_p = 0
      n_treshold_salt_p = 0
      DO i = 1, (mesh_fesom%nl-1) * myDim_nod2D 
          ! temp: set to -1.895 Celsius where smaller than that
          IF ( ens_p(i+offset(id% temp),col) < -1.895D0) n_treshold_temp_p = n_treshold_temp_p + 1 ! counter
          ens_p(i+offset(id% temp),col)= max(ens_p(i+offset(id% temp),col),-1.895D0)               ! apply correction
          ! salt: set to null where negative
          IF ( ens_p(i+offset(id% salt),col) < 0.0) n_treshold_salt_p = n_treshold_salt_p + 1      ! counter
          ens_p(i+offset(id% salt),col)= max(ens_p(i+offset(id% salt),col),0.)                     ! apply correction
      END DO
      
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

  DEALLOCATE(svals, eof_p, omega)
  ENDIF ! this_is_pdaf_restart

END SUBROUTINE init_ens_pdaf
  
