MODULE mod_atmos_ens_stochasticity

! Description: Adds stochastic synoptic variability to the atmospheric
! forcing fields for an ensemble of atmospheric forcings.

! Routines in this module:
! --- init_atmos_ens_stochasticity()
! --- add_atmos_ens_stochasticity(istep)
! --- init_atmos_stochasticity_output
! --- write_atmos_stochasticity_output(istep)

! Covariance file contains statistical information on all 9 atmospheric
! forcing fields. Which of these fields shall be perturbed is set
! in namelist.fesom.pdaf.


  USE mod_parallel_pdaf, &
       ONLY: mype_filter, mype_model, mype_world, &
             COMM_filter, filterpe, task_id, COMM_model
  USE mod_assim_pdaf, &
       ! dimensions:
       ONLY: dim_ens, dim_state_p, &
       ! netCDF file:
       path_atm_cov, forget
  USE g_clock, &
       ONLY: cyearnew, cyearold
  USE g_PARSUP, &
       ONLY: myDim_nod2D, MPI_DOUBLE_PRECISION, MPIerr, eDim_nod2D
  USE g_sbf, &
	ONLY: atmdata, &
	      i_xwind, i_ywind, i_humi, &
	      i_qsr, i_qlw, i_tair, i_prec, i_mslp, i_snow
  USE g_comm_auto
  USE g_config, &
    ONLY: step_per_day
       
  IMPLICIT NONE
  
  INCLUDE 'netcdf.inc'

  INTEGER                 :: cnt                     ! counters
  REAL(8),ALLOCATABLE, save  :: eof_p_cvrfile(:,:)      ! Matrix of eigenvectors of covariance matrix (all fields read from file)
  REAL(8),ALLOCATABLE, save  :: eof_p(:,:)              ! Matrix of eigenvectors of covariance matrix (only activated fields)
  REAL(8),ALLOCATABLE, save  :: svals(:)                ! Singular values
  REAL(8),ALLOCATABLE        :: omega(:,:)           ! Transformation matrix Omega
  REAL(8),ALLOCATABLE        :: omega_v(:)           ! Transformation vector for local ensemble member
  INTEGER                 :: rank                    ! Rank stored in cov.-file
  CHARACTER(len=4)        :: mype_string             ! String for process rank
  CHARACTER(len=110)      :: filename                ! Name of covariance netCDF file
  INTEGER,          save  :: nfields_cvrfile         ! Number of atmospheric forcing fields in covariance netCDF file
  INTEGER,          save  :: nfields                 ! Number of activated atmospheric forcing fields
  REAL,ALLOCATABLE        :: perturbation(:)         ! Vector containing perturbation field for local ensemble member
  
  REAL,ALLOCATABLE, save  :: perturbation_humi (:)   ! Final perturbations for each variable
  REAL,ALLOCATABLE, save  :: perturbation_prec (:)
  REAL,ALLOCATABLE, save  :: perturbation_snow (:)
  REAL,ALLOCATABLE, save  :: perturbation_mslp (:)
  REAL,ALLOCATABLE, save  :: perturbation_qlw  (:)
  REAL,ALLOCATABLE, save  :: perturbation_qsr  (:)
  REAL,ALLOCATABLE, save  :: perturbation_tair (:)
  REAL,ALLOCATABLE, save  :: perturbation_xwind(:)
  REAL,ALLOCATABLE, save  :: perturbation_ywind(:)
  
  REAL,ALLOCATABLE, save  :: atmdata_debug(:,:)
  
  TYPE field_ids
     INTEGER :: humi
     INTEGER :: prec
     INTEGER :: snow
     INTEGER :: mslp
     INTEGER :: qlw 
     INTEGER :: qsr
     INTEGER :: tair
     INTEGER :: xwind
     INTEGER :: ywind
  END TYPE field_ids
  
  ! Type variable holding field IDs in atmospheric state vector
  TYPE(field_ids)    , save :: id_atm               ! field IDs of perturbed fields
  TYPE(field_ids)    , save :: id_cvrf              ! field IDs of all fields in covariance file
  
  INTEGER,ALLOCATABLE, save :: atm_offset(:)        ! offset of perturbed fields in atmospheric state vector
  INTEGER,ALLOCATABLE, save :: atm_offset_cvrf(:)   ! offset of fields hold in covariance fields
  
  CHARACTER(len=200) :: fname_atm  ! filename to write out atmospheric stochasticity

LOGICAL :: disturb_xwind           ! which atmospheric fields to be perturbed
LOGICAL :: disturb_ywind           ! (set in namelist)
LOGICAL :: disturb_humi
LOGICAL :: disturb_qlw
LOGICAL :: disturb_qsr
LOGICAL :: disturb_tair
LOGICAL :: disturb_prec
LOGICAL :: disturb_snow
LOGICAL :: disturb_mslp

LOGICAL :: atmos_stochasticity_ON   ! if any atmospheric fields to be perturbed

REAL :: varscale_wind = 0.1         ! scaling factors
REAL :: varscale_tair = 1.0
REAL :: varscale_humi = 0.5
REAL :: varscale_qlw  = 0.5

LOGICAL :: write_atmos_st = .false. ! wether to protocol the perturbed atmospheric fields,
                                    ! i.e. writing output at every time step

REAL :: stable_rmse = 0 ! (ocean temperature) ensemble spread after 16 months of assimilation 

CONTAINS

! ************************************
! ************************************
! *** init_atmos_ens_stochasticity ***
! ************************************
! ************************************

SUBROUTINE init_atmos_ens_stochasticity()

IMPLICIT NONE

! Local variables:
INTEGER :: s, i                          ! Counters
INTEGER :: ncstat(50)                    ! Status flag for netCDF commands
INTEGER :: fileid                        ! netCDF file handle
INTEGER :: id_dim, &
           id_eof, id_svals              ! handles for netCDF commands
INTEGER :: startv(2), countv(2)          ! specifier for netCDF commands
INTEGER :: dim_p_file                    ! state dimension read from cov.-file

! **********************
! *** INITIALIZATION ***
! **********************

IF (mype_world==0) THEN
WRITE(*,*) 'FESOM-PDAF: Init atmospheric perturbation.'
END IF

! Count number of activated fields that shall be perturbed:
nfields = COUNT((/ disturb_humi,  &
                   disturb_prec,  &
                   disturb_snow,  &
                   disturb_mslp,  &
                   disturb_qlw,   &
                   disturb_qsr,   &
                   disturb_tair,  &
                   disturb_xwind, &
                   disturb_ywind /))

!~ ! Debugging:                   
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'NFIELDS', nfields
                   
ALLOCATE(eof_p(nfields * myDim_nod2D, dim_ens-1))

! Nine fields in covariance file:
nfields_cvrfile = 9
ALLOCATE(eof_p_cvrfile(nfields_cvrfile * myDim_nod2D, dim_ens-1))

ALLOCATE(svals(dim_ens-1))

! Composition of atmospheric state vector as in covariance file:
id_cvrf% humi = 1
id_cvrf% prec = 2
id_cvrf% snow = 3
id_cvrf% mslp = 4
id_cvrf% qlw  = 5
id_cvrf% qsr  = 6
id_cvrf% tair = 7
id_cvrf% xwind= 8
id_cvrf% ywind= 9

! Composition of atmospheric state vector, only activated fields:
cnt = 1
IF (disturb_xwind) THEN; id_atm% xwind = cnt; cnt=cnt+1; ENDIF
IF (disturb_ywind) THEN; id_atm% ywind = cnt; cnt=cnt+1; ENDIF
IF (disturb_humi ) THEN; id_atm% humi  = cnt; cnt=cnt+1; ENDIF
IF (disturb_qlw  ) THEN; id_atm% qlw   = cnt; cnt=cnt+1; ENDIF
IF (disturb_qsr  ) THEN; id_atm% qsr   = cnt; cnt=cnt+1; ENDIF
IF (disturb_tair ) THEN; id_atm% tair  = cnt; cnt=cnt+1; ENDIF
IF (disturb_prec ) THEN; id_atm% prec  = cnt; cnt=cnt+1; ENDIF
IF (disturb_snow ) THEN; id_atm% snow  = cnt; cnt=cnt+1; ENDIF
IF (disturb_mslp ) THEN; id_atm% mslp  = cnt; cnt=cnt+1; ENDIF

! Offset of each field in state vector spaced equally (all fields are surface fields):
! As in covariance file:
ALLOCATE(atm_offset_cvrf(nfields_cvrfile))
DO i=1,nfields_cvrfile
	atm_offset_cvrf(i)=(i-1)*myDim_nod2D +1
END DO
! Only activated fields:
ALLOCATE(atm_offset(nfields))
DO i=1,nfields
	atm_offset(i)=(i-1)*myDim_nod2D +1
END DO

write(mype_string,'(i4.4)') mype_model

! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

filename = TRIM(path_atm_cov)//'cov_'//TRIM(mype_string)//'.nc'

s = 1
ncstat(s) = NF_OPEN(trim(filename), NF_NOWRITE, fileid)

IF (mype_world==0) THEN
WRITE(*,*) 'FESOM-PDAF: Reading atm. covariance data from netCDF ', trim(filename)
END IF


IF (ncstat(1) /= NF_NOERR) THEN
   WRITE(*, *) 'FESOM-PDAF: NetCDF error in opening atm. covariance file, no.', i
   STOP
END IF


! Read size of state vector
s = 1
ncstat(s) = NF_INQ_DIMID(fileid, 'dim_state', id_dim)
s = s + 1
ncstat(s) = NF_INQ_DIMLEN(fileid, id_dim, dim_p_file)

! Read rank stored in file
s = s + 1
ncstat(s) = NF_INQ_DIMID(fileid, 'rank', id_dim)
s = s + 1
ncstat(s) = NF_INQ_DIMLEN(fileid, id_dim, rank)

DO i = 1,  s
     IF (ncstat(i) /= NF_NOERR) THEN
          WRITE(*, *) 'FESOM-PDAF: NetCDF error in reading dimensions from atm. covariance file, no.', i
          STOP
     ENDIF
END DO

checkdim: IF (dim_p_file /= (nfields_cvrfile*myDim_nod2D)) THEN
     WRITE(*,*) 'FESOM-PDAF: Dimensions inconsistent reading covariance of atmospheric forcing.'
     STOP
ENDIF checkdim

! Inquire IDs for mean state, singular vectors and values
s = 1
ncstat(s) = NF_INQ_VARID(fileid, 'V', id_eof)
s = s + 1
ncstat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

! EOF and singular values
startv(2) = 1
countv(2) = dim_ens-1
startv(1) = 1
countv(1) = nfields_cvrfile * myDim_nod2D
s = s + 1
ncstat(s) = NF_GET_VARA_DOUBLE(fileid, id_eof, startv, countv, eof_p_cvrfile)

s = s + 1
ncstat(s) = NF_GET_VARA_DOUBLE(fileid, id_svals, 1, dim_ens-1, svals)

s = s + 1
ncstat(s) = nf_close(fileid)

DO i = 1,  s
        IF (ncstat(i) /= NF_NOERR) THEN
             WRITE(*, *) 'FESOM-PDAF: NetCDF error in reading atm. covariance file, no.', i
             STOP
        ENDIF
END DO
        
! EOF for activated fields only:
IF(disturb_xwind) eof_p (atm_offset(id_atm% xwind) : atm_offset(id_atm% xwind)+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% xwind) : atm_offset_cvrf(id_cvrf% xwind)+myDim_nod2D, :)
IF(disturb_ywind) eof_p (atm_offset(id_atm% ywind) : atm_offset(id_atm% ywind)+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% ywind) : atm_offset_cvrf(id_cvrf% ywind)+myDim_nod2D, :)
IF(disturb_humi ) eof_p (atm_offset(id_atm% humi ) : atm_offset(id_atm% humi )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% humi ) : atm_offset_cvrf(id_cvrf% humi )+myDim_nod2D, :)
IF(disturb_qlw  ) eof_p (atm_offset(id_atm% qlw  ) : atm_offset(id_atm% qlw  )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% qlw  ) : atm_offset_cvrf(id_cvrf% qlw  )+myDim_nod2D, :)
IF(disturb_qsr  ) eof_p (atm_offset(id_atm% qsr  ) : atm_offset(id_atm% qsr  )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% qsr  ) : atm_offset_cvrf(id_cvrf% qsr  )+myDim_nod2D, :)
IF(disturb_tair ) eof_p (atm_offset(id_atm% tair ) : atm_offset(id_atm% tair )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% tair ) : atm_offset_cvrf(id_cvrf% tair )+myDim_nod2D, :)
IF(disturb_prec ) eof_p (atm_offset(id_atm% prec ) : atm_offset(id_atm% prec )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% prec ) : atm_offset_cvrf(id_cvrf% prec )+myDim_nod2D, :)
IF(disturb_snow ) eof_p (atm_offset(id_atm% snow ) : atm_offset(id_atm% snow )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% snow ) : atm_offset_cvrf(id_cvrf% snow )+myDim_nod2D, :)
IF(disturb_mslp ) eof_p (atm_offset(id_atm% mslp ) : atm_offset(id_atm% mslp )+myDim_nod2D, :)  =  eof_p_cvrfile (atm_offset_cvrf(id_cvrf% mslp ) : atm_offset_cvrf(id_cvrf% mslp )+myDim_nod2D, :)

!~ ! Debugging output:
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'xwind', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% xwind) : atm_offset_cvrf(id_cvrf% xwind)+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'ywind', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% ywind) : atm_offset_cvrf(id_cvrf% ywind)+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'humi ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% humi ) : atm_offset_cvrf(id_cvrf% humi )+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'qlw  ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% qlw  ) : atm_offset_cvrf(id_cvrf% qlw  )+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'qsr  ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% qsr  ) : atm_offset_cvrf(id_cvrf% qsr  )+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'tair ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% tair ) : atm_offset_cvrf(id_cvrf% tair )+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'prec ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% prec ) : atm_offset_cvrf(id_cvrf% prec )+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'snow ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% snow ) : atm_offset_cvrf(id_cvrf% snow )+2, :)
!~ IF (mype_world==0) write(*,*) 'disturb_atmos_debug ', 'eof_p_cvr ', 'mslp ', eof_p_cvrfile (atm_offset_cvrf(id_cvrf% mslp ) : atm_offset_cvrf(id_cvrf% mslp )+2, :)

!~ IF ((mype_world==0) .and. (disturb_xwind )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'xwind', eof_p (atm_offset(id_atm% xwind) : atm_offset(id_atm% xwind)+2, :)
!~ IF ((mype_world==0) .and. (disturb_ywind )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'ywind', eof_p (atm_offset(id_atm% ywind) : atm_offset(id_atm% ywind)+2, :)
!~ IF ((mype_world==0) .and. (disturb_humi  )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'humi ', eof_p (atm_offset(id_atm% humi ) : atm_offset(id_atm% humi )+2, :)
!~ IF ((mype_world==0) .and. (disturb_qlw   )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'qlw  ', eof_p (atm_offset(id_atm% qlw  ) : atm_offset(id_atm% qlw  )+2, :)
!~ IF ((mype_world==0) .and. (disturb_qsr   )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'qsr  ', eof_p (atm_offset(id_atm% qsr  ) : atm_offset(id_atm% qsr  )+2, :)
!~ IF ((mype_world==0) .and. (disturb_tair  )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'tair ', eof_p (atm_offset(id_atm% tair ) : atm_offset(id_atm% tair )+2, :)
!~ IF ((mype_world==0) .and. (disturb_prec  )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'prec ', eof_p (atm_offset(id_atm% prec ) : atm_offset(id_atm% prec )+2, :)
!~ IF ((mype_world==0) .and. (disturb_snow  )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'snow ', eof_p (atm_offset(id_atm% snow ) : atm_offset(id_atm% snow )+2, :)
!~ IF ((mype_world==0) .and. (disturb_mslp  )) write(*,*) 'disturb_atmos_debug ', 'eof_p ', 'mslp ', eof_p (atm_offset(id_atm% mslp ) : atm_offset(id_atm% mslp )+2, :)

deallocate(eof_p_cvrfile)

END SUBROUTINE




! ***********************************
! ***********************************
! *** add_atmos_ens_stochasticity ***
! ***********************************
! ***********************************

SUBROUTINE add_atmos_ens_stochasticity(istep)

USE mod_assim_pdaf, &
    ONLY: this_is_pdaf_restart

IMPLICIT NONE

! Arguments:
INTEGER, INTENT(in)    :: istep

! Local variables:
INTEGER :: row, col                      ! counters
REAL :: fac                              ! Square-root of dim_ens or dim_ens-1
REAL :: arc, varscale                    ! autoregression coefficient and scaling factor
CHARACTER(len=3) :: istep_string

ALLOCATE(omega(dim_ens, dim_ens-1))
ALLOCATE(omega_v(dim_ens-1))
ALLOCATE(perturbation(nfields * myDim_nod2D))

! set parameters:
varscale = 10 ! 10: blowup at early time-step (3)
arc      = 1/REAL(step_per_day)

! ****************************************
! *** Generate ensemble of atm. states ***
! ****************************************

IF (mype_model==0) THEN

   IF (istep==0) WRITE (*,'(a,8x,a)') 'FESOM-PDAF','generate atm. state ensemble'

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
   
   Omega_v = Omega(task_id,:)
   
END IF

! Clean up:
DEALLOCATE(omega)

! rank 0 within the model communicator (i.e. model rank 0)
CALL MPI_Bcast(Omega_v, dim_ens-1, MPI_DOUBLE_PRECISION, 0, &
	 COMM_model, MPIerr)

fac = varscale * SQRT(REAL(dim_ens-1)) ! varscale: scaling factor for ensemble variance
perturbation = 0.0

!    ____           =====   _______    ____
!    pert = fac * ( eof_p * omega_v) + null

CALL DGEMV('n', nfields*myDim_nod2D, dim_ens-1, fac, eof_p, nfields*myDim_nod2D, omega_v, 1, 0, perturbation, 1) ! dgemv: matrix-vector multiplication

IF (istep==1) THEN

	IF(disturb_xwind) ALLOCATE(perturbation_xwind (myDim_nod2D + eDim_nod2D))
	IF(disturb_ywind) ALLOCATE(perturbation_ywind (myDim_nod2D + eDim_nod2D))
	IF(disturb_humi ) ALLOCATE(perturbation_humi  (myDim_nod2D + eDim_nod2D))
	IF(disturb_qlw  ) ALLOCATE(perturbation_qlw   (myDim_nod2D + eDim_nod2D))
	IF(disturb_qsr  ) ALLOCATE(perturbation_qsr   (myDim_nod2D + eDim_nod2D))
	IF(disturb_tair ) ALLOCATE(perturbation_tair  (myDim_nod2D + eDim_nod2D))
	IF(disturb_prec ) ALLOCATE(perturbation_prec  (myDim_nod2D + eDim_nod2D))
	IF(disturb_snow ) ALLOCATE(perturbation_snow  (myDim_nod2D + eDim_nod2D))
	IF(disturb_mslp ) ALLOCATE(perturbation_mslp  (myDim_nod2D + eDim_nod2D))

	IF (this_is_pdaf_restart) THEN

		CALL read_atmos_stochasticity_restart()

	ELSE

		IF(disturb_xwind) perturbation_xwind = 0.0
		IF(disturb_ywind) perturbation_ywind = 0.0
		IF(disturb_humi ) perturbation_humi  = 0.0
		IF(disturb_qlw  ) perturbation_qlw   = 0.0
		IF(disturb_qsr  ) perturbation_qsr   = 0.0
		IF(disturb_tair ) perturbation_tair  = 0.0
		IF(disturb_prec ) perturbation_prec  = 0.0
		IF(disturb_snow ) perturbation_snow  = 0.0
		IF(disturb_mslp ) perturbation_mslp  = 0.0

ENDIF


!~ ! debugging output:
!~ ALLOCATE(atmdata_debug(nfields,myDim_nod2D))
!~ IF(disturb_xwind) atmdata_debug(id_atm% xwind,:) = atmdata(i_xwind,:myDim_nod2D)
!~ IF(disturb_ywind) atmdata_debug(id_atm% ywind,:) = atmdata(i_ywind,:myDim_nod2D)
!~ IF(disturb_humi ) atmdata_debug(id_atm% humi ,:) = atmdata(i_humi ,:myDim_nod2D)
!~ IF(disturb_qlw  ) atmdata_debug(id_atm% qlw  ,:) = atmdata(i_qlw  ,:myDim_nod2D)
!~ IF(disturb_qsr  ) atmdata_debug(id_atm% qsr  ,:) = atmdata(i_qsr  ,:myDim_nod2D)
!~ IF(disturb_tair ) atmdata_debug(id_atm% tair ,:) = atmdata(i_tair ,:myDim_nod2D)
!~ IF(disturb_prec ) atmdata_debug(id_atm% prec ,:) = atmdata(i_prec ,:myDim_nod2D)
!~ IF(disturb_snow ) atmdata_debug(id_atm% snow ,:) = atmdata(i_snow ,:myDim_nod2D)
!~ IF(disturb_mslp ) atmdata_debug(id_atm% mslp ,:) = atmdata(i_mslp ,:myDim_nod2D)


END IF

! autoregressive: next perturbation from last perturbation and new stochastic element
! perturb[n+1] = (1 - arc) * perturb[n] + arc * varscale * random
IF(disturb_xwind) perturbation_xwind ( :myDim_nod2D) = (1-arc) * perturbation_xwind ( :myDim_nod2D) + arc * varscale_wind * perturbation(atm_offset(id_atm% xwind) : atm_offset(id_atm% xwind) +myDim_nod2D)
IF(disturb_ywind) perturbation_ywind ( :myDim_nod2D) = (1-arc) * perturbation_ywind ( :myDim_nod2D) + arc * varscale_wind * perturbation(atm_offset(id_atm% ywind) : atm_offset(id_atm% ywind) +myDim_nod2D)
IF(disturb_humi ) perturbation_humi  ( :myDim_nod2D) = (1-arc) * perturbation_humi  ( :myDim_nod2D) + arc * varscale_humi * perturbation(atm_offset(id_atm% humi ) : atm_offset(id_atm% humi ) +myDim_nod2D)
IF(disturb_qlw  ) perturbation_qlw   ( :myDim_nod2D) = (1-arc) * perturbation_qlw   ( :myDim_nod2D) + arc * varscale_qlw  * perturbation(atm_offset(id_atm% qlw  ) : atm_offset(id_atm% qlw  ) +myDim_nod2D)
IF(disturb_qsr  ) perturbation_qsr   ( :myDim_nod2D) = (1-arc) * perturbation_qsr   ( :myDim_nod2D) + arc                 * perturbation(atm_offset(id_atm% qsr  ) : atm_offset(id_atm% qsr  ) +myDim_nod2D)
IF(disturb_tair ) perturbation_tair  ( :myDim_nod2D) = (1-arc) * perturbation_tair  ( :myDim_nod2D) + arc * varscale_tair * perturbation(atm_offset(id_atm% tair ) : atm_offset(id_atm% tair ) +myDim_nod2D)
IF(disturb_prec ) perturbation_prec  ( :myDim_nod2D) = (1-arc) * perturbation_prec  ( :myDim_nod2D) + arc                 * perturbation(atm_offset(id_atm% prec ) : atm_offset(id_atm% prec ) +myDim_nod2D)
IF(disturb_snow ) perturbation_snow  ( :myDim_nod2D) = (1-arc) * perturbation_snow  ( :myDim_nod2D) + arc                 * perturbation(atm_offset(id_atm% snow ) : atm_offset(id_atm% snow ) +myDim_nod2D)
IF(disturb_mslp ) perturbation_mslp  ( :myDim_nod2D) = (1-arc) * perturbation_mslp  ( :myDim_nod2D) + arc                 * perturbation(atm_offset(id_atm% mslp ) : atm_offset(id_atm% mslp ) +myDim_nod2D)

! fill external nodes:
IF(disturb_xwind) CALL exchange_nod( perturbation_xwind)
IF(disturb_ywind) CALL exchange_nod( perturbation_ywind)
IF(disturb_humi ) CALL exchange_nod( perturbation_humi)
IF(disturb_qlw  ) CALL exchange_nod( perturbation_qlw)
IF(disturb_qsr  ) CALL exchange_nod( perturbation_qsr)
IF(disturb_tair ) CALL exchange_nod( perturbation_tair)
IF(disturb_prec ) CALL exchange_nod( perturbation_prec)
IF(disturb_snow ) CALL exchange_nod( perturbation_snow)
IF(disturb_mslp ) CALL exchange_nod( perturbation_mslp)

!~ ! debugging output:
!~ IF ((mype_world==0) .and. (disturb_xwind )) write(*,*) 'disturb_atmos_debug ', 'perturbation_xwind', perturbation_xwind(:2)
!~ IF ((mype_world==0) .and. (disturb_ywind )) write(*,*) 'disturb_atmos_debug ', 'perturbation_ywind', perturbation_ywind(:2)
!~ IF ((mype_world==0) .and. (disturb_humi  )) write(*,*) 'disturb_atmos_debug ', 'perturbation_humi ', perturbation_humi (:2)
!~ IF ((mype_world==0) .and. (disturb_qlw   )) write(*,*) 'disturb_atmos_debug ', 'perturbation_qlw  ', perturbation_qlw  (:2)
!~ IF ((mype_world==0) .and. (disturb_qsr   )) write(*,*) 'disturb_atmos_debug ', 'perturbation_qsr  ', perturbation_qsr  (:2)
!~ IF ((mype_world==0) .and. (disturb_tair  )) write(*,*) 'disturb_atmos_debug ', 'perturbation_tair ', perturbation_tair (:2)
!~ IF ((mype_world==0) .and. (disturb_prec  )) write(*,*) 'disturb_atmos_debug ', 'perturbation_prec ', perturbation_prec (:2)
!~ IF ((mype_world==0) .and. (disturb_snow  )) write(*,*) 'disturb_atmos_debug ', 'perturbation_snow ', perturbation_snow (:2)
!~ IF ((mype_world==0) .and. (disturb_mslp  )) write(*,*) 'disturb_atmos_debug ', 'perturbation_mslp ', perturbation_mslp (:2)


! add perturbation to atmospheric fields:
IF (disturb_xwind) atmdata(i_xwind,:) = atmdata(i_xwind,:) + perturbation_xwind
IF (disturb_ywind) atmdata(i_ywind,:) = atmdata(i_ywind,:) + perturbation_ywind
IF (disturb_humi)  atmdata(i_humi ,:) = atmdata(i_humi ,:) + perturbation_humi 
IF (disturb_qlw)   atmdata(i_qlw  ,:) = atmdata(i_qlw  ,:) + perturbation_qlw  
IF (disturb_tair)  atmdata(i_tair ,:) = atmdata(i_tair ,:) + perturbation_tair 
IF (disturb_prec)  atmdata(i_prec ,:) = atmdata(i_prec ,:) + perturbation_prec 
IF (disturb_snow)  atmdata(i_snow ,:) = atmdata(i_snow ,:) + perturbation_snow 
IF (disturb_mslp)  atmdata(i_mslp ,:) = atmdata(i_mslp ,:) + perturbation_mslp
IF (disturb_qsr)   atmdata(i_qsr  ,:) = atmdata(i_qsr  ,:) + perturbation_qsr  

! corrections:
! night: shortwave is zero.
IF(disturb_qsr) THEN
  WHERE(atmdata(i_qsr,:) /= 0)
  atmdata(i_qsr  ,:) = atmdata(i_qsr  ,:) + perturbation_qsr
  ENDWHERE
ENDIF
! rain, snow, humidity, downwelling shortwave and longwave radiation \
! must not be negative:
IF(disturb_prec) THEN
  WHERE(atmdata(i_prec,:) <0 )
  atmdata(i_prec,:)=0
  ENDWHERE
ENDIF
IF(disturb_snow) THEN
  WHERE(atmdata(i_snow,:) <0 )
  atmdata(i_snow,:)=0
  ENDWHERE
ENDIF
IF(disturb_humi) THEN
  WHERE(atmdata(i_humi,:) <0 )
  atmdata(i_humi,:)=0
  ENDWHERE
  WHERE(atmdata(i_humi,:) >1 )
  atmdata(i_humi,:)=1
  ENDWHERE
ENDIF
IF(disturb_qlw) THEN
  WHERE(atmdata(i_qlw,:) <0 )
  atmdata(i_qlw,:)=0
  ENDWHERE
ENDIF
IF(disturb_qsr) THEN
  WHERE(atmdata(i_qsr,:) <0 )
  atmdata(i_qsr,:)=0
  ENDWHERE
ENDIF

!~ ! debugging output:
!~ IF ((mype_world==0) .and. (disturb_xwind )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_xwind ,:2)', atmdata(i_xwind ,:2)
!~ IF ((mype_world==0) .and. (disturb_ywind )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_ywind ,:2)', atmdata(i_ywind ,:2)
!~ IF ((mype_world==0) .and. (disturb_humi  )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_humi  ,:2)', atmdata(i_humi  ,:2)
!~ IF ((mype_world==0) .and. (disturb_qlw   )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_qlw   ,:2)', atmdata(i_qlw   ,:2)
!~ IF ((mype_world==0) .and. (disturb_qsr   )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_tair  ,:2)', atmdata(i_tair  ,:2)
!~ IF ((mype_world==0) .and. (disturb_tair  )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_prec  ,:2)', atmdata(i_prec  ,:2)
!~ IF ((mype_world==0) .and. (disturb_prec  )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_snow  ,:2)', atmdata(i_snow  ,:2)
!~ IF ((mype_world==0) .and. (disturb_snow  )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_mslp  ,:2)', atmdata(i_mslp  ,:2)
!~ IF ((mype_world==0) .and. (disturb_mslp  )) write(*,*) 'disturb_atmos_debug ', 'atmdata(i_qsr   ,:2)', atmdata(i_qsr   ,:2)


DEALLOCATE(perturbation,omega_v)

END SUBROUTINE



! ***************************************
! ***************************************
! *** init_atmos_stochasticity_output ***
! ***************************************
! ***************************************

SUBROUTINE init_atmos_stochasticity_output

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_PARSUP, &
         ONLY: myDim_nod2D
    USE mod_assim_pdaf, &
         ONLY: DAoutput_path
         
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
    
    ! Local variables:
    INTEGER :: s ! auxiliary status counter
    INTEGER :: i ! counter
    INTEGER :: fileid ! ID of netCDF file
    INTEGER :: stat(500)                    ! auxiliary: status array
    INTEGER :: dimID_n2D   ! dimension ID: nodes
    INTEGER :: dimID_step  ! dimension ID: time steps
    INTEGER :: dimID_dummy ! dummy dimension
    INTEGER :: varID_xwind
    INTEGER :: varID_ywind
    INTEGER :: varID_humi 
    INTEGER :: varID_qlw  
    INTEGER :: varID_qsr  
    INTEGER :: varID_tair 
    INTEGER :: varID_prec 
    INTEGER :: varID_snow 
    INTEGER :: varID_mslp
    INTEGER :: varID_restart_xwind
    INTEGER :: varID_restart_ywind
    INTEGER :: varID_restart_humi 
    INTEGER :: varID_restart_qlw  
    INTEGER :: varID_restart_qsr  
    INTEGER :: varID_restart_tair 
    INTEGER :: varID_restart_prec 
    INTEGER :: varID_restart_snow 
    INTEGER :: varID_restart_mslp
    INTEGER :: varID_restart_rmse, varID_restart_forget
    INTEGER :: dimarray(2)
    
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'_'//cyearnew//'.nc'

IF (mype_world==0) THEN
WRITE(*,*) 'FESOM-PDAF: cyearold, cyearnew', cyearold, cyearnew
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize netCDF file to protocol atmospheric stochasticity:', fname_atm
END IF

s = 1
stat(s) = NF_CREATE(TRIM(fname_atm),0,fileid)
s = s+1

! --- define dimensions:
stat(s) = NF_DEF_DIM(fileid,'myDim_nod2D', myDim_nod2D, dimID_n2D)
s = s+1
stat(s) = NF_DEF_DIM(fileid,'step', NF_UNLIMITED, dimId_step)
s = s+1
stat(s) = NF_DEF_DIM(fileid,'dimensionless', 1, dimID_dummy)
s = s+1

DO i = 1,  s - 1
     IF (stat(i) /= NF_NOERR) &
     WRITE(*, *) 'NetCDF error in defining dimensions in atmos. netCDF file, no.', i
END DO

! --- define variables:
dimarray(1) = dimID_n2D
dimarray(2) = dimID_step
s = 1

! perturbed atmospheric forcing fields
IF (disturb_xwind) stat(s) = NF_DEF_VAR(fileid, 'xwind', NF_FLOAT, 2, dimarray(1:2), varID_xwind)
IF (disturb_xwind) s = s+1                                           
IF (disturb_ywind) stat(s) = NF_DEF_VAR(fileid, 'ywind', NF_FLOAT, 2, dimarray(1:2), varID_ywind)
IF (disturb_ywind) s = s+1                                             
IF (disturb_humi ) stat(s) = NF_DEF_VAR(fileid, 'humi' , NF_FLOAT, 2, dimarray(1:2), varID_humi )
IF (disturb_humi ) s = s+1                                             
IF (disturb_qlw  ) stat(s) = NF_DEF_VAR(fileid, 'qlw'  , NF_FLOAT, 2, dimarray(1:2), varID_qlw  )
IF (disturb_qlw  ) s = s+1                                             
IF (disturb_qsr  ) stat(s) = NF_DEF_VAR(fileid, 'qsr'  , NF_FLOAT, 2, dimarray(1:2), varID_qsr  )
IF (disturb_qsr  ) s = s+1                                             
IF (disturb_tair ) stat(s) = NF_DEF_VAR(fileid, 'tair' , NF_FLOAT, 2, dimarray(1:2), varID_tair )
IF (disturb_tair ) s = s+1                                             
IF (disturb_prec ) stat(s) = NF_DEF_VAR(fileid, 'prec' , NF_FLOAT, 2, dimarray(1:2), varID_prec )
IF (disturb_prec ) s = s+1                                             
IF (disturb_snow ) stat(s) = NF_DEF_VAR(fileid, 'snow' , NF_FLOAT, 2, dimarray(1:2), varID_snow )
IF (disturb_snow ) s = s+1                                             
IF (disturb_mslp ) stat(s) = NF_DEF_VAR(fileid, 'mslp' , NF_FLOAT, 2, dimarray(1:2), varID_mslp )
IF (disturb_mslp ) s = s+1

! perturbation for restart:
IF (disturb_xwind) stat(s) = NF_DEF_VAR(fileid, 'restart_xwind', NF_FLOAT, 1, dimID_n2D, varID_restart_xwind)
IF (disturb_xwind) s = s+1
IF (disturb_ywind) stat(s) = NF_DEF_VAR(fileid, 'restart_ywind', NF_FLOAT, 1, dimID_n2D, varID_restart_ywind)
IF (disturb_ywind) s = s+1
IF (disturb_humi ) stat(s) = NF_DEF_VAR(fileid, 'restart_humi' , NF_FLOAT, 1, dimID_n2D, varID_restart_humi )
IF (disturb_humi ) s = s+1
IF (disturb_qlw  ) stat(s) = NF_DEF_VAR(fileid, 'restart_qlw'  , NF_FLOAT, 1, dimID_n2D, varID_restart_qlw  )
IF (disturb_qlw  ) s = s+1
IF (disturb_qsr  ) stat(s) = NF_DEF_VAR(fileid, 'restart_qsr'  , NF_FLOAT, 1, dimID_n2D, varID_restart_qsr  )
IF (disturb_qsr  ) s = s+1
IF (disturb_tair ) stat(s) = NF_DEF_VAR(fileid, 'restart_tair' , NF_FLOAT, 1, dimID_n2D, varID_restart_tair )
IF (disturb_tair ) s = s+1
IF (disturb_prec ) stat(s) = NF_DEF_VAR(fileid, 'restart_prec' , NF_FLOAT, 1, dimID_n2D, varID_restart_prec )
IF (disturb_prec ) s = s+1
IF (disturb_snow ) stat(s) = NF_DEF_VAR(fileid, 'restart_snow' , NF_FLOAT, 1, dimID_n2D, varID_restart_snow )
IF (disturb_snow ) s = s+1
IF (disturb_mslp ) stat(s) = NF_DEF_VAR(fileid, 'restart_mslp' , NF_FLOAT, 1, dimID_n2D, varID_restart_mslp )
IF (disturb_mslp ) s = s+1

! save stable ensemble spread (computed during month 16) and forgetting factor for restart:
stat(s) = NF_DEF_VAR(fileid, 'restart_rmse'   , NF_FLOAT, 1, dimID_dummy, varID_restart_rmse)
s = s+1
stat(s) = NF_DEF_VAR(fileid, 'restart_forget' , NF_FLOAT, 1, dimID_dummy, varID_restart_forget)
s = s+1

DO i = 1,  s - 1
     IF (stat(i) /= NF_NOERR) THEN
     WRITE(*,*) 'NetCDF error in defining variables in atmos. netCDF file, no.', i
     WRITE(*,*)  NF_STRERROR(stat(i))
     STOP
     END IF
END DO

stat(s) = NF_ENDDEF(fileid) 
s = s + 1
stat(1) = NF_CLOSE(fileid)

IF (stat(1) /= NF_NOERR) THEN
   WRITE(*, *) 'NetCDF error in closing atmos. netCDF file'
END IF

IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'netCDF file to protocol atmospheric stochasticity has been initialized'
END IF

END SUBROUTINE



! ****************************************
! ****************************************
! *** write_atmos_stochasticity_output ***
! ****************************************
! ****************************************

SUBROUTINE write_atmos_stochasticity_output(istep)

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_PARSUP, &
         ONLY: myDim_nod2D
    USE mod_assim_pdaf, &
         ONLY: DAoutput_path, step_null
         
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
    
    ! Arguments:
    INTEGER, INTENT(in)    :: istep
    
    ! Local variables:
    INTEGER :: s ! auxiliary status counter
    INTEGER :: i ! counter
    INTEGER :: fileid ! ID of netCDF file
    INTEGER :: stat(500)                    ! auxiliary: status array
    INTEGER :: dimID_n2D  ! dimension ID: nodes
    INTEGER :: dimID_step ! dimension ID: time steps
    INTEGER :: varID_xwind
    INTEGER :: varID_ywind
    INTEGER :: varID_humi 
    INTEGER :: varID_qlw  
    INTEGER :: varID_qsr  
    INTEGER :: varID_tair 
    INTEGER :: varID_prec 
    INTEGER :: varID_snow 
    INTEGER :: varID_mslp
    
    INTEGER :: posvec(2) ! write position in netCDF file
    INTEGER :: nmbvec(2) ! write dimension in netCDF file
    
    
IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Write atmospheric stochasticity to netCDF.'
END IF
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'_'//cyearnew//'.nc'

s=1
stat(s) = NF_OPEN(TRIM(fname_atm), NF_WRITE, fileid)
IF (stat(s) /= NF_NOERR) STOP 'error opening atmospheric stochasticity netCDF'

! ----- inquire variable IDs:

s=1
IF (disturb_xwind) stat(s) = NF_INQ_VARID( fileid, 'xwind', varID_xwind )
IF (disturb_xwind) s=s+1
IF (disturb_ywind) stat(s) = NF_INQ_VARID( fileid, 'ywind', varID_ywind )
IF (disturb_ywind) s=s+1
IF (disturb_humi ) stat(s) = NF_INQ_VARID( fileid, 'humi' , varID_humi  )
IF (disturb_humi ) s=s+1
IF (disturb_qlw  ) stat(s) = NF_INQ_VARID( fileid, 'qlw'  , varID_qlw   )
IF (disturb_qlw  ) s=s+1
IF (disturb_qsr  ) stat(s) = NF_INQ_VARID( fileid, 'qsr'  , varID_qsr   )
IF (disturb_qsr  ) s=s+1
IF (disturb_tair ) stat(s) = NF_INQ_VARID( fileid, 'tair' , varID_tair  )
IF (disturb_tair ) s=s+1
IF (disturb_prec ) stat(s) = NF_INQ_VARID( fileid, 'prec' , varID_prec  )
IF (disturb_prec ) s=s+1
IF (disturb_snow ) stat(s) = NF_INQ_VARID( fileid, 'snow' , varID_snow  )
IF (disturb_snow ) s=s+1
IF (disturb_mslp ) stat(s) = NF_INQ_VARID( fileid, 'mslp' , varID_mslp  )
IF (disturb_mslp ) s=s+1

DO i = 1, s - 1
  IF (stat(i) /= NF_NOERR) &
  WRITE(*, *) 'NetCDF error inquiring atmospheric stochasticity variable IDs, no.', i
END DO

! --- write variables:
posvec = (/ 1,           istep+step_null /)
nmbvec = (/ myDim_nod2D, 1               /)

s=1
IF (disturb_xwind) stat(s) = NF_PUT_VARA_REAL( fileid, varid_xwind, posvec, nmbvec, REAL(atmdata(i_xwind,:myDim_nod2D),4))
IF (disturb_xwind) s=s+1                                                                              
IF (disturb_ywind) stat(s) = NF_PUT_VARA_REAL( fileid, varid_ywind, posvec, nmbvec, REAL(atmdata(i_ywind,:myDim_nod2D),4))
IF (disturb_ywind) s=s+1                                                                              
IF (disturb_humi)  stat(s) = NF_PUT_VARA_REAL( fileid, varid_humi,  posvec, nmbvec, REAL(atmdata(i_humi ,:myDim_nod2D),4))
IF (disturb_humi)  s=s+1                                                                              
IF (disturb_qlw)   stat(s) = NF_PUT_VARA_REAL( fileid, varid_qlw,   posvec, nmbvec, REAL(atmdata(i_qlw  ,:myDim_nod2D),4))
IF (disturb_qlw)   s=s+1                                                                              
IF (disturb_qsr)   stat(s) = NF_PUT_VARA_REAL( fileid, varid_qsr,   posvec, nmbvec, REAL(atmdata(i_qsr  ,:myDim_nod2D),4))
IF (disturb_qsr)   s=s+1                                                                              
IF (disturb_tair)  stat(s) = NF_PUT_VARA_REAL( fileid, varid_tair,  posvec, nmbvec, REAL(atmdata(i_tair ,:myDim_nod2D),4))
IF (disturb_tair)  s=s+1                                                                              
IF (disturb_prec)  stat(s) = NF_PUT_VARA_REAL( fileid, varid_prec,  posvec, nmbvec, REAL(atmdata(i_prec ,:myDim_nod2D),4))
IF (disturb_prec)  s=s+1                                                                              
IF (disturb_snow)  stat(s) = NF_PUT_VARA_REAL( fileid, varid_snow,  posvec, nmbvec, REAL(atmdata(i_snow ,:myDim_nod2D),4))
IF (disturb_snow)  s=s+1                                                                              
IF (disturb_mslp)  stat(s) = NF_PUT_VARA_REAL( fileid, varid_mslp,  posvec, nmbvec, REAL(atmdata(i_mslp ,:myDim_nod2D),4))
IF (disturb_mslp)  s=s+1

DO i = 1, s - 1
  IF (stat(i) /= NF_NOERR) THEN
  WRITE(*, *) 'NetCDF error writing atmospheric stochasticity variable IDs, no.', i
  STOP
  END IF
END DO

stat(1) = NF_CLOSE(fileid)

  IF (stat(1) /= NF_NOERR) THEN
     WRITE(*, *) 'NetCDF error in closing NetCDF file'
  END IF

END SUBROUTINE


! *****************************************
! *****************************************
! *** write_atmos_stochasticity_restart ***
! *****************************************
! *****************************************

SUBROUTINE write_atmos_stochasticity_restart()

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_PARSUP, &
         ONLY: myDim_nod2D
    USE mod_assim_pdaf, &
         ONLY: DAoutput_path
         
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
    
    ! Local variables:
    INTEGER :: s ! auxiliary status counter
    INTEGER :: i ! counter
    INTEGER :: fileid ! ID of netCDF file
    INTEGER :: stat(500)  ! auxiliary: status array
    INTEGER :: dimID_n2D  ! dimension ID: nodes
    INTEGER :: dimID_step ! dimension ID: time steps
    INTEGER :: dimID_dummy ! dummy dimension
    INTEGER :: varID_restart_xwind
    INTEGER :: varID_restart_ywind
    INTEGER :: varID_restart_humi 
    INTEGER :: varID_restart_qlw  
    INTEGER :: varID_restart_qsr  
    INTEGER :: varID_restart_tair 
    INTEGER :: varID_restart_prec 
    INTEGER :: varID_restart_snow 
    INTEGER :: varID_restart_mslp
    INTEGER :: varID_restart_rmse, varID_restart_forget
    
IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Write atmospheric stochasticity restart to netCDF at the end.'
END IF
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'_'//cyearnew//'.nc'

s=1
stat(s) = NF_OPEN(TRIM(fname_atm), NF_WRITE, fileid)
IF (stat(s) /= NF_NOERR) STOP 'error opening atmospheric stochasticity netCDF at the end'

! ----- inquire variable IDs:

s=1
 IF (disturb_xwind) stat(s) = NF_INQ_VARID( fileid, 'restart_xwind', varID_restart_xwind )
 IF (disturb_xwind) s=s+1
 IF (disturb_ywind) stat(s) = NF_INQ_VARID( fileid, 'restart_ywind', varID_restart_ywind )
 IF (disturb_ywind) s=s+1
 IF (disturb_humi ) stat(s) = NF_INQ_VARID( fileid, 'restart_humi' , varID_restart_humi  )
 IF (disturb_humi ) s=s+1
 IF (disturb_qlw  ) stat(s) = NF_INQ_VARID( fileid, 'restart_qlw'  , varID_restart_qlw   )
 IF (disturb_qlw  ) s=s+1
 IF (disturb_qsr  ) stat(s) = NF_INQ_VARID( fileid, 'restart_qsr'  , varID_restart_qsr   )
 IF (disturb_qsr  ) s=s+1
 IF (disturb_tair ) stat(s) = NF_INQ_VARID( fileid, 'restart_tair' , varID_restart_tair  )
 IF (disturb_tair ) s=s+1
 IF (disturb_prec ) stat(s) = NF_INQ_VARID( fileid, 'restart_prec' , varID_restart_prec  )
 IF (disturb_prec ) s=s+1
 IF (disturb_snow ) stat(s) = NF_INQ_VARID( fileid, 'restart_snow' , varID_restart_snow  )
 IF (disturb_snow ) s=s+1
 IF (disturb_mslp ) stat(s) = NF_INQ_VARID( fileid, 'restart_mslp' , varID_restart_mslp  )
 IF (disturb_mslp ) s=s+1
 
 stat(s) = NF_INQ_VARID( fileid, 'restart_rmse' ,   varID_restart_rmse )
 s=s+1
 stat(s) = NF_INQ_VARID( fileid, 'restart_forget' , varID_restart_forget )
 s=s+1

DO i = 1, s - 1
  IF (stat(i) /= NF_NOERR) &
  WRITE(*, *) 'NetCDF error inquiring atmospheric stochasticity restart variable IDs at the end, no.', i
END DO

s=1
 IF (disturb_xwind) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_xwind, REAL(perturbation_xwind(:myDim_nod2D),4))
 IF (disturb_xwind) s=s+1
 IF (disturb_ywind) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_ywind, REAL(perturbation_ywind(:myDim_nod2D),4))
 IF (disturb_ywind) s=s+1
 IF (disturb_humi ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_humi,  REAL(perturbation_humi (:myDim_nod2D),4))
 IF (disturb_humi ) s=s+1
 IF (disturb_qlw  ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_qlw,   REAL(perturbation_qlw  (:myDim_nod2D),4))
 IF (disturb_qlw  ) s=s+1
 IF (disturb_qsr  ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_qsr,   REAL(perturbation_qsr  (:myDim_nod2D),4))
 IF (disturb_qsr  ) s=s+1
 IF (disturb_tair ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_tair,  REAL(perturbation_tair (:myDim_nod2D),4))
 IF (disturb_tair ) s=s+1
 IF (disturb_prec ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_prec,  REAL(perturbation_prec (:myDim_nod2D),4))
 IF (disturb_prec ) s=s+1
 IF (disturb_snow ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_snow,  REAL(perturbation_snow (:myDim_nod2D),4))
 IF (disturb_snow ) s=s+1
 IF (disturb_mslp ) stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_mslp,  REAL(perturbation_mslp (:myDim_nod2D),4))
 IF (disturb_mslp ) s=s+1
 
 stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_rmse,    REAL(stable_rmse,4))
 s=s+1
 stat(s) = NF_PUT_VAR_REAL( fileid, varid_restart_forget,  REAL(forget     ,4))
 s=s+1

DO i = 1, s - 1
  IF (stat(i) /= NF_NOERR) THEN
  WRITE(*, *) 'NetCDF error writing atmospheric stochasticity restart variable IDs at the end, no.', i
  STOP
  END IF
END DO

stat(1) = NF_CLOSE(fileid)

  IF (stat(1) /= NF_NOERR) THEN
     WRITE(*, *) 'NetCDF error in closing atmospheric stochasticity NetCDF file at the end'
  END IF

END SUBROUTINE

! ****************************************
! ****************************************
! *** read_atmos_stochasticity_restart ***
! ****************************************
! ****************************************

SUBROUTINE read_atmos_stochasticity_restart()

    USE g_config, &
         ONLY: runid, ResultPath
    USE g_PARSUP, &
         ONLY: myDim_nod2D
    USE mod_assim_pdaf, &
         ONLY: DAoutput_path
         
    IMPLICIT NONE

    INCLUDE 'netcdf.inc'
    
    ! Local variables:
    INTEGER :: s ! auxiliary status counter
    INTEGER :: i ! counter
    INTEGER :: fileid ! ID of netCDF file
    INTEGER :: stat(500)  ! auxiliary: status array
    INTEGER :: dimID_n2D  ! dimension ID: nodes
    INTEGER :: dimID_step ! dimension ID: time steps
    INTEGER :: dimID_dummy! dimension ID: dummy
    INTEGER :: varID_restart_xwind
    INTEGER :: varID_restart_ywind
    INTEGER :: varID_restart_humi 
    INTEGER :: varID_restart_qlw  
    INTEGER :: varID_restart_qsr  
    INTEGER :: varID_restart_tair 
    INTEGER :: varID_restart_prec 
    INTEGER :: varID_restart_snow 
    INTEGER :: varID_restart_mslp
    INTEGER :: varID_restart_rmse
    INTEGER :: varID_restart_forget
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'_'//cyearold//'.nc'

IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Read atmospheric stochasticity at restart from netCDF:', fname_atm
END IF

s=1
stat(s) = NF_OPEN(TRIM(fname_atm), NF_WRITE, fileid)
IF (stat(s) /= NF_NOERR) STOP 'error opening atmospheric stochasticity netCDF at restart'

! ----- inquire variable IDs:

s=1
IF (disturb_xwind) stat(s) = NF_INQ_VARID( fileid, 'restart_xwind', varID_restart_xwind )
IF (disturb_xwind) s=s+1
IF (disturb_ywind) stat(s) = NF_INQ_VARID( fileid, 'restart_ywind', varID_restart_ywind )
IF (disturb_ywind) s=s+1
IF (disturb_humi ) stat(s) = NF_INQ_VARID( fileid, 'restart_humi' , varID_restart_humi  )
IF (disturb_humi ) s=s+1
IF (disturb_qlw  ) stat(s) = NF_INQ_VARID( fileid, 'restart_qlw'  , varID_restart_qlw   )
IF (disturb_qlw  ) s=s+1
IF (disturb_qsr  ) stat(s) = NF_INQ_VARID( fileid, 'restart_qsr'  , varID_restart_qsr   )
IF (disturb_qsr  ) s=s+1
IF (disturb_tair ) stat(s) = NF_INQ_VARID( fileid, 'restart_tair' , varID_restart_tair  )
IF (disturb_tair ) s=s+1
IF (disturb_prec ) stat(s) = NF_INQ_VARID( fileid, 'restart_prec' , varID_restart_prec  )
IF (disturb_prec ) s=s+1
IF (disturb_snow ) stat(s) = NF_INQ_VARID( fileid, 'restart_snow' , varID_restart_snow  )
IF (disturb_snow ) s=s+1
IF (disturb_mslp ) stat(s) = NF_INQ_VARID( fileid, 'restart_mslp' , varID_restart_mslp  )
IF (disturb_mslp ) s=s+1

stat(s) = NF_INQ_VARID( fileid, 'restart_rmse' , varID_restart_rmse )
s=s+1
stat(s) = NF_INQ_VARID( fileid, 'restart_forget' , varID_restart_forget )
s=s+1

DO i = 1, s - 1
  IF (stat(i) /= NF_NOERR) &
  WRITE(*, *) 'NetCDF error inquiring atmospheric stochasticity restart variable IDs at restart, no.', i
END DO

s=1
IF (disturb_xwind) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_xwind, perturbation_xwind(:myDim_nod2D)) ! make use of automatic type conversion
IF (disturb_xwind) s=s+1                                                                                       ! reading real4 values from file into an array of doubles
IF (disturb_ywind) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_ywind, perturbation_ywind(:myDim_nod2D))
IF (disturb_ywind) s=s+1
IF (disturb_humi ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_humi,  perturbation_humi (:myDim_nod2D))
IF (disturb_humi ) s=s+1
IF (disturb_qlw  ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_qlw,   perturbation_qlw  (:myDim_nod2D))
IF (disturb_qlw  ) s=s+1
IF (disturb_qsr  ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_qsr,   perturbation_qsr  (:myDim_nod2D))
IF (disturb_qsr  ) s=s+1
IF (disturb_tair ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_tair,  perturbation_tair (:myDim_nod2D))
IF (disturb_tair ) s=s+1
IF (disturb_prec ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_prec,  perturbation_prec (:myDim_nod2D))
IF (disturb_prec ) s=s+1
IF (disturb_snow ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_snow,  perturbation_snow (:myDim_nod2D))
IF (disturb_snow ) s=s+1
IF (disturb_mslp ) stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_mslp,  perturbation_mslp (:myDim_nod2D))
IF (disturb_mslp ) s=s+1

stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_rmse,   stable_rmse)
s=s+1
stat(s) = NF_GET_VAR_DOUBLE( fileid, varid_restart_forget, forget)
s=s+1

DO i = 1, s - 1
  IF (stat(i) /= NF_NOERR) THEN
  WRITE(*, *) 'NetCDF error reading atmospheric stochasticity restart variables at restart, no.', i
  STOP
  END IF
END DO

stat(1) = NF_CLOSE(fileid)

  IF (stat(1) /= NF_NOERR) THEN
     WRITE(*, *) 'NetCDF error in closing atmospheric stochasticity NetCDF file at restart'
  END IF

END SUBROUTINE

END MODULE mod_atmos_ens_stochasticity
