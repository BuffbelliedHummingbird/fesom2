MODULE mod_atmos_ens_stochasticity

! Description: Adds stochastic synoptic variability to the atmospheric
! forcing fields for an ensemble of atmospheric forcings.

  USE mod_parallel_pdaf, &
       ONLY: mype_filter, mype_model, mype_world, &
             COMM_filter, filterpe, task_id, COMM_model
  USE mod_assim_pdaf, &
       ! dimensions:
       ONLY: dim_ens, dim_state_p, &
       ! netCDF file:
       path_atm_cov
       
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

  REAL,ALLOCATABLE, save  :: eof_p(:,:)              ! Matrix of eigenvectors of covariance matrix
  REAL,ALLOCATABLE, save  :: svals(:)                ! Singular values
  REAL(8),ALLOCATABLE        :: omega(:,:)           ! Transformation matrix Omega
  REAL(8),ALLOCATABLE        :: omega_v(:)           ! Transformation vector for local ensemble member
  INTEGER                 :: rank                    ! Rank stored in cov.-file
  CHARACTER(len=4)        :: mype_string             ! String for process rank
  CHARACTER(len=110)      :: filename                ! Name of covariance netCDF file
  INTEGER,          save  :: nfields                 ! Number of atmospheric forcing fields
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
  TYPE(field_ids)    , save :: id_atm
  INTEGER,ALLOCATABLE, save :: atm_offset(:)


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
WRITE(*,*) 'FESOM-PDAF: Init ensemble of atmospheric forcings from covariance matrix'
END IF

nfields = 9

ALLOCATE(eof_p(nfields * myDim_nod2D, dim_ens-1))
ALLOCATE(svals(dim_ens-1))

! Composition of atmospheric state vector (as in covariance file):
id_atm% humi = 1
id_atm% prec = 2
id_atm% snow = 3
id_atm% mslp = 4
id_atm% qlw  = 5
id_atm% qsr  = 6
id_atm% tair = 7
id_atm% xwind= 8
id_atm% ywind= 9

ALLOCATE(atm_offset(nfields))
DO i=1,nfields
	atm_offset(i)=(i-1)*myDim_nod2D
END DO

write(mype_string,'(i4.4)') mype_model

! *************************************************
! *** Initialize initial state and covar matrix ***
! *************************************************

filename = TRIM(path_atm_cov)//'cov_'//TRIM(mype_string)//'.nc'

s = 1
ncstat(s) = NF_OPEN(trim(filename), NF_NOWRITE, fileid)

DO i = 1,  s
     IF (ncstat(i) /= NF_NOERR) THEN
        WRITE(*, *) 'FESOM-PDAF: NetCDF error in opening atm. covariance file, no.', i
        STOP
     END IF
END DO

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

checkdim: IF (dim_p_file /= (nfields*myDim_nod2D)) THEN
     WRITE(*,*) 'FESOM-PDAF: Dimensions inconsistent reading covariance of atmospheric forcing.'
     STOP
ENDIF checkdim

IF (mype_world == 0) WRITE (*,'(8x,a)') 'FESOM-PDAF: Read atm. covariance matrix'

! Inquire IDs for mean state, singular vectors and values
s = 1
ncstat(s) = NF_INQ_VARID(fileid, 'V', id_eof)
s = s + 1
ncstat(s) = NF_INQ_VARID(fileid, 'sigma', id_svals)

! EOF and singular values
startv(2) = 1
countv(2) = dim_ens-1
startv(1) = 1
countv(1) = nfields * myDim_nod2D
s = s + 1
ncstat(s) = NF_GET_VARA_DOUBLE(fileid, id_eof, startv, countv, eof_p)

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
END SUBROUTINE




! ***********************************
! ***********************************
! *** add_atmos_ens_stochasticity ***
! ***********************************
! ***********************************

SUBROUTINE add_atmos_ens_stochasticity(istep) ! ens_p?

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

ALLOCATE(perturbation_xwind (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_ywind (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_humi  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_qlw   (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_qsr   (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_tair  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_prec  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_snow  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_mslp  (myDim_nod2D + eDim_nod2D))

perturbation_xwind = 0.0
perturbation_ywind = 0.0
perturbation_humi  = 0.0
perturbation_qlw   = 0.0
perturbation_qsr   = 0.0
perturbation_tair  = 0.0
perturbation_prec  = 0.0
perturbation_snow  = 0.0
perturbation_mslp  = 0.0

!~ ! debugging output:
!~ ALLOCATE(atmdata_debug(nfields,myDim_nod2D))
!~ atmdata_debug(id_atm% xwind,:) = atmdata(i_xwind,:myDim_nod2D)
!~ atmdata_debug(id_atm% ywind,:) = atmdata(i_ywind,:myDim_nod2D)
!~ atmdata_debug(id_atm% humi ,:) = atmdata(i_humi ,:myDim_nod2D)
!~ atmdata_debug(id_atm% qlw  ,:) = atmdata(i_qlw  ,:myDim_nod2D)
!~ atmdata_debug(id_atm% qsr  ,:) = atmdata(i_qsr  ,:myDim_nod2D)
!~ atmdata_debug(id_atm% tair ,:) = atmdata(i_tair ,:myDim_nod2D)
!~ atmdata_debug(id_atm% prec ,:) = atmdata(i_prec ,:myDim_nod2D)
!~ atmdata_debug(id_atm% snow ,:) = atmdata(i_snow ,:myDim_nod2D)
!~ atmdata_debug(id_atm% mslp ,:) = atmdata(i_mslp ,:myDim_nod2D)

END IF

! autoregressive: next perturbation from last perturbation and new stochastic element
perturbation_xwind ( :myDim_nod2D) = (1-arc) * perturbation_xwind ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% xwind) +1 : atm_offset(id_atm% xwind) +myDim_nod2D)
perturbation_ywind ( :myDim_nod2D) = (1-arc) * perturbation_ywind ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% ywind) +1 : atm_offset(id_atm% ywind) +myDim_nod2D)
perturbation_humi  ( :myDim_nod2D) = (1-arc) * perturbation_humi  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% humi ) +1 : atm_offset(id_atm% humi ) +myDim_nod2D)
perturbation_qlw   ( :myDim_nod2D) = (1-arc) * perturbation_qlw   ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% qlw  ) +1 : atm_offset(id_atm% qlw  ) +myDim_nod2D)
perturbation_qsr   ( :myDim_nod2D) = (1-arc) * perturbation_qsr   ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% qsr  ) +1 : atm_offset(id_atm% qsr  ) +myDim_nod2D)
perturbation_tair  ( :myDim_nod2D) = (1-arc) * perturbation_tair  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% tair ) +1 : atm_offset(id_atm% tair ) +myDim_nod2D)
perturbation_prec  ( :myDim_nod2D) = (1-arc) * perturbation_prec  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% prec ) +1 : atm_offset(id_atm% prec ) +myDim_nod2D)
perturbation_snow  ( :myDim_nod2D) = (1-arc) * perturbation_snow  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% snow ) +1 : atm_offset(id_atm% snow ) +myDim_nod2D)
perturbation_mslp  ( :myDim_nod2D) = (1-arc) * perturbation_mslp  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% mslp ) +1 : atm_offset(id_atm% mslp ) +myDim_nod2D)


!~ IF ((ANY( istep == (/ 1,33,65,97,129,161,193,225,257,289 /) ))) THEN
write(istep_string,'(i3.3)') istep
open (mype_world+1, file = 'atmdist_'//mype_string//'_'//istep_string//'.out')
write(mype_world+1,*) perturbation_xwind( :myDim_nod2D)
write(mype_world+1,*) perturbation_ywind( :myDim_nod2D)
write(mype_world+1,*) perturbation_humi ( :myDim_nod2D)
write(mype_world+1,*) perturbation_qlw  ( :myDim_nod2D)
write(mype_world+1,*) perturbation_qsr  ( :myDim_nod2D)
write(mype_world+1,*) perturbation_tair ( :myDim_nod2D)
write(mype_world+1,*) perturbation_prec ( :myDim_nod2D)
write(mype_world+1,*) perturbation_snow ( :myDim_nod2D)
write(mype_world+1,*) perturbation_mslp ( :myDim_nod2D)
close(mype_world+1)
!~ ENDIF

! fill external nodes:
CALL exchange_nod( perturbation_xwind)
CALL exchange_nod( perturbation_ywind)
CALL exchange_nod( perturbation_humi)
CALL exchange_nod( perturbation_qlw)
CALL exchange_nod( perturbation_qsr)
CALL exchange_nod( perturbation_tair)
CALL exchange_nod( perturbation_prec)
CALL exchange_nod( perturbation_snow)
CALL exchange_nod( perturbation_mslp)

! add perturabtion to atmospheric fields:
atmdata(i_xwind,:) = atmdata(i_xwind,:) + perturbation_xwind
atmdata(i_ywind,:) = atmdata(i_ywind,:) + perturbation_ywind
atmdata(i_humi ,:) = atmdata(i_humi ,:) + perturbation_humi 
atmdata(i_qlw  ,:) = atmdata(i_qlw  ,:) + perturbation_qlw  
atmdata(i_tair ,:) = atmdata(i_tair ,:) + perturbation_tair 
atmdata(i_prec ,:) = atmdata(i_prec ,:) + perturbation_prec 
atmdata(i_snow ,:) = atmdata(i_snow ,:) + perturbation_snow 
atmdata(i_mslp ,:) = atmdata(i_mslp ,:) + perturbation_mslp

! night: shortwave is zero.
WHERE(atmdata(i_qsr,:) /= 0)
atmdata(i_qsr  ,:) = atmdata(i_qsr  ,:) + perturbation_qsr
ENDWHERE
! rain, snow, humidity, downwelling shortwave and longwave radiation \
! must not be negative:
WHERE(atmdata(i_prec,:) <0 )
atmdata(i_prec,:)=0
ENDWHERE
WHERE(atmdata(i_snow,:) <0 )
atmdata(i_snow,:)=0
ENDWHERE
WHERE(atmdata(i_humi,:) <0 )
atmdata(i_humi,:)=0
ENDWHERE
WHERE(atmdata(i_humi,:) >1 )
atmdata(i_humi,:)=1
ENDWHERE
WHERE(atmdata(i_qlw,:) <0 )
atmdata(i_qlw,:)=0
ENDWHERE
WHERE(atmdata(i_qsr,:) <0 )
atmdata(i_qsr,:)=0
ENDWHERE

!~ xwind xwind xwind xwind
!~ ywind ywind ywind ywind
!~ humi  humi  humi  humi 
!~ qlw   qlw   qlw   qlw  
!~ qsr   qsr   qsr   qsr  
!~ tair  tair  tair  tair 
!~ prec  prec  prec  prec 
!~ snow  snow  snow  snow 
!~ mslp  mslp  mslp  mslp 

DEALLOCATE(perturbation,omega_v)

END SUBROUTINE

END MODULE mod_atmos_ens_stochasticity
