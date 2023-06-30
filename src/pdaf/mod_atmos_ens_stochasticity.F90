MODULE mod_atmos_ens_stochasticity

! Description: Adds stochastic synoptic variability to the atmospheric
! forcing fields for an ensemble of atmospheric forcings.

! Routines in this module:
! --- init_atmos_ens_stochasticity()
! --- add_atmos_ens_stochasticity(istep)
! --- init_atmos_stochasticity_output
! --- write_atmos_stochasticity_output(istep)


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
  
  CHARACTER(len=200) :: fname_atm ! filename to write out atmospheric stochasticity

LOGICAL :: disturb_xwind
LOGICAL :: disturb_ywind
LOGICAL :: disturb_humi
LOGICAL :: disturb_qlw
LOGICAL :: disturb_qsr
LOGICAL :: disturb_tair
LOGICAL :: disturb_prec
LOGICAL :: disturb_snow
LOGICAL :: disturb_mslp

LOGICAL :: atmos_stochasticity_ON

REAL :: varscale_wind
REAL :: varscale_tair

LOGICAL :: write_atmos_st = .false.


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

IF (mype_world==0) THEN
WRITE(*,*) 'Reading atm. covariance data from netCDF ', trim(filename)
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

ALLOCATE(perturbation_xwind (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_ywind (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_humi  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_qlw   (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_qsr   (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_tair  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_prec  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_snow  (myDim_nod2D + eDim_nod2D))
ALLOCATE(perturbation_mslp  (myDim_nod2D + eDim_nod2D))

IF (this_is_pdaf_restart) THEN

CALL read_atmos_stochasticity_restart()

ELSE

perturbation_xwind = 0.0
perturbation_ywind = 0.0
perturbation_humi  = 0.0
perturbation_qlw   = 0.0
perturbation_qsr   = 0.0
perturbation_tair  = 0.0
perturbation_prec  = 0.0
perturbation_snow  = 0.0
perturbation_mslp  = 0.0

ENDIF


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
perturbation_xwind ( :myDim_nod2D) = (1-arc) * perturbation_xwind ( :myDim_nod2D) + arc * varscale_wind * perturbation(atm_offset(id_atm% xwind) +1 : atm_offset(id_atm% xwind) +myDim_nod2D)
perturbation_ywind ( :myDim_nod2D) = (1-arc) * perturbation_ywind ( :myDim_nod2D) + arc * varscale_wind * perturbation(atm_offset(id_atm% ywind) +1 : atm_offset(id_atm% ywind) +myDim_nod2D)
perturbation_humi  ( :myDim_nod2D) = (1-arc) * perturbation_humi  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% humi ) +1 : atm_offset(id_atm% humi ) +myDim_nod2D)
perturbation_qlw   ( :myDim_nod2D) = (1-arc) * perturbation_qlw   ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% qlw  ) +1 : atm_offset(id_atm% qlw  ) +myDim_nod2D)
perturbation_qsr   ( :myDim_nod2D) = (1-arc) * perturbation_qsr   ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% qsr  ) +1 : atm_offset(id_atm% qsr  ) +myDim_nod2D)
perturbation_tair  ( :myDim_nod2D) = (1-arc) * perturbation_tair  ( :myDim_nod2D) + arc * varscale_tair * perturbation(atm_offset(id_atm% tair ) +1 : atm_offset(id_atm% tair ) +myDim_nod2D)
perturbation_prec  ( :myDim_nod2D) = (1-arc) * perturbation_prec  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% prec ) +1 : atm_offset(id_atm% prec ) +myDim_nod2D)
perturbation_snow  ( :myDim_nod2D) = (1-arc) * perturbation_snow  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% snow ) +1 : atm_offset(id_atm% snow ) +myDim_nod2D)
perturbation_mslp  ( :myDim_nod2D) = (1-arc) * perturbation_mslp  ( :myDim_nod2D) + arc * perturbation(atm_offset(id_atm% mslp ) +1 : atm_offset(id_atm% mslp ) +myDim_nod2D)

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

! add perturbation to atmospheric fields:
IF (disturb_xwind) atmdata(i_xwind,:) = atmdata(i_xwind,:) + perturbation_xwind
IF (disturb_ywind) atmdata(i_ywind,:) = atmdata(i_ywind,:) + perturbation_ywind
IF (disturb_humi)  atmdata(i_humi ,:) = atmdata(i_humi ,:) + perturbation_humi 
IF (disturb_qlw)   atmdata(i_qlw  ,:) = atmdata(i_qlw  ,:) + perturbation_qlw  
IF (disturb_tair)  atmdata(i_tair ,:) = atmdata(i_tair ,:) + perturbation_tair 
IF (disturb_prec)  atmdata(i_prec ,:) = atmdata(i_prec ,:) + perturbation_prec 
IF (disturb_snow)  atmdata(i_snow ,:) = atmdata(i_snow ,:) + perturbation_snow 
IF (disturb_mslp)  atmdata(i_mslp ,:) = atmdata(i_mslp ,:) + perturbation_mslp

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
    INTEGER :: dimID_n2D ! dimension ID: nodes
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
    INTEGER :: varID_restart_xwind
    INTEGER :: varID_restart_ywind
    INTEGER :: varID_restart_humi 
    INTEGER :: varID_restart_qlw  
    INTEGER :: varID_restart_qsr  
    INTEGER :: varID_restart_tair 
    INTEGER :: varID_restart_prec 
    INTEGER :: varID_restart_snow 
    INTEGER :: varID_restart_mslp
    INTEGER :: dimarray(2)
    
    
IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Initialize netCDF file to protocol atmospheric stochasticity'
END IF
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'.nc'

s = 1
stat(s) = NF_CREATE(TRIM(fname_atm),0,fileid)
s = s+1

! --- define dimensions:
stat(s) = NF_DEF_DIM(fileid,'myDim_nod2D', myDim_nod2D, dimID_n2D)
s = s+1
stat(s) = NF_DEF_DIM(fileid,'step', NF_UNLIMITED, dimId_step)
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
    INTEGER :: dimID_n2D ! dimension ID: nodes
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
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'.nc'

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
    INTEGER :: varID_restart_xwind
    INTEGER :: varID_restart_ywind
    INTEGER :: varID_restart_humi 
    INTEGER :: varID_restart_qlw  
    INTEGER :: varID_restart_qsr  
    INTEGER :: varID_restart_tair 
    INTEGER :: varID_restart_prec 
    INTEGER :: varID_restart_snow 
    INTEGER :: varID_restart_mslp
    
IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Write atmospheric stochasticity restart to netCDF at the end.'
END IF
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'.nc'

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
    INTEGER :: varID_restart_xwind
    INTEGER :: varID_restart_ywind
    INTEGER :: varID_restart_humi 
    INTEGER :: varID_restart_qlw  
    INTEGER :: varID_restart_qsr  
    INTEGER :: varID_restart_tair 
    INTEGER :: varID_restart_prec 
    INTEGER :: varID_restart_snow 
    INTEGER :: varID_restart_mslp
    
IF (mype_world==0) THEN
WRITE (*, '(/a, 1x, a)') 'FESOM-PDAF', 'Read atmospheric stochasticity at restart from netCDF.'
END IF
    
! --- open file:
fname_atm = TRIM(DAoutput_path)//'atmos_'//mype_string//'.nc'

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
