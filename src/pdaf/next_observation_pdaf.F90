! *****************************
! -----------------------------
! *** next_observation_pdaf ***
! -----------------------------
SUBROUTINE next_observation_pdaf(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The subroutine is called before each forecast phase
! by PDAF\_get\_state. It has to initialize the number 
! of time steps until the next available observation 
! (nsteps) and the current model time (time). In 
! addition the exit flag (exit) has to be initialized.
! It indicates if the data assimilation process is 
! completed such that the ensemble loop in the model 
! routine can be exited.
!
! The routine is called by all filter processes. 
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_model, task_id
  USE mod_assim_pdaf, ONLY: delt_obs_ocn, step_null, assim_time
  USE recom_config, ONLY: secondsperday
  USE g_clock, &
     ONLY: timenew, daynew

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step in PDAF
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  IF (stepnow==step_null) THEN
      ! at start, one assimilation step at first step of day 1
      nsteps=1
      ! at start, last step of day 1 for double assim.-calls:
      ! nsteps=2*delt_obs_ocn-1
      ! at start, last step of day 1 for single assim.-calls:
      ! nsteps=delt_obs_ocn
      assim_time = INT( REAL(nsteps) / REAL(delt_obs_ocn) * REAL(secondsperday))
      
      IF (mype_model==0 .AND. task_id==1) &
      WRITE (*,'(a,1x,a,i8,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
            'FESOM-PDAF','next_observation_pdaf; at step_null, nsteps = ', nsteps, &
            'stepnow = ', stepnow, 'Day', daynew, 'Time', FLOOR(timenew/3600.0),'h',INT(MOD(timenew,3600.0)/60.0),'min'
      
  ELSE
      ! daily assimilation steps during model time loop
      nsteps=delt_obs_ocn
      
      IF (mype_model==0 .AND. task_id==1) &
      WRITE (*,'(a,1x,a,i8,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
            'FESOM-PDAF','next_observation_pdaf; nsteps = ', nsteps, &
            'stepnow = ', stepnow, 'Day', daynew, 'Time', FLOOR(timenew/3600.0),'h',INT(MOD(timenew,3600.0)/60.0),'min'
  ENDIF


! *********************************
! *** Set current physical time ***
! *********************************

  time = 1.0

! *********************
! *** Set exit flag ***
! *********************

  doexit = 0

END SUBROUTINE next_observation_pdaf


! ***************************************
! ---------------------------------------
! *** next_observation_pdaf_ncalls2_1 ***
! ---------------------------------------
SUBROUTINE next_observation_pdaf_ncalls2_1(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! Alternative routine for repreated calls to PDAFomi_assimilate_XXX
! This routine is for the first call.
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_model, task_id
  USE mod_assim_pdaf, ONLY: delt_obs_ocn, step_null, assim_time
  USE recom_config, ONLY: secondsperday
  USE g_clock, &
     ONLY: timenew, daynew

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step in PDAF
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  ! first of two calls: assimilate now and at the next call
  nsteps=1
  
  IF (mype_model==0 .AND. task_id==1) &
  WRITE (*,'(a,1x,a,i8,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
        'FESOM-PDAF','next_observation_pdaf_ncalls2_1; nsteps = ', nsteps, &
        'stepnow = ', stepnow, 'Day', daynew, 'Time', FLOOR(timenew/3600.0),'h',INT(MOD(timenew,3600.0)/60.0),'min'

! *********************************
! *** Set current physical time ***
! *********************************

  time = 1.0

! *********************
! *** Set exit flag ***
! *********************

  doexit = 0

END SUBROUTINE next_observation_pdaf_ncalls2_1


! ***************************************
! ---------------------------------------
! *** next_observation_pdaf_ncalls2_2 ***
! ---------------------------------------
SUBROUTINE next_observation_pdaf_ncalls2_2(stepnow, nsteps, doexit, time)

! !DESCRIPTION:
! Alternative routine for repreated calls to PDAFomi_assimilate_XXX
! This routine is for the first call.
!
! !USES:
  USE mod_parallel_pdaf, ONLY: mype_model, task_id
  USE mod_assim_pdaf, ONLY: delt_obs_ocn, step_null, assim_time
  USE recom_config, ONLY: secondsperday
  USE g_clock, &
     ONLY: timenew, daynew

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in)  :: stepnow  ! Number of the current time step in PDAF
  INTEGER, INTENT(out) :: nsteps   ! Number of time steps until next obs
  INTEGER, INTENT(out) :: doexit   ! Whether to exit forecasting (1 for exit)
  REAL, INTENT(out)    :: time     ! Current model (physical) time

! !CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_next_obs)


! *************************************************************
! *** Determine number of time steps until next observation ***
! *************************************************************

  ! second of two calls: assimilate now and after one day
  nsteps=2*delt_obs_ocn-1
  
  IF (mype_model==0 .AND. task_id==1) &
  WRITE (*,'(a,1x,a,i8,1x,a,1x,i5,1x,a,1x,i3,1x,a,1x,i2,a,1x,i2,a)') &
        'FESOM-PDAF','next_observation_pdaf_ncalls2_2; nsteps = ', nsteps, &
        'stepnow = ', stepnow, 'Day', daynew, 'Time', FLOOR(timenew/3600.0),'h',INT(MOD(timenew,3600.0)/60.0),'min'

! *********************************
! *** Set current physical time ***
! *********************************

  time = 1.0

! *********************
! *** Set exit flag ***
! *********************

  doexit = 0

END SUBROUTINE next_observation_pdaf_ncalls2_2
