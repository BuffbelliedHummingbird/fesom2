! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE prepoststep_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This variant is used with the simplified interface of
! PDAF. In this case, the name of the routine is defined
! within PDAF. This routine just calls the prepoststep
! routine corresponding to the selected filter algorithm.
!
! !REVISION HISTORY:
! 2010-07 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022    - Frauke       - Adapted for FESOM2.0
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, npes_filter, COMM_filter, writepe, mype_world
  USE mod_assim_pdaf, & ! Variables for assimilation
       ONLY: step_null, filtertype, dim_lag, eff_dim_obs, loctype, &
       offset, proffiles_o, var_p, std_p, state_fcst, &
       monthly_state_f, monthly_state_a, write_monthly_mean,mon_snapshot_mem, &
       num_day_in_month, endday_of_month_in_year, startday_of_month_in_year, &
       depth_excl, depth_excl_no, this_is_pdaf_restart, mesh_fesom, &
       dim_fields, dim_fields_glob, offset, offset_glob, nfields, id, &
       timemean, monthly_state_m, delt_obs_ocn
  USE g_PARSUP, &
       ONLY: MPI_DOUBLE_PRECISION, MPI_SUM, MPIerr, MPI_STATUS_SIZE, &
       MPI_INTEGER, MPI_MAX, MPI_MIN, mydim_nod2d, MPI_COMM_FESOM, &
       myList_edge2D, myDim_edge2D, myList_nod2D
  USE o_ARRAYS, ONLY: hnode_new
!~   USE output_pdaf, &
!~        ONLY: write_da, write_netcdf_pdaf, write_netcdf_pdaf_ens, &
!~        write_pos_da, write_ens_snapshot, write_pos_da_ens
  USE mod_nc_out_routines, &
       ONLY: netCDF_out
  USE mod_nc_out_variables, &
       ONLY: write_pos_da
  USE obs_TSprof_EN4_pdafomi, &
       ONLY: assim_o_en4_t, assim_o_en4_s, prof_exclude_diff, mean_temp_p
  USE obs_sst_pdafomi, &
       ONLY: assim_o_sst, sst_exclude_ice, sst_exclude_diff, &
             mean_ice_p, mean_sst_p
  USE obs_sss_smos_pdafomi, &
        ONLY: assim_o_sss, sss_exclude_ice, sss_exclude_diff, &
              mean_sss_p
  USE g_clock, &
        ONLY: dayold, yearold, check_fleapyr,daynew,yearnew
  USE mod_assim_pdaf, &
        ONLY: debug_id_nod2

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step
                                     ! (When the routine is called before
                                     ! the analysis, -step is provided.)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL, INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is not generally not initialized in the case of SEIK.
  ! It can be used freely here.
  REAL, INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL, INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

! CALLING SEQUENCE:
! Called by: PDAF_get_state      (as U_prepoststep)
! Called by: PDAF_X_update       (as U_prepoststep)

! *** Local variables ***
  INTEGER :: i, j,k, member, ed    ! Counters
  REAL :: invdim_ens                ! Inverse ensemble size
  REAL :: invdim_ensm1              ! Inverse of ensemble size minus 1 
  REAL :: rmse_p(nfields)                ! PE-local estimated rms errors
  REAL :: rmse(nfields)                  ! Global estimated rms errors
  REAL :: diffm                  ! temporary 
  CHARACTER(len=1) :: typestr       ! Character indicating call type
  REAL :: min_eff_dim_obs, max_eff_dim_obs       ! Stats on effective observation dimensions
  REAL :: min_eff_dim_obs_g, max_eff_dim_obs_g   ! Stats on effective observation dimensions
  REAL :: sum_eff_dim_obs, avg_eff_dim_obs_g     ! Stats on effective observation dimensions
  INTEGER :: dim_nod2D_ice_p        ! PE-local ice node dimension
  LOGICAL :: now_to_write_monthly, flag_print
  INTEGER :: month_iter, whichmonth,fleap,kk,row,nod
  REAL :: tmp_a, tmp_f
  INTEGER, allocatable :: count_lim_salt0_g(:) , count_lim_salt0_p(:)       ! Count how many excessively large updates are limited to treshold
  INTEGER, allocatable :: count_lim_absvel_g(:), count_lim_absvel_p(:)      ! Count how many excessively large updates are limited to treshold
  INTEGER, allocatable :: count_lim_ssh_g(:)   , count_lim_ssh_p(:)         ! Count how many excessively large updates are limited to treshold
  INTEGER, allocatable :: count_lim_tempM2_g(:), count_lim_tempM2_p(:)      ! Count how many excessively large updates are limited to treshold
  
  character(len=100) :: fmt ! format specifier
  
  INTEGER :: myDebug_id(1)
  
  LOGICAL :: debugging_monthlymean = .false.
  
  ! generate the format specifier
  ! fmt = '(a,5x,a,'
  ! do i = 1, dim_ens
  ! fmt = trim(fmt) // 'i7,1x'
  ! end do
  ! fmt = trim(fmt) // ')'

! **********************
! *** INITIALIZATION ***
! **********************

  IF (mype_filter==0) THEN
     IF (step-step_null==0) THEN
       IF (.not.(this_is_pdaf_restart)) THEN
        WRITE (*,'(a, i7,3x,a)') 'FESOM-PDAF', step,'Analyze initial state ensemble'
        WRITE (typestr,'(a1)') 'i'
       ELSE
        WRITE (*,*) 'FESOM-PDAF: This is a PDAF restart.'
       END IF
     ELSE IF (step>0) THEN
        WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Analyze assimilated state ensemble'
        WRITE (typestr,'(a1)') 'a'
     ELSE IF (step<0) THEN
        WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Analyze forecast state ensemble'
        WRITE (typestr,'(a1)') 'f'
     END IF
  END IF ! IF (mype_filter==0)

  IF (.not. ALLOCATED(var_p)) ALLOCATE(var_p(dim_p))
  IF (.not. ALLOCATED(state_fcst)) ALLOCATE(state_fcst(dim_p,dim_ens))
  IF (.not. ALLOCATED(std_p)) ALLOCATE(std_p(dim_p))
  
  IF (.not. ALLOCATED(monthly_state_a))    ALLOCATE(monthly_state_a(dim_p))
  IF (.not. ALLOCATED(monthly_state_f))    ALLOCATE(monthly_state_f(dim_p))
  IF (.not. ALLOCATED(monthly_state_m))    ALLOCATE(monthly_state_m(dim_p))

  
  IF (.not. ALLOCATED(count_lim_salt0_g))  ALLOCATE(count_lim_salt0_g  (dim_ens))
  IF (.not. ALLOCATED(count_lim_salt0_p))  ALLOCATE(count_lim_salt0_p  (dim_ens))
  IF (.not. ALLOCATED(count_lim_absvel_g)) ALLOCATE(count_lim_absvel_g (dim_ens))
  IF (.not. ALLOCATED(count_lim_absvel_p)) ALLOCATE(count_lim_absvel_p (dim_ens))
  IF (.not. ALLOCATED(count_lim_ssh_g))    ALLOCATE(count_lim_ssh_g    (dim_ens))
  IF (.not. ALLOCATED(count_lim_ssh_p))    ALLOCATE(count_lim_ssh_p    (dim_ens))
  IF (.not. ALLOCATED(count_lim_tempM2_g)) ALLOCATE(count_lim_tempM2_g (dim_ens))
  IF (.not. ALLOCATED(count_lim_tempM2_p)) ALLOCATE(count_lim_tempM2_p (dim_ens))

  ! Initialize numbers
  rmse_p = 0.0
  rmse   = 0.0
  invdim_ens = 1.0 / REAL(dim_ens)  
  invdim_ensm1 = 1.0 / REAL(dim_ens-1)  


! ****************************
! *** Perform pre/poststep ***
! ****************************


   ! *** Correction ***
   
    IF ((step - step_null) < 0) THEN
    ! *** store forecast state fields temporarily to compare with analysis afterwards ***
    state_fcst = ens_p
    
    ELSE IF ((step - step_null) > 0) THEN
   ! *** correcting assimilated state fields ***
   
   ! *** salinity must be > 0 ***
   count_lim_salt0_p = 0
   
   DO member = 1, dim_ens
      DO i = 1, dim_fields(id% salt)
           IF (ens_p(i+ offset(id% salt),member) < 0.0D0) THEN
               ens_p(i+ offset(id% salt),member) = 0.0D0
               
               count_lim_salt0_p(member) = count_lim_salt0_p(member)+1
           END IF
      END DO
   END DO
   
   CALL MPI_Allreduce(count_lim_salt0_p, count_lim_salt0_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   !IF (mype_filter == 0) &
   !         WRITE(*, fmt) 'FESOM-PDAF', &
   !         '--- Updated salinity limited to zero: ', (count_lim_salt0_g(member), member = 1, dim_ens)
   
     
   ! *** Analysed horizontal velocities must be less than a doubling of forecast velocity ***
!~    count_lim_absvel_p = 0
   
!~    DO member = 1, dim_ens
!~    DO i = 1, dim_fields(id% u)
   
!~       ! absolute velocities
!~       tmp_a =   ens_p(i+ offset(id% u),member) * ens_p(i+ offset(id% u),member) &
!~               + ens_p(i+ offset(id% v),member) * ens_p(i+ offset(id% v),member)
                    
!~       tmp_f =   state_fcst(i+ offset(id% u),member) * state_fcst(i+ offset(id% u),member) &
!~               + state_fcst(i+ offset(id% v),member) * state_fcst(i+ offset(id% v),member)
                    
!~       IF (tmp_a > 4.0*tmp_f) THEN
!~          !ens_p(i+ offset(id% u),member) = state_fcst(i+ offset(id% u),member) * 2.0 ! *tmp_f/max(tmp_a,0.0001)
!~          !ens_p(i+ offset(id% v),member) = state_fcst(i+ offset(id% v),member) * 2.0 ! *tmp_f/max(tmp_a,0.0001)
         
!~          !count_lim_absvel_p(member) = count_lim_absvel_p(member)+1
         
!~       END IF
!~    END DO
!~    END DO
   
!~    CALL MPI_Allreduce(count_lim_absvel_p, count_lim_absvel_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   !IF (mype_filter == 0) &
   !         WRITE(*, fmt) 'FESOM-PDAF', &
   !         '--- Velocity updates limited to doubling: ', (count_lim_absvel_g(member), member = 1, dim_ens)
   
   
   ! *** SSH state update must be <= 2*sigma
   count_lim_ssh_p = 0
   
   DO member = 1, dim_ens
   DO i = offset(id% SSH)+1, offset(id% SSH)+dim_fields(id% SSH)
       diffm = ens_p(i,member) - state_fcst(i,member)
       IF (ABS(diffm) > 2.0*std_p(i)) THEN
           ens_p(i,member) = state_fcst(i,member) + SIGN(2.0*std_p(i),diffm)
           
           count_lim_ssh_p(member) = count_lim_ssh_p(member)+1
       END IF
   END DO
   END DO
   
   CALL MPI_Allreduce(count_lim_ssh_p, count_lim_ssh_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   !IF (mype_filter == 0) &
   !         WRITE(*, fmt) 'FESOM-PDAF', &
   !         '--- SSH updates limited to 2x standard deviation: ', (count_lim_ssh_g(member), member = 1, dim_ens)


   ! *** temperature must be > -2 degC ***
   count_lim_tempM2_p = 0
   
   DO member = 1, dim_ens
      DO i = 1, dim_fields(id% salt)
           IF (ens_p(i+ offset(id% temp),member) < -2.0) THEN
               ens_p(i+ offset(id% temp),member) = -2.0
               
               count_lim_tempM2_p(member) = count_lim_tempM2_p(member)+1
           END IF
      END DO
   END DO
   
   CALL MPI_Allreduce(count_lim_tempM2_p, count_lim_tempM2_g, dim_ens, MPI_INTEGER, MPI_SUM, COMM_filter, MPIerr)
   !IF (mype_filter == 0) &
   !         WRITE(*, fmt) 'FESOM-PDAF', &
   !         '--- Updated temperature limited to -2 degC: ', (count_lim_tempM2_g(member), member = 1, dim_ens)

   END IF ! Corrections


  ! *** Compute mean state
  IF (mype_filter==0) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', '--- compute ensemble mean'
  
  ! Local: 
  state_p = 0.0
  DO member = 1, dim_ens
     DO i = 1, dim_p
        state_p(i) = state_p(i) + ens_p(i,member)
     END DO
  END DO
  state_p(:) = invdim_ens * state_p(:)


! *********************************************************************
! *** Store ensemble mean values for observation exclusion criteria ***
! *********************************************************************

  IF ((assim_o_sst .OR. assim_o_sss) .AND. step<0) THEN
! (forecast phase)
!
     ! sea-ice concentration
     IF (sst_exclude_ice .OR. sss_exclude_ice) THEN 
        IF (ALLOCATED(mean_ice_p)) DEALLOCATE(mean_ice_p)
        ALLOCATE (mean_ice_p(dim_fields(id% a_ice)))
        mean_ice_p = state_p(offset(id% a_ice)+ 1 : &
                             offset(id% a_ice)+ dim_fields(id% a_ice))
     END IF

     ! SST
     IF (sst_exclude_ice .OR. sst_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_sst_p)) DEALLOCATE(mean_sst_p)
        ALLOCATE (mean_sst_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sst_p(i) = state_p(offset(id% temp) + (i-1) * (mesh_fesom%nl-1) + 1)
        END DO
     END IF
     
     ! SSS
     IF (sss_exclude_ice .OR. sss_exclude_diff > 0.0) THEN
        IF (ALLOCATED(mean_sss_p)) DEALLOCATE(mean_sss_p)
        ALLOCATE (mean_sss_p(myDim_nod2D))
        DO i = 1, myDim_nod2D
          mean_sss_p(i) = state_p(offset(id% salt) + (i-1) * (mesh_fesom%nl-1) + 1)
        END DO
     END IF

  END IF

  ! 3D temperature field
  IF (        (assim_o_en4_t .OR. assim_o_en4_s) &
      .AND. &
              (step < 0) &
      .AND. &
              (prof_exclude_diff > 0.0)) THEN

     ! Store mean temperature for profile assimilation
     IF (ALLOCATED(mean_temp_p)) DEALLOCATE(mean_temp_p)
     ALLOCATE (mean_temp_p(myDim_nod2D))
     mean_temp_p = state_p(offset(id%temp)+1 : offset(id%temp)+dim_fields(id%temp))

  END IF


! *****************************************************************
! *** Compute estimate RMS ensemble spread for different fields ***
! *****************************************************************

  ! *** Compute local sampled variances of ensemble state vector ***

  var_p(:) = 0.0D0
  DO member = 1, dim_ens
     DO j = 1, dim_p
        var_p(j) = var_p(j) + &
             (ens_p(j, member) - state_p(j)) * (ens_p(j, member) - state_p(j))
     END DO
  END DO
  var_p(:) = invdim_ensm1 * var_p(:)

  ! IF ((step - step_null) < 0) ! if analysis: compute STD
  std_p = SQRT(var_p)     

  ! *** Compute RMS ensemble spread ***

  DO j = 1, nfields
       DO i = 1, dim_fields(j)
          rmse_p(j) = rmse_p(j) + var_p(i + offset(j))
       END DO
       rmse_p(j) = rmse_p(j) / REAL(dim_fields_glob(j))
  END DO

  ! Global sum of RMS errors
  CALL MPI_Allreduce (rmse_p, rmse, nfields, MPI_DOUBLE_PRECISION, MPI_SUM, &
       MPI_COMM_FESOM, MPIerr)

  rmse = SQRT(rmse)

  ! Display RMS errors
  IF (mype_filter==0) THEN
     WRITE (*,'(a, 10x,a)') &
          'FESOM-PDAF', 'RMS error according to sampled covariance:'
     WRITE (*,'(a,7x,    a14, a14,a14,a14, a14,  a14,            /a, 10x,66a)') &
          'FESOM-PDAF', 'SSH','U','V','W','temp','salt', &
          'FESOM-PDAF', ('-',i=1,66)
     WRITE (*,'(a,10x,  6es14.4, 1x,a5,a1,/a, 10x,66a)') &
          'FESOM-PDAF', rmse(1), rmse(2), rmse(3), rmse(4), rmse(5), rmse(6), 'RMSe-', typestr,&
          'FESOM-PDAF', ('-',i=1,66)
  END IF


! ***************************************************************
! *** Compute statistics for effective observation dimensions ***
! ***************************************************************

  IF (loctype==1 .AND. step>0) THEN

     max_eff_dim_obs = 0.0
     min_eff_dim_obs = 1.0e16
     sum_eff_dim_obs = 0.0

     DO i = 1, myDim_nod2D
        IF (eff_dim_obs(i) > max_eff_dim_obs) max_eff_dim_obs = eff_dim_obs(i)
        IF (eff_dim_obs(i) < min_eff_dim_obs) min_eff_dim_obs = eff_dim_obs(i)
        sum_eff_dim_obs = sum_eff_dim_obs + eff_dim_obs(i)
     END DO
     IF (npes_filter>1) THEN
        CALL MPI_Reduce(sum_eff_dim_obs, avg_eff_dim_obs_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
             0, COMM_filter, MPIerr)
        CALL MPI_Reduce(max_eff_dim_obs, max_eff_dim_obs_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
             0, COMM_filter, MPIerr)
        CALL MPI_Reduce(min_eff_dim_obs, min_eff_dim_obs_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
             0, COMM_filter, MPIerr)
     ELSE
        ! This is a work around for working with nullmpi.F90
        avg_eff_dim_obs_g = sum_eff_dim_obs
        min_eff_dim_obs_g = min_eff_dim_obs
        max_eff_dim_obs_g = max_eff_dim_obs
     END IF

     IF (mype_filter==0) THEN
        avg_eff_dim_obs_g = avg_eff_dim_obs_g / REAL(mesh_fesom%nod2d)

        WRITE (*, '(a, 8x, a)') &
             'FESOM-PDAF', '--- Effective observation dimensions for local analysis:'
        WRITE (*, '(a, 12x, a, f12.2)') &
             'FESOM-PDAF', 'min. effective observation dimension:       ', min_eff_dim_obs_g
        WRITE (*, '(a, 12x, a, f12.2)') &
             'FESOM-PDAF', 'max. effective observation dimension:       ', max_eff_dim_obs_g
        WRITE (*, '(a, 12x, a, f12.2)') &
             'FESOM-PDAF', 'avg. effective observation dimension:       ', avg_eff_dim_obs_g
     END IF
  END IF
  
! ***************************
! *** Compute daily means ***
! ***************************

  IF ((step - step_null) > 0) THEN
     timemean = timemean + state_p / delt_obs_ocn
  ENDIF

! *****************************
! *** Compute monthly means ***
! *****************************

  debugging_monthlymean = .false.

  IF (write_monthly_mean) THEN
  
      now_to_write_monthly = .FALSE.
      call check_fleapyr(yearold, fleap)
      
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'dayold: ', dayold
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'daynew: ', daynew
      
      ! set "now_to_write_monthly" at last day of month
      find_month: DO month_iter=1,12
        if (dayold <= endday_of_month_in_year(fleap,month_iter)) THEN
            whichmonth = month_iter
            if (dayold == endday_of_month_in_year(fleap,month_iter)) THEN
               now_to_write_monthly = .TRUE.
            endif
            exit find_month
        endif
      ENDDO find_month
      
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'whichmonth: ', whichmonth
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'now_to_write_monthly: ', now_to_write_monthly

      ! reset monthly_state to zero at first day of month
      monthly_state_p_0: DO month_iter=1,12
        if (dayold == startday_of_month_in_year(fleap,month_iter)) THEN
          IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'resetting monthly mean to zero.'
          IF ((step - step_null) > 0) THEN
          ! *** assimilated state fields ***
            monthly_state_a= 0.0D0
            monthly_state_m= 0.0D0
          ELSE IF ((step - step_null) < 0) THEN
          ! *** forecasted state fields ***
            monthly_state_f= 0.0D0
          END IF
          exit monthly_state_p_0
        endif
      ENDDO monthly_state_p_0
      
      IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'num_day_in_month: ', num_day_in_month(fleap,whichmonth)
      
      ! include daily state into monthly mean
      IF ((step - step_null) > 0) THEN
      ! *** analyzed state fields ***
        IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'adding to mean.'
        monthly_state_a = monthly_state_a + state_p /num_day_in_month(fleap,whichmonth)
        monthly_state_m = monthly_state_m + timemean/num_day_in_month(fleap,whichmonth)
      ELSE IF ((step - step_null) < 0) THEN
      ! *** forecasted state fields ***
        IF (debugging_monthlymean .and. mype_filter==0) WRITE(*,*) 'step ', step, 'adding to mean.'
        monthly_state_f = monthly_state_f + state_p/num_day_in_month(fleap,whichmonth)
      END IF

   ENDIF
   
! **************************
! *** Write output files ***
! **************************

! *** write initial state fields ***
IF ((step - step_null)==0 .and. ( .not. this_is_pdaf_restart)) THEN
      CALL netCDF_out('i',daynew,step,state_p,ens_p,rmse)
ENDIF

IF (.not. write_monthly_mean) THEN
  ! daily output
  write_pos_da = daynew
  IF ((step - step_null) < 0) THEN
        ! *** write forecast state fields ***
        CALL netCDF_out('f',daynew,step,state_p,ens_p,rmse)
  ELSE IF ((step - step_null) > 0) THEN
        ! *** write analysis state fields ***
        CALL netCDF_out('a',daynew,step,state_p, ens_p,rmse)
        CALL netCDF_out('m',daynew,step,timemean,ens_p,rmse)
  END IF

ELSEIF (write_monthly_mean .and. now_to_write_monthly) THEN
  ! monthly output
  write_pos_da = whichmonth
  IF ((step - step_null) < 0) THEN
        ! *** write forecast state fields ***
        CALL netCDF_out('f',write_pos_da,step,monthly_state_f,ens_p,rmse)
  ELSE IF ((step - step_null) > 0) THEN
        ! *** write analysis state fields ***
        CALL netCDF_out('a',write_pos_da,step,monthly_state_a,ens_p,rmse)
        CALL netCDF_out('m',write_pos_da,step,monthly_state_m,ens_p,rmse)
  END IF
END IF ! write_monthly_mean

! ********************
! *** finishing up ***
! ********************
  IF ((step - step_null) < 0) DEALLOCATE(var_p)
  IF ((step - step_null) > 0) DEALLOCATE(std_p)
!~   DEALLOCATE(monthly_state_a, monthly_state_f)
  
!~ ENDIF ! this_is_pdaf_restart



END SUBROUTINE prepoststep_pdaf
