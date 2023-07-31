SUBROUTINE init_dim_obs_f_proffile_pdaf(startyr, endyr)

   ! !DESCRIPTION:
   ! Create distributed observation files for EN4 data. Each process reads
   ! the original EN4 data, selects those observations that belong to the 
   ! process and pass the quality checks. Then the observations are written
   ! into one netCDF file per process including the node information required
   ! for the observation operator.
   !
   ! !REVISION HISTORY:
   ! 2019-03 - Qi Tang      - Initial code based on init_dim_obs_f_prof_pdaf
   ! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
   ! 2022-05 - Frauke B     - Adapted for FESOM2.1
   !
   ! !USES:
   USE mod_parallel_pdaf, &
      ONLY: mype_filter, COMM_filter, abort_parallel
   USE mod_assim_pdaf, &
      ONLY: path_obs_rawprof, file_rawprof_prefix,  &
            file_rawprof_suffix, &
            num_day_in_month, &
            mesh_fesom
            ! mesh_fesom % coord_nod2D        ! vertex coordinates in radian measure
            ! mesh_fesom % nlevels            ! number of levels at (below) elem     considering bottom topography
            ! mesh_fesom % nlevels_nod2D      ! number of levels at (below) vertices considering bottom topography
            ! mesh_fesom % nl                 ! number of levels not considering bottom topography
            ! mesh_fesom % elem2D_nodes(:,e)  ! 3 nodes of element e
            ! mesh_fesom % Z                  ! mid-layer depths not considering bottom topography
            ! mesh_fesom % zbar               ! upper layer bounds not considering bottom topography
   USE obs_TSprof_EN4_pdafomi, &
      ONLY: assim_o_en4_t, assim_o_en4_s, rms_obs_T, rms_obs_S, &
            path_obs_prof, file_prof_prefix, file_prof_suffix, &
            prof_exclude_diff, file_syntobs_prof, bias_obs_prof
   USE g_PARSUP, &
      ONLY: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_SUM,  &
            MPIerr, myDim_nod2D, MPI_COMM_FESOM
   USE o_ARRAYS, &
      ONLY: z_n, zbar_n
   USE o_param, &
      ONLY: rad
   use g_clock, &
      ONLY: fleapyear, check_fleapyr
   USE IFPORT, &
      ONLY: MAKEDIRQQ

     IMPLICIT NONE
     INCLUDE 'netcdf.inc'
   
!    !ARGUMENTS:
      INTEGER, INTENT(in)    :: startyr, endyr      
   
   ! !CALLING SEQUENCE:
   ! Called by: init_dim_obs_f_pdaf
   !EOP
   
     ! Local variables
     CHARACTER(len=120) :: outpath, outfile, ncfile_out, ncid_out
     INTEGER :: steps
     INTEGER :: month_day (12)
     CHARACTER(len=150) :: attstr
     CHARACTER(len=10) :: mype_string
     INTEGER :: dimid_time, dimid_two, dimid_three, dimid_nobs_temp, dimid_nobs_sal, dimid_ndepth
     INTEGER :: id_time, id_nprof_temp, id_nprof_sal
     INTEGER :: id_nflagged_temp, id_nflagged_sal
     INTEGER :: id_nobs_temp, id_nobs_sal, id_depth_temp, id_depth_sal
     INTEGER :: id_temp_out, id_sal_out
     INTEGER :: id_node_temp, id_node_sal, id_layer_temp, id_layer_sal
     INTEGER :: id_coord_temp, id_coord_sal, id_offset_temp, id_offset_sal
     INTEGER :: startv, countv, startv2(2), countv2(2)
     INTEGER :: i                            ! Counters
     INTEGER :: stat(100)                    ! Status for NetCDF functions
     INTEGER :: dimids(2)
     INTEGER :: offset_temp, offset_sal, cnt_day
     INTEGER :: month_index
     
     INTEGER :: ncid_in, id_nprof, id_nlevels, id_lat, id_lon, id_judate, id_refjudate   ! IDs
     INTEGER :: id_temp, id_sal, id_depth   ! IDs for variables
     INTEGER :: id_position_qc, id_pro_temp_qc, id_pro_sal_qc, id_temp_qc, id_sal_qc     ! IDs for quality check 
     INTEGER :: n_prof, n_levels, s, nprof_day, nprof_day_p     ! Counters
     INTEGER :: j, k, nn, cnt_prof_temp, cnt_prof_sal           ! Counters
     REAL, ALLOCATABLE :: obs_coord(:,:), obs_coord_day(:,:)    ! Coordinates for the observations
     REAL, ALLOCATABLE :: judate(:)                             ! Julian date
     CHARACTER(len=14) :: refjudate_char                        ! Reference Julian date
     REAL :: refjudate, JD, JT                                  ! Julian date
     REAL:: sum_temp, sum_sal
     INTEGER :: Y, M, D, L, N, YEAR, MONTH, DAY, HOUR, MINUTE, SECOND   !Gregorian date
     INTEGER, ALLOCATABLE :: index_obs(:)
     INTEGER, ALLOCATABLE :: index_ele_day_p(:),index_pro_day_p(:)
     INTEGER :: element, ilayer, nlayer, cnt_cal
     INTEGER :: cnt_temp, cnt_sal, cnt_temp_g, cnt_sal_g
     REAL(4), ALLOCATABLE :: temp(:,:), sal(:,:), depth(:,:)                ! Values read from file for full month
     REAL, ALLOCATABLE :: temp_day(:,:), sal_day(:,:), depth_day(:,:)       ! values for a certain day
     CHARACTER(len=1), ALLOCATABLE :: pos_qc(:)                             ! position quality check read from file
     CHARACTER(len=1), ALLOCATABLE :: pro_temp_qc(:), pro_sal_qc(:)         ! Variables for quality check
     CHARACTER(len=1), ALLOCATABLE :: pos_qc_day(:)                         ! Variables for quality check
     CHARACTER(len=1), ALLOCATABLE :: pro_temp_qc_day(:), pro_sal_qc_day(:) ! Variables for quality check
     CHARACTER(len=1), ALLOCATABLE :: temp_qc(:,:), sal_qc(:,:)             ! Variables for quality check
     CHARACTER(len=1), ALLOCATABLE :: temp_qc_day(:,:), sal_qc_day(:,:)     ! Variables for quality check 
     CHARACTER(len=1), PARAMETER :: one='1'
     INTEGER :: n2d_prof(3)
     INTEGER, ALLOCATABLE :: n2d_temp(:,:), n2d_sal(:,:)
     INTEGER, ALLOCATABLE :: nl1_temp(:), nl1_sal(:)
     REAL, ALLOCATABLE :: depth_temp (:), temp_obs(:)
     REAL, ALLOCATABLE :: depth_sal (:), sal_obs(:)
     REAL, ALLOCATABLE :: obs_coord_day_temp(:,:), obs_coord_day_sal(:,:)
     INTEGER :: d_in_month
     INTEGER :: flag_n2d                               ! flag for checking n2d_sal/n2d_temp is within mydim_nod2d
     INTEGER :: cnt_en2d_sal
     INTEGER :: iprof_day,iprof_day_p
     CHARACTER(len=100) :: prof_file = ''              ! Complete name of profile observation file without path
     INTEGER :: cnt_pro_temp_g, cnt_pro_sal_g          ! Total number of profiles at one day
     integer :: cnt_flagged_temp_p, cnt_flagged_sal_p  ! PE-local number of flagged profiles
     integer :: cnt_flagged_temp, cnt_flagged_sal      ! Global number of flagged profiles
     integer :: yyyy
     character (4) :: cyyyy
     LOGICAL :: FEXIST=.FALSE.
     INTEGER :: RES
     LOGICAL :: debugging

     DO yyyy = startyr, endyr
      ! **********************
      ! *** Initialization ***
      ! **********************

        ! Define numbers
        call check_fleapyr(yyyy, fleapyear)
        IF (fleapyear==1) then
             steps = 366
             month_day = num_day_in_month(1,:)
        else
             steps = 365
             month_day = num_day_in_month(0,:)
        endif
      
        ! Create the PE_local netCDF output file
        ! Path to and name stub of output files
        WRITE(cyyyy,'(i4.4)') yyyy      

        outpath = TRIM(path_obs_prof)//cyyyy//'/'
        IF (mype_filter==0) THEN
            INQUIRE(DIRECTORY=TRIM(outpath),EXIST=FEXIST)
            IF(.not.FEXIST) THEN
                res=MAKEDIRQQ(TRIM(outpath))
                IF (res) THEN
                    WRITE (*,*) 'Profile created',outpath
                ELSE
                    WRITE (*,*) 'Failed to create directory:',outpath
                END IF
            ENDIF
        ENDIF
        FEXIST = .FALSE.
        CALL MPI_BARRIER(MPI_COMM_FESOM, MPIerr)

        ncfile_out = TRIM(outpath)//TRIM(file_prof_prefix)

        WRITE(mype_string,'(i4.4)') mype_filter      
      
        IF (mype_filter==0) THEN
           WRITE (*,'(a,5x,a)')   'FESOM-PDAF', 'Generate PE-local EN4 profile observations'
           WRITE (*,'(a,5x,a,a)') 'FESOM-PDAF', '--- Write distributed fields to files: ', &
                                   TRIM(ncfile_out)//'XXXX.nc'
        END IF
   
   ! *********************************************
   ! *** Generate PE_local profile netCDF file ***
   ! *********************************************
   
     ! Initialization
     s = 1
     stat(s) = NF_CREATE(TRIM(ncfile_out)//TRIM(mype_string)//'.nc', NF_netcdf4, ncid_out)
     attstr  = 'daily EN4 profile observations from year '//cyyyy
     s = s + 1
     stat(s) = NF_PUT_ATT_TEXT(ncid_out, NF_GLOBAL, 'title', LEN_TRIM(attstr), TRIM(attstr))
     ! Define dimensions
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'time', steps, dimid_time)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'two', 2, dimid_two)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'three', 3, dimid_three)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'n_obs_temp', NF_UNLIMITED, dimid_nobs_temp)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'n_obs_sal', NF_UNLIMITED, dimid_nobs_sal)
     s = s + 1
     stat(s) = NF_DEF_DIM(ncid_out, 'n_depth', mesh_fesom%nl - 1, dimid_ndepth)  
   
     DO i = 1,  s - 1
        IF (stat(i) /= NF_NOERR) then
           WRITE(*, *) 'NetCDF error in profile dimension definitions, no.', i, mype_filter
           WRITE (*,*) NF_STRERROR(stat(s))
        end IF
     END DO
   
     ! define variables   ------ 1 dimension : time, nprof, nobs, depth, layer, temperature, salinity
     s = 1
     stat(s) = NF_DEF_VAR(ncid_out, 'day', NF_INT, 1, dimid_time, id_time)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nprof_temp', NF_INT, 1, dimid_time, id_nprof_temp)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nprof_sal', NF_INT, 1, dimid_time, id_nprof_sal)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nobs_temp', NF_INT, 1, dimid_time, id_nobs_temp)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nobs_sal', NF_INT, 1, dimid_time, id_nobs_sal)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nflagged_temp',NF_INT, 1, dimid_time, id_nflagged_temp)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'nflagged_sal',NF_INT, 1, dimid_time, id_nflagged_sal)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'offset_temp', NF_INT, 1, dimid_time, id_offset_temp)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'offset_sal', NF_INT, 1, dimid_time, id_offset_sal)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'depth_temp', NF_DOUBLE, 1, dimid_nobs_temp, id_depth_temp)
     s = s + 1
     stat(s) = NF_DEF_VAR(ncid_out, 'depth_sal', NF_DOUBLE, 1, dimid_nobs_sal, id_depth_sal)
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'temperature', NF_DOUBLE, 1, dimid_nobs_temp, id_temp_out)
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'salinity', NF_DOUBLE, 1, dimid_nobs_sal, id_sal_out)
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'nl1_temp', NF_INT, 1, dimid_nobs_temp, id_layer_temp)
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'nl1_sal', NF_INT, 1, dimid_nobs_sal, id_layer_sal)
   
     ! define variables   ------ 2 dimensions : node, coordinate
   
     dimids(1) = dimid_three
     dimids(2) = dimid_nobs_temp
     
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'n2d_temp', NF_INT, 2, dimids, id_node_temp)
   
     dimids(1) = dimid_two
     
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'coordinate_temp', NF_DOUBLE, 2, dimids, id_coord_temp)
   
     dimids(1) = dimid_three
     dimids(2) = dimid_nobs_sal
     
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'n2d_sal', NF_INT, 2, dimids, id_node_sal)
   
     dimids(1) = dimid_two
   
     s = s + 1
     stat(s)=NF_DEF_VAR(ncid_out, 'coordinate_sal', NF_DOUBLE, 2, dimids, id_coord_sal)
     
     s = s + 1
     stat(s) = NF_ENDDEF(ncid_out)
   
     DO i = 1, s
        IF (stat(i) /= NF_NOERR) THEN
           WRITE(*, *) 'NetCDF error in profile variable definitions, no.', i
           WRITE (*,*) NF_STRERROR(stat(s))
        END IF
     END DO
   
   
   ! **********************************************************************************
   ! *** Loop through monthly original EN4 file and generate and write local fields ***
   ! **********************************************************************************
     
   
     offset_temp=1
     offset_sal=1  
     cnt_day=1
     
     LOOP_month: DO month_index = 1, 12
   
        ! *******************************************
        ! *** Read profile data for current month ***
        ! *******************************************
   
        WRITE (prof_file,'(A,I2.2,A)') TRIM(file_rawprof_prefix)//cyyyy, month_index, TRIM(file_rawprof_suffix)
   
        ! Open the EN4 file
        s = 1
        stat(s) = NF_OPEN(TRIM(TRIM(path_obs_rawprof)//TRIM(prof_file)), NF_NOWRITE, ncid_in)
        s = s + 1
   
        ! *** Get dimensions ***
   
        stat(s) = NF_INQ_DIMID(ncid_in, 'N_PROF', id_nprof)
        s = s + 1
        stat(s) = NF_INQ_DIMLEN(ncid_in, id_nprof, n_prof)
        s = s + 1
        stat(s) = NF_INQ_DIMID(ncid_in, 'N_LEVELS', id_nlevels)
        s = s + 1
        stat(s) = NF_INQ_DIMLEN(ncid_in, id_nlevels, n_levels)
   
        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading dimensions, no.', i
        END DO
   
        ! *** Real Julian date ***
   
        ! Allocate array which stores the julian date
        ALLOCATE(judate(n_prof))
   
        ! Get variable -- Julian date IDs
   
        s = 1
        stat(s) = NF_INQ_VARID(ncid_in, 'JULD', id_judate)
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'REFERENCE_DATE_TIME', id_refjudate)
   
        ! Read Julian date
   
        s = s + 1
        stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_judate, judate)
   
        ! Read reference Julian date
   
        s = s + 1
        stat(s) = NF_GET_VAR_TEXT(ncid_in, id_refjudate, refjudate_char)
   
        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading julian date, no.', i
        END DO
   
        READ(refjudate_char(1:4),*) Y
        READ(refjudate_char(5:6),*) M
        READ(refjudate_char(7:8),*) D
        READ(refjudate_char(9:10),*) HOUR
        READ(refjudate_char(11:12),*) MINUTE
        READ(refjudate_char(13:14),*) SECOND
   
        ! Convert reference Gregorian date into Julian date
   
        refjudate =   (D-32075 + 1461*(Y+4800+(M-14)/12)/4 &
                    + 367*(M-2-(M-14)/12*12)/12 &
                    - 3*((Y+4900+(M-14)/12)/100)/4) &
                    + ((HOUR + (MINUTE + SECOND/60.D0)/60.D0)/24.D0)
   
   
        ! ***************************************
        ! *** Read coordinates and profile QC ***
        ! ***************************************
   
        ! *** Coordinates ***
   
        ALLOCATE(obs_coord(n_prof, 2))
   
        s = 1
        stat(s) = NF_INQ_VARID(ncid_in, 'LATITUDE', id_lat)
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'LONGITUDE', id_lon)
        s = s + 1
        stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_lat, obs_coord(:, 2))
        s = s + 1
        stat(s) = NF_GET_VAR_DOUBLE(ncid_in, id_lon, obs_coord(:, 1))
   
   
        ! *** Position quality check information ***
   
        ALLOCATE(pos_qc(n_prof))
   
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'POSITION_QC', id_position_qc)
        s = s + 1
        stat(s) = NF_GET_VAR_TEXT(ncid_in, id_position_qc, pos_qc)
   
        DO i = 1,  s
           IF (stat(i) /= NF_NOERR) then
              WRITE(*, *) 'NetCDF error in reading coordinates and position QC , no.', i
              WRITE (*,*) NF_STRERROR(stat(s))
           end IF
        END DO
   
   
        ! *********************************************
        ! *** Read temperature, salinity, and depth ***
        ! *********************************************
   
        ! *** Depth ***
   
        ALLOCATE(depth(n_levels, n_prof))
   
        s = 1
        stat(s) = NF_INQ_VARID(ncid_in, 'DEPH_CORRECTED', id_depth)
        s = s + 1
        stat(s) = NF_GET_VAR_REAL(ncid_in, id_depth, depth)
   
      
        ! *** Temperature ***
   
        ALLOCATE(temp(n_levels, n_prof))
        ALLOCATE(pro_temp_qc(n_prof))
        ALLOCATE(temp_qc(n_levels, n_prof))
   
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'POTM_CORRECTED', id_temp)
        s = s + 1
        stat(s) = NF_GET_VAR_REAL(ncid_in, id_temp, temp)
   
        ! Read temperature profile quality check information
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'PROFILE_POTM_QC', id_pro_temp_qc)
        s = s + 1
        stat(s) = NF_GET_VAR_TEXT(ncid_in, id_pro_temp_qc, pro_temp_qc)
   
        ! Read level-temperature quality check information
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'POTM_CORRECTED_QC', id_temp_qc)
        s = s + 1
        stat(s) = NF_GET_VAR_TEXT(ncid_in, id_temp_qc, temp_qc)
   
   
        ! *** Salinity ***
   
        ALLOCATE(sal(n_levels, n_prof))
        ALLOCATE(pro_sal_qc(n_prof))
        ALLOCATE(sal_qc(n_levels, n_prof))
   
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'PSAL_CORRECTED', id_sal)
        s = s + 1
        stat(s) = NF_GET_VAR_REAL(ncid_in, id_sal, sal)
   
        ! Read salinity profile quality check information
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'PROFILE_PSAL_QC', id_pro_sal_qc)
        s = s + 1
        stat(s) = NF_GET_VAR_TEXT(ncid_in, id_pro_sal_qc, pro_sal_qc)
   
        ! Read level-salinity quality check information
        s = s + 1
        stat(s) = NF_INQ_VARID(ncid_in, 'PSAL_CORRECTED_QC', id_sal_qc)
        s = s + 1
        stat(s) = NF_GET_VAR_TEXT(ncid_in, id_sal_qc, sal_qc)
   
        s = s + 1
        stat(s) = NF_CLOSE(ncid_in)
   
        DO k = 1,  s
           IF (stat(k) /= NF_NOERR) &
                WRITE(*, *) 'NetCDF error in reading fields and quality control, no.', i
        END DO
   
   
        ! **************************************************
        ! *** Process profiles for all days of the month ***
        ! **************************************************
   
        ALLOCATE(index_obs(n_prof))
   
        ! Loop over day
        LOOP_day: DO d_in_month = 1, month_day(month_index)
   
   
           IF (mype_filter==0) WRITE (*,'(a,5x,a,i3,a,i3)') &
                'FESOM-PDAF', 'Process day:', d_in_month, ' month', month_index
   
           
           ! ************************************
           ! *** Find profiles at current day ***
           ! ************************************
   
           nprof_day = 0
   
           LOOP_prof: DO i = 1, n_prof
   
              ! *** Convert Julian date into Gregorian date ***
   
              JD = judate(i) + refjudate
              L = INT(JD)+68569
              N = 4*L/146097
              L = L-(146097*N+3)/4
              Y = 4000*(L+1)/1461001
              L = L-1461*Y/4+31
              M = 80*L/2447
              D = L-2447*M/80
              L = M/11
              M = M+2-12*L
              Y = 100*(N-49)+Y+L
   
              YEAR = Y
              MONTH = M
              DAY = D
              JT = DMOD(JD,1.D0)*24.D0
              HOUR = INT(JT)
              JT = DMOD(JT,1.D0)*60.D0
              MINUTE = INT(JT)
              JT = DMOD(JT,1.D0)*60.D0
              SECOND = NINT(JT)
   
              IF (SECOND == 60) THEN
                 SECOND = SECOND-60
                 MINUTE = MINUTE+1
              END IF
   
              ! Check if profile is for current day
              IF (DAY == d_in_month) THEN
                 nprof_day = nprof_day + 1
                 index_obs(nprof_day) = i
              END IF
   
           END DO LOOP_prof
   
   
           ! *********************************************
           ! *** Initialize arrays for the current day ***
           ! *********************************************
     
           ! *** Coordinates ***
   
           ALLOCATE(obs_coord_day(2, nprof_day))
   
           DO i = 1, nprof_day
              ! convert to radians and transpose the coordinate array
              obs_coord_day(1, i) = obs_coord(index_obs(i),1)*rad
              obs_coord_day(2, i) = obs_coord(index_obs(i),2)*rad
           END DO
   
   
           ! *** Position quality check information ***
   
           ALLOCATE(pos_qc_day(nprof_day))
   
           DO i = 1, nprof_day
              pos_qc_day(i) = pos_qc(index_obs(i))
           END DO
   
   
           ! *** Depth ***
   
           ALLOCATE(depth_day(nprof_day,n_levels))
   
           DO i = 1, nprof_day
              depth_day(i,:) = REAL(depth(:,index_obs(i)),8)
           END DO
   
   
           ! *** Temperature ***
   
           ALLOCATE(temp_day(nprof_day,n_levels))
           ALLOCATE(pro_temp_qc_day(nprof_day))
           ALLOCATE(temp_qc_day(nprof_day,n_levels))
   
           ! Temperature values
           DO i = 1, nprof_day
              temp_day(i,:) = REAL(temp(:,index_obs(i)),8)
           END DO
   
           ! Temperature profile quality check
           DO i = 1, nprof_day
              pro_temp_qc_day(i) = pro_temp_qc(index_obs(i))
           END DO
   
           ! level-temperature quality check
           DO i = 1, nprof_day
              temp_qc_day(i,:) = temp_qc(:,index_obs(i))
           END DO
   
   
           ! *** Salinity ***
   
           ALLOCATE(sal_day(nprof_day,n_levels))
           ALLOCATE(pro_sal_qc_day(nprof_day))
           ALLOCATE(sal_qc_day(nprof_day,n_levels))
   
           ! Salinity values
           DO i = 1, nprof_day
              sal_day(i,:) = REAL(sal(:,index_obs(i)),8)
           END DO
           
           ! Salinity profile quality check
           DO i = 1, nprof_day
              pro_sal_qc_day(i) = pro_sal_qc(index_obs(i))
           END DO
   
           ! level-salinity quality check information
           DO i = 1, nprof_day
              sal_qc_day(i,:) = sal_qc(:,index_obs(i))
           END DO
   
   
           ! ******************************************************
           ! *** Find profiles at current day that are PE-local ***
           ! ******************************************************
           
           ! Check and find out if these profiles belong to the sub-domain: 
           ! If yes, store the index and the corresponding element indices.
           ALLOCATE(index_ele_day_p(nprof_day))
           ALLOCATE(index_pro_day_p(nprof_day))
   
           nprof_day_p = 0
           DO iprof_day = 1, nprof_day
              CALL point_in_triangle_pdaf(element,obs_coord_day(:,iprof_day))
              IF (element /= 0) THEN
                 nprof_day_p = nprof_day_p +1
                 index_pro_day_p(nprof_day_p) = iprof_day 
                 index_ele_day_p(nprof_day_p) = element
              END IF
           END DO
   
   
           ! **************************************************************
           ! *** Quality check for PE-local profiles of the current day ***
           ! **************************************************************
           
           cnt_prof_temp = 0
           cnt_prof_sal  = 0
           cnt_temp      = 0
           cnt_sal       = 0
           cnt_en2d_sal  = 0
           cnt_flagged_temp_p = 0
           cnt_flagged_sal_p  = 0
   
           IF (nprof_day_p > 0) THEN
   
              ! Allocate arrays for observations passing quality checks
   
              ALLOCATE(n2d_temp           (3, nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(n2d_sal            (3, nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(nl1_temp           (   nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(nl1_sal            (   nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(depth_temp         (   nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(depth_sal          (   nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(temp_obs           (   nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(sal_obs            (   nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(obs_coord_day_temp (2, nprof_day_p* (mesh_fesom% nl-1)))
              ALLOCATE(obs_coord_day_sal  (2, nprof_day_p* (mesh_fesom% nl-1)))
              
              ! Filling up from 1 to cnt_temp:
              
              ! n2d_temp:            Index refering to model node
              ! nl1_temp:            Index refering to model layer
              ! depth_temp:          Depth of the model layer
              ! temp_obs:            Quality-checked temperature averaged on model depth
              ! obs_coord_day_temp:  Observation coordinate
   
              LOOP_profile: DO iprof_day_p = 1, nprof_day_p
              
                 ! Debugging:
                 
                 IF ( (obs_coord_day(1,index_pro_day_p(iprof_day_p)) > -69.3  *rad) .AND. &
                      (obs_coord_day(1,index_pro_day_p(iprof_day_p)) < -69.2  *rad) .AND. &
                      (obs_coord_day(2,index_pro_day_p(iprof_day_p)) >  34.112*rad) .AND. &
                      (obs_coord_day(2,index_pro_day_p(iprof_day_p)) <  34.14 *rad))        THEN
                      
                      debugging = .TRUE.
                 ELSE
                 
                      debugging = .FALSE.
                 END IF
                 
				 !IF (debugging) WRITE(*,*) 'proffile-debug', 'Z_n ', Z_n
			     !IF (debugging) WRITE(*,*) 'proffile-debug', 'zbar_n ', zbar_n
				 IF (debugging) WRITE(*,*) 'proffile-debug ', 'depth_day(index_pro_day_p(iprof_day_p),:)', depth_day(index_pro_day_p(iprof_day_p),:)

                 
                 IF (debugging) WRITE(*,*) 'proffile-debug ', 'iprof_day_p', iprof_day_p
                 IF (debugging) WRITE(*,*) 'proffile-debug ', 'obs_coord_day(1,index_pro_day_p(iprof_day_p))/rad', obs_coord_day(1,index_pro_day_p(iprof_day_p))/rad
                 IF (debugging) WRITE(*,*) 'proffile-debug ', 'obs_coord_day(2,index_pro_day_p(iprof_day_p))/rad', obs_coord_day(2,index_pro_day_p(iprof_day_p))/rad
   
                 ! *******************
                 ! *** Temperature ***
                 ! *******************
   
                 ! Profile check
                 IF ((pos_qc_day     (index_pro_day_p(iprof_day_p)) == one) .AND. &
                     (pro_temp_qc_day(index_pro_day_p(iprof_day_p)) == one))      THEN
                     
                    ! Find nodes
                    flag_n2d = 0
                    n2d_prof(:) = mesh_fesom% elem2D_nodes(:,index_ele_day_p(iprof_day_p))
                    DO s = 1, 3
                    ! Reject profiles on halo
                       IF (n2d_prof(s) > myDim_nod2D) THEN
                          flag_n2d = 1
                          cnt_flagged_temp_p = cnt_flagged_temp_p + 1
                       ENDIF
                    END DO
   
                    ! Number of model layers at profiles (considering topography)
                    IF (flag_n2d == 0) THEN
                       nlayer = MIN( mesh_fesom% nlevels_nod2D(n2d_prof(1))-1, &
                                     mesh_fesom% nlevels_nod2D(n2d_prof(2))-1, &
                                     mesh_fesom% nlevels_nod2D(n2d_prof(3))-1)
                    ELSE
                       nlayer = 0
                    END IF
                  
                    nlayer= min(nlayer,mesh_fesom%nl-1)
                    
                    IF (debugging) WRITE(*,*) 'proffile-debug ', 'mesh_fesom% nlevels_nod2D(n2d_prof(1))', mesh_fesom% nlevels_nod2D(n2d_prof(1))
                    IF (debugging) WRITE(*,*) 'proffile-debug ', 'mesh_fesom% nlevels_nod2D(n2d_prof(2))', mesh_fesom% nlevels_nod2D(n2d_prof(2))
                    IF (debugging) WRITE(*,*) 'proffile-debug ', 'mesh_fesom% nlevels_nod2D(n2d_prof(3))', mesh_fesom% nlevels_nod2D(n2d_prof(3))
                    IF (debugging) WRITE(*,*) 'proffile-debug ', 'nlayer', nlayer

                    IF (nlayer > 0) THEN
   
                       cnt_prof_temp = cnt_prof_temp + 1
   
                       DO ilayer = 1, nlayer
                          cnt_temp = cnt_temp + 1
                          
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'cnt_temp', cnt_temp
                          
                          DO s = 1, 3 
                                n2d_temp (s,cnt_temp) = n2d_prof(s)
                          END DO
                          
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'n2d_temp (1,cnt_temp)', n2d_temp (1,cnt_temp)
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'n2d_temp (2,cnt_temp)', n2d_temp (2,cnt_temp)
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'n2d_temp (3,cnt_temp)', n2d_temp (3,cnt_temp)
   
                          ! Coordinate of the observation
                          obs_coord_day_temp(:, cnt_temp) = obs_coord_day(:,index_pro_day_p(iprof_day_p))
                          
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'obs_coord_day_temp(1, cnt_temp)/rad', obs_coord_day_temp(1, cnt_temp)/rad
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'obs_coord_day_temp(2, cnt_temp)/rad', obs_coord_day_temp(2, cnt_temp)/rad
                          
                          
                          ! Model layer
                          nl1_temp (cnt_temp) = ilayer
                          
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'nl1_temp (cnt_temp)', nl1_temp (cnt_temp)
                            
                          ! Depth of the model layer (mid between upper and lower layer bounds)
                          depth_temp(cnt_temp) = mesh_fesom% Z(ilayer)
                          
                          IF (debugging) WRITE(*,*) 'proffile-debug ', 'depth_temp(cnt_temp)', depth_temp(cnt_temp)

                          ! Level quality check
                          cnt_cal  = 0
                          sum_temp = 0.0
                          DO nn = 1, n_levels
                             ! Check if nn'th observation depth is on model layer (ilayer)
                             IF (temp_qc_day(index_pro_day_p(iprof_day_p),nn) == one) THEN

                                IF ((ABS(depth_day(index_pro_day_p(iprof_day_p),nn)) >  ABS(mesh_fesom% zbar(ilayer  )))  .and.  &
                                    (ABS(depth_day(index_pro_day_p(iprof_day_p),nn)) <= ABS(mesh_fesom% zbar(ilayer+1)))         ) THEN
                                   cnt_cal  = cnt_cal + 1
                                   sum_temp = sum_temp + temp_day(index_pro_day_p(iprof_day_p),nn)
                                   
                                   IF (debugging) WRITE(*,*) 'proffile-debug ', 'depth_day(index_pro_day_p(iprof_day_p),nn)', depth_day(index_pro_day_p(iprof_day_p),nn)
                                   IF (debugging) WRITE(*,*) 'proffile-debug ', 'zbar_n(ilayer  )', zbar_n(ilayer  )
                                   IF (debugging) WRITE(*,*) 'proffile-debug ', 'zbar_n(ilayer+1)', zbar_n(ilayer+1)
                                   
                                END IF
                             END IF
                          END DO
                          ! Mean of nn observations on model layer (ilayer)
                          IF (cnt_cal /= 0) THEN
                             temp_obs(cnt_temp) = sum_temp / cnt_cal
                             
                             IF (debugging) WRITE(*,*) 'proffile-debug ', 'temp_obs(cnt_temp)', temp_obs(cnt_temp) 
                             
                          ELSE
                             cnt_temp = cnt_temp - 1
                          END IF
                       END DO !  DO ilayer = 1, nlayer
                    END IF ! IF (nlayer > 0)
                 END IF !  IF ((pos_qc_day(index_pro_day_p(iprof_day_p)) == one) ...
   
   
                 ! ****************
                 ! *** Salinity ***
                 ! ****************
   
                 ! profile check
                 IF ((pos_qc_day(index_pro_day_p(iprof_day_p)) == one) &
                      .AND. (pro_sal_qc_day(index_pro_day_p(iprof_day_p)) == one)) THEN
   
                    flag_n2d = 0
                    n2d_prof(:) = mesh_fesom%elem2D_nodes(:,index_ele_day_p(iprof_day_p))
                    DO s = 1, 3
                       IF (n2d_prof(s) > myDim_nod2D) then
                          flag_n2d = 1
                          cnt_flagged_sal_p = cnt_flagged_sal_p + 1
                       ENDIF
                    END DO

                    IF (flag_n2d == 0) THEN
                       nlayer = MIN(mesh_fesom%nlevels_nod2D(n2d_prof(1))-1, &
                                    mesh_fesom%nlevels_nod2D(n2d_prof(2))-1, &
                                    mesh_fesom%nlevels_nod2D(n2d_prof(3))-1)
                    ELSE
                       nlayer = 0
                       cnt_en2d_sal = cnt_en2d_sal + 1
                    END IF

                    nlayer= min(nlayer,mesh_fesom%nl-1)

                    IF (nlayer > 0) THEN
   
                       cnt_prof_sal = cnt_prof_sal + 1
   
                       DO ilayer = 1, nlayer
                          cnt_sal = cnt_sal + 1
                          DO s = 1, 3
                               n2d_sal (s, cnt_sal) = n2d_prof(s)
                          END DO
                          nl1_sal (cnt_sal) = ilayer
   
                          ! Coordinate of the observations
                          obs_coord_day_sal(:, cnt_sal) = obs_coord_day(:,index_pro_day_p(iprof_day_p))
   
                          ! depth of the qq-th layer
                          depth_sal (cnt_sal) = mesh_fesom% Z(ilayer)
   
                          ! level quality check
                          cnt_cal = 0
                          sum_sal = 0.0
                          DO nn = 1, n_levels
                             IF (sal_qc_day(index_pro_day_p(iprof_day_p),nn) == one) THEN
                              IF ((ABS(depth_day(index_pro_day_p(iprof_day_p),nn)) >  ABS(mesh_fesom% zbar(ilayer))) .and.  &
                              (ABS(depth_day(index_pro_day_p(iprof_day_p),nn))     <= ABS(mesh_fesom% zbar(ilayer+1)))) THEN
                                   cnt_cal = cnt_cal + 1
                                   sum_sal = sum_sal + sal_day(index_pro_day_p(iprof_day_p),nn)
                                END IF
                             END IF
                          END DO
                          IF (cnt_cal /= 0) THEN
                             sal_obs(cnt_sal) = sum_sal / cnt_cal
                          ELSE
                             cnt_sal = cnt_sal - 1
                          END IF
                       END DO ! DO ilayer = 1, nlayer
                    END IF ! IF (nlayer > 0)
                 END IF ! IF ((pos_qc_day(index_pro_day_p(iprof_day_p)) == a)
              END DO LOOP_profile
   
           END IF ! IF (nprof_day_p > 0)
   
   
           ! *** Get total number of profiles for temperature and salt ***
   
           ! Global number of temperature profiles
           CALL MPI_Allreduce (cnt_prof_temp, cnt_pro_temp_g, 1, MPI_INTEGER, MPI_SUM, &
                COMM_filter, MPIerr)
   
           ! Global number of salinity profiles
           CALL MPI_Allreduce (cnt_prof_sal, cnt_pro_sal_g, 1, MPI_INTEGER, MPI_SUM, &
                COMM_filter, MPIerr)
   
           ! Global number of temperature observations
           CALL MPI_Allreduce (cnt_temp, cnt_temp_g, 1, MPI_INTEGER, MPI_SUM, &
                COMM_filter, MPIerr)
   
           ! Global number of salinity observations
           CALL MPI_Allreduce (cnt_sal, cnt_sal_g, 1, MPI_INTEGER, MPI_SUM, &
                COMM_filter, MPIerr)
   
           ! Global number of flagged profiles
           CALL MPI_Allreduce (cnt_flagged_temp_p, cnt_flagged_temp, 1, MPI_INTEGER, MPI_SUM, &
                COMM_filter, MPIerr)
           CALL MPI_Allreduce (cnt_flagged_sal_p, cnt_flagged_sal, 1, MPI_INTEGER, MPI_SUM, &
                COMM_filter, MPIerr)
   
           ! Write number of profiles
           IF (mype_filter==0) THEN
              WRITE (*,'(a,5x,a,i2,a,i2,x,a,2i7,x,a,2i7)') &
                   'FESOM-PDAF', '--- Date ', d_in_month,'.', month_index, &
                   'number of profiles: T, S', cnt_pro_temp_g, cnt_pro_sal_g, &
                   'observations:', cnt_temp_g, cnt_sal_g
   
              WRITE (*,'(a,5x,a,i2,a,i2,x,a,2i7)') &
                   'FESOM-PDAF', '--- Date ', d_in_month,'.', month_index, &
                   'number of flagged nodes:', cnt_flagged_temp, cnt_flagged_sal
           END IF
   
   
   
           ! ****************************************
           ! *** Write variables into netCDF file ***
           ! ****************************************
   
           ! Write 1-D variable   
           
           startv = cnt_day
           countv = 1
   
           s = 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_time, startv, countv, cnt_day)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_nprof_temp, startv, countv, cnt_prof_temp)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_nprof_sal, startv, countv, cnt_prof_sal)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_nflagged_temp, startv, countv, cnt_flagged_temp_p)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_nflagged_sal, startv, countv, cnt_flagged_sal_p)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_nobs_temp, startv, countv, cnt_temp)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_nobs_sal, startv, countv, cnt_sal)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_offset_temp, startv, countv, offset_temp)
           s = s + 1
           stat(s) = NF_PUT_VARA_INT(ncid_out, id_offset_sal, startv, countv, offset_sal)  
   
           startv = offset_temp
           countv = cnt_temp
   
           IF (cnt_temp > 0) THEN
              s = s + 1
              stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_depth_temp, startv, countv, depth_temp(1:cnt_temp))
              s = s + 1
              stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_temp_out, startv, countv, temp_obs(1:cnt_temp))
              s = s + 1
              stat(s) = NF_PUT_VARA_INT(ncid_out, id_layer_temp, startv, countv, nl1_temp(1:cnt_temp))
           END IF
   
           startv = offset_sal
           countv = cnt_sal
   
           IF (cnt_sal >0) THEN
              s = s + 1
              stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_depth_sal, startv, countv, depth_sal(1:cnt_sal))
              s = s + 1
              stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_sal_out, startv, countv, sal_obs(1:cnt_sal))
              s = s + 1
              stat(s) = NF_PUT_VARA_INT(ncid_out, id_layer_sal, startv, countv, nl1_sal(1:cnt_sal))
           END IF
   
           DO i = 1, s
              IF (stat(i) /= NF_NOERR) &
                   WRITE(*, *) 'NetCDF error in writing profile variable, no.', i
           END DO
   
           ! Write 2-D variables
   
           startv2(1) = 1
           countv2(1) = 3       
           startv2(2) = offset_temp
           countv2(2) = cnt_temp
   
           s = 0
           IF (cnt_temp > 0) THEN
              s = s + 1
              stat(s) = NF_PUT_VARA_INT(ncid_out, id_node_temp, startv2, countv2, n2d_temp(:,1:cnt_temp))
   
              countv2(1) = 2
   
              s = s + 1
              stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_coord_temp, startv2, countv2, obs_coord_day_temp(:,1:cnt_temp))
   
              if (stat(s)/=NF_NOERR) write (*,*) NF_STRERROR(stat(s))
           END IF
   
           startv2(1) = 1
           countv2(1) = 3       
           startv2(2) = offset_sal
           countv2(2) = cnt_sal
   
           IF (cnt_sal >0) THEN
              s = s + 1
              stat(s) = NF_PUT_VARA_INT(ncid_out, id_node_sal, startv2, countv2, n2d_sal(:,1:cnt_sal))
   
              countv2(1) = 2
   
              s = s + 1
              stat(s) = NF_PUT_VARA_DOUBLE(ncid_out, id_coord_sal, startv2, countv2, obs_coord_day_sal(:,1:cnt_sal))
           END IF
      
           DO i = 1, s
              IF (stat(i) /= NF_NOERR) &
                   WRITE(*, *) 'NetCDF error in writing profile variables, no.', i, mype_filter
           END DO
   
           cnt_day = cnt_day + 1
           offset_temp = offset_temp + cnt_temp
           offset_sal = offset_sal + cnt_sal
   
   
           ! **********************
           ! *** Daily clean up ***
           ! **********************
   
           DEALLOCATE(obs_coord_day, depth_day, pos_qc_day)
           DEALLOCATE(temp_day, pro_temp_qc_day, temp_qc_day)
           DEALLOCATE(sal_day, pro_sal_qc_day, sal_qc_day)
           DEALLOCATE(index_pro_day_p,index_ele_day_p)
   
           IF (nprof_day_p > 0) THEN
              DEALLOCATE(obs_coord_day_temp,n2d_temp,nl1_temp)
              DEALLOCATE(obs_coord_day_sal,n2d_sal,nl1_sal)
              DEALLOCATE(depth_temp, depth_sal, temp_obs, sal_obs)
           END IF
     
        END DO LOOP_day
   
        ! Deallocate monthly arrays
   
        DEALLOCATE(obs_coord, pos_qc, depth)
        DEALLOCATE(temp, pro_temp_qc, temp_qc)
        DEALLOCATE(sal, pro_sal_qc, sal_qc)
        DEALLOCATE(judate)
        DEALLOCATE(index_obs)
   
     END DO LOOP_month
   
   
   ! ********************
   ! *** Finishing up ***
   ! ********************
   
     s = s + 1
     stat(s) = NF_CLOSE(ncid_out)
   
   END DO
   CALL MPI_BARRIER(COMM_filter,MPIerr)
   CALL abort_parallel
   END SUBROUTINE init_dim_obs_f_proffile_pdaf
   
