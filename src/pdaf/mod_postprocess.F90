MODULE mod_postprocess

use mod_assim_pdaf, &
    only: nlmax, id, mesh_fesom, dim_state_p, offset, &
    DAoutput_path, dim_fields
use mod_parallel_pdaf, &
    only: writepe, COMM_filter, MPIerr, mype_filter
use mod_nc_out_variables
use mod_nc_out_routines, &
    only: check
use obs_o2_comf_pdafomi, &
    only: assim_o_o2_comf
use obs_o2_argo_pdafomi, &
    only: assim_o_o2_argo
use obs_o2_merged_pdafomi, &
    only: assim_o_o2_merged
use obs_n_merged_pdafomi, &
    only: assim_o_n_merged

use g_config, &
    only: dt
use g_clock, &
    only: clock, yearnew, daynew, timenew, cyearnew, day_in_month, &
    month
use g_parsup, &
    only: myDim_nod2D, eDim_nod2D
use g_comm_auto, &
    only: broadcast_nod

use netcdf

IMPLICIT NONE
SAVE

include 'mpif.h'

! postprocessing settings
logical :: isPP = .false.                   ! .true. - run postprocessing code instead of model simulation
integer :: yearPP             
integer, parameter :: FRRN=1                ! free run
integer, parameter :: ASML=2                ! assimilation run

! PDAF simulation output
INTEGER :: fidphyday(2), fidbgcday(2), &    ! IDs of netCDF files
           fidphymon(2), fidbgcmon(2)
integer :: dimid_time
CHARACTER(len=200) :: pathsim(2)            ! path
CHARACTER(len=200) :: fnameread             ! filename
integer :: ndaysim_tmp(2), ndaysim          ! number of days
real, allocatable :: state_p(:,:)           ! model state
real, allocatable :: forc_p(:)              ! forecast ASML

! model observation comparison
INTEGER              :: dim_obs_f        ! number of observations
INTEGER              :: stats            ! status flag for PDAF gather operation
REAL,    ALLOCATABLE :: ostate_f(:,:)    ! observed model state
REAL,    ALLOCATABLE :: diff_f(:,:)      ! model-observation difference
REAL,    ALLOCATABLE :: impr_f(:)        ! improvement
REAL,    ALLOCATABLE :: diffAF_f(:)      ! difference ASML-FREE
REAL,    ALLOCATABLE :: oforc_f(:)       ! observed forecast ASML

! observation types
integer, parameter :: iDIC=1, iAlk=2, ipCO2s=3, iO2comf=4, iNcomf=5, &
                      iOargo=6, iNargo=7, iOmerged=8, iNmerged=9
integer, parameter :: nobs=9         ! number of observation types
character(len=20)  :: cobstype(nobs) ! filename
integer            :: fidobs(nobs)   ! file ID

! model-observation difference exclusion
real :: o2_merged_excl_absolutePP, o2_merged_excl_relativePP
real :: n_merged_excl_absolutePP,  n_merged_excl_relativePP

! output
integer, allocatable :: offset_day(:,:)       ! specifies where to write
REAL, parameter      :: fill_value = -999.0
INTEGER, PARAMETER   :: to_nf90_byte = SELECTED_INT_KIND(2)

! general
integer :: iday, isim, ifield, i, n, nz, s, d  ! counters


CONTAINS

! ************
! *** MAIN ***
! ************

SUBROUTINE doPP(nsteps)

   implicit none

   ! arguments
   integer, intent(inout) :: nsteps ! number of model time steps
   
   if (writepe) then
        write(*,*) '===================================================='
        write(*,*) '========= PDAF POSTPROCESSING OBSERVATIONS ========='
        write(*,*) '===================================================='
   endif
   
   ! init
   cobstype(iDIC    ) = 'DIC_GLODAP'
   cobstype(iAlk    ) = 'ALK_GLODAP'
   cobstype(iO2comf ) = 'OXY_COMFORT'
   cobstype(ipCO2s  ) = 'pCO2_SOCAT'
   cobstype(iNcomf  ) = 'DIN_COMFORT'
   cobstype(iOargo  ) = 'OXY_ARGO'
   cobstype(iNargo  ) = 'DIN_ARGO'
   cobstype(iOmerged) = 'OXY_MERGED'
   cobstype(iNmerged) = 'DIN_MERGED'
   
   ! init time step and clock
   nsteps = 0 ! set model steps to be run in main model loop to zero
   dt = 86400 ! set time step to daily to process daily observations
   yearnew = yearPP
   daynew  = 1
   timenew = 0
   write(cyearnew,'(i4)') yearnew
   call clock
   
   ! reset observation exclusion criteria
   call reset_exclusion_limits
   
   ! open model file
   if (writepe) then
      do isim=1,2
         call netCDF_openfile('bgc','day',isim,fidbgcday(isim))
      enddo
   endif
   
   ! length of simulation
   if (writepe) then
      do isim=1,2
         call check(nf90_inq_dimid(fidbgcday(isim),'time', dimid_time))
         call check(nf90_inquire_dimension(fidbgcday(isim), dimid_time, len=ndaysim_tmp(isim)))
      enddo
   ndaysim=MIN(ndaysim_tmp(FRRN),ndaysim_tmp(ASML))
   endif
   call MPI_Bcast(ndaysim,1, MPI_INTEGER, 0, COMM_filter, MPIerr)
   allocate(offset_day(nobs,ndaysim))
   offset_day(:,1)=1
   
   ! create output files
   if (writepe) then
      call create_ncfile(cobstype(iDIC    ), fidobs(iDIC    ))
      call create_ncfile(cobstype(iAlk    ), fidobs(iAlk    ))
      ! removed: call create_ncfile(cobstype(iO2comf ), fidobs(iO2comf ))
      call create_ncfile(cobstype(ipCO2s  ), fidobs(ipCO2s  ))
      ! removed: call create_ncfile(cobstype(iNcomf  ), fidobs(iNcomf  ))
      ! removed: call create_ncfile(cobstype(iOargo  ), fidobs(iOargo  ))
      ! removed: call create_ncfile(cobstype(iNargo  ), fidobs(iNargo  ))
      call create_ncfile(cobstype(iOmerged), fidobs(iOmerged))
      call create_ncfile(cobstype(iNmerged), fidobs(iNmerged))
   endif
   
   ! init simulation data
   allocate(state_p(2,dim_state_p))
   allocate(forc_p(dim_state_p))
   ! removed: IF (assim_o_o2_comf .or. assim_o_o2_merged .or. assim_o_o2_argo .or. &
   ! removed:    assim_o_n_comf  .or. assim_o_n_merged  .or. assim_o_n_argo )     &
   ! removed:    allocate(forc_p(dim_state_p))
   
   ! daily time loop
   do iday=1,ndaysim
      ! message
      if (writepe) WRITE (*,'(a,5x,a,1x,i2,a1,i2,a1,i4,1x,f5.2,a)') &
            'FESOM-PDAF', 'Postprocessing at', &
            day_in_month, '.', month, '.', yearnew, timenew/3600.0, 'h'
      ! read simulation data
      call netCDF_getstate(FRRN,daynew)
      call netCDF_getstate(ASML,daynew)
      ! use for InnoOmit / observation exclusion criteria; strictly, these are not forecast data.
      forc_p(:)=state_p(ASML,:)
      ! removed: read forecast for InnoOmit / observation exclusion criteria
      ! removed: call netCDF_getforc (ASML,daynew)
      
      ! call observation modules
      call PP_DIC_GLODAP()
      call PP_Alk_GLODAP()
      ! call PP_O2_COMFORT()
      call PP_PCO2_SOCAT()
      ! call PP_DIN_COMFORT()
      ! call PP_O2_ARGO()
      ! call PP_DIN_ARGO()
      call PP_O2_MERGED()
      call PP_DIN_MERGED()
      
      call clock
   enddo
   
   ! clean up
   if (writepe) then
      ! close simulation files
      do isim=1,2
         call check(NF90_close(fidbgcday(isim)))
      enddo
      ! close output files
      call check(nf90_close(fidobs(iDIC    )))
      call check(nf90_close(fidobs(iAlk    )))
      ! call check(nf90_close(fidobs(iO2comf )))
      call check(nf90_close(fidobs(ipCO2s  )))
      ! call check(nf90_close(fidobs(iNcomf  )))
      ! call check(nf90_close(fidobs(iNargo  )))
      ! call check(nf90_close(fidobs(iOargo  )))
      call check(nf90_close(fidobs(iOmerged)))
      call check(nf90_close(fidobs(iNmerged)))
   endif
   deallocate(state_p)
   IF (allocated(forc_p)) deallocate(forc_p)

END SUBROUTINE doPP


! ******************************
! *** reset exclusion limits ***
! ******************************
! observations which are excluded in the observation model due to model-observation differences,
! nevertheless should, for evaluation purposes, be included in the postprocessing

SUBROUTINE reset_exclusion_limits()
   use obs_o2_merged_pdafomi, &
       only: o2_merged_excl_absolute, o2_merged_excl_relative
   use obs_n_merged_pdafomi, &
       only: n_merged_excl_absolute, n_merged_excl_relative
       
   save

   ! save original limits
   o2_merged_excl_absolutePP  = o2_merged_excl_absolute
   o2_merged_excl_relativePP  = o2_merged_excl_relative
   n_merged_excl_absolutePP   = n_merged_excl_absolute
   n_merged_excl_relativePP   = n_merged_excl_relative
   
   ! reset limits used in observation module to zero
   o2_merged_excl_absolute = 0
   o2_merged_excl_relative = 0
   n_merged_excl_absolute  = 0
   n_merged_excl_relative  = 0
   
END SUBROUTINE reset_exclusion_limits


! ***********************
! *** netCDF_openfile ***
! ***********************
! simulation data

SUBROUTINE netCDF_openfile(typ,freq,sim,fid)
   ! Arguments
   CHARACTER(len=3), intent (in)  :: typ   ! phy or bgc
   CHARACTER(len=3), intent (in)  :: freq  ! day or mon
   INTEGER,          intent (in)  :: sim   ! FRRN or ASML
   INTEGER,          intent (out) :: fid   ! file ID
          
   fnameread = TRIM(pathsim(sim))//'fesom-recom-pdaf.'//typ//'.'//cyearnew//'.'//freq//'.nc'
   WRITE (*,'(a,5x,a,1x,a)') &
            'FESOM-PDAF', 'Read simulation output from file', &
             fnameread
   ! open file
   call check(NF90_OPEN(trim(fnameread),NF90_NOWRITE,fid))
   
END SUBROUTINE netCDF_openfile


! ***********************
! *** netCDF_getstate ***
! ***********************
! simulation data

SUBROUTINE netCDF_getstate(sim,ReadAtDay)
   ! Arguments
   INTEGER,      intent (in)  :: sim                ! FRRN or ASML
   INTEGER,      intent (in)  :: ReadAtDay
   ! Local variables:
   REAL, allocatable          :: myData2(:)         ! Temporary array for pe-local surface fields
   REAL, allocatable          :: myData3(:,:)       ! Temporary array for pe-local 3D-fields
   REAL(kind=4), allocatable  :: data2_g(:)         ! Temporary array for global surface fields
   REAL(kind=4), allocatable  :: data3_g(:,:)       ! Temporary array for global 3D-fields
   INTEGER, parameter         :: nfields_obs = 7
   INTEGER                    :: observedfields(nfields_obs)
   INTEGER                    :: varid
   
   observedfields = (/ id% PhyChl,  &
                       id% DiaChl,  &
                       id% DIC,     &   
                       id% Alk,     &   
                       id% pCO2s,   & 
                       id% O2,      &    
                       id% DIN      &
                     /)   
   
   do i=1,nfields_obs
      ifield=observedfields(i)
      if (sfields(ifield)%ndims==2) then
         ! 3D fields
         allocate(data3_g(mesh_fesom%nod2D,nlmax))
         allocate(myData3(nlmax, myDim_nod2D+eDim_nod2D))
         ! read global state
         if (writepe) then
            call check(nf90_inq_varid(fidbgcday(sim), TRIM(sfields(ifield)%variable)//'_mm', varid))
            call check(nf90_get_var(fidbgcday(sim),varid,data3_g,                &
                                    start=(/                1,     1, ReadAtDay /), &
                                    count=(/ mesh_fesom%nod2D, nlmax,      1    /) ))
         endif
         ! broadcast state
         call broadcast_nod(myData3,REAL(TRANSPOSE(data3_g),8))
         ! add to state vector
         do n=1, myDim_nod2D
            do nz=1, nlmax
               s = (n-1) * (nlmax) + nz
               state_p(sim,offset(ifield)+s) = myData3(nz,n)
            enddo ! nz
         enddo ! n
         deallocate(data3_g,myData3)
      else
         ! surface fields
         allocate(data2_g(mesh_fesom%nod2D))
         allocate(myData2(myDim_nod2D+eDim_nod2D))
         if (writepe) then
            call check(nf90_inq_varid(fidbgcday(sim), TRIM(sfields(ifield)%variable)//'_mm', varid))
            call check(nf90_get_var(fidbgcday(sim),varid,data2_g,          &
                                    start=(/                1, ReadAtDay /), &
                                    count=(/ mesh_fesom%nod2D,         1 /) ))
         endif
         ! broadcast state
         call broadcast_nod(myData2,REAL(data2_g,8))
         ! add to state vector
         do n=1, myDim_nod2D
               state_p(sim,offset(ifield)+n) = myData2(n)
         enddo ! n
         deallocate(data2_g,myData2)
      endif
   enddo ! i, ifield
END SUBROUTINE netCDF_getstate


! ***********************
! *** netCDF_getforc  ***
! ***********************
! simulation data

SUBROUTINE netCDF_getforc(sim,ReadAtDay)
   ! Arguments
   INTEGER,      intent (in)  :: sim                ! FRRN or ASML
   INTEGER,      intent (in)  :: ReadAtDay
   ! Local variables:
   REAL, allocatable          :: myData2(:)         ! Temporary array for pe-local surface fields
   REAL, allocatable          :: myData3(:,:)       ! Temporary array for pe-local 3D-fields
   REAL(kind=4), allocatable  :: data2_g(:)         ! Temporary array for global surface fields
   REAL(kind=4), allocatable  :: data3_g(:,:)       ! Temporary array for global 3D-fields
   INTEGER                    :: nfields_obs
   INTEGER, allocatable       :: observedfields(:)
   INTEGER                    :: varid
   
   
   IF ((assim_o_o2_comf .or. assim_o_o2_merged .or. assim_o_o2_argo) .and. &
       (assim_o_n_comf  .or. assim_o_n_merged  .or. assim_o_n_argo ))      &
       then
          nfields_obs = 2
          allocate(observedfields(nfields_obs))
          observedfields = (/ id% O2, id% DIN /)
   ELSEIF ((assim_o_o2_comf .or. assim_o_o2_merged .or. assim_o_o2_argo))  &
       then
          nfields_obs = 1
          allocate(observedfields(nfields_obs))
          observedfields = (/ id% O2 /)
   ELSEIF ((assim_o_n_comf .or. assim_o_n_merged .or. assim_o_n_argo))  &
       then
          nfields_obs = 1
          allocate(observedfields(nfields_obs))
          observedfields = (/ id% DIN /)
   ELSE
       nfields_obs = 0
       allocate(observedfields(nfields_obs))
   ENDIF
   
   do i=1,nfields_obs
      ifield=observedfields(i)
      if (sfields(ifield)%ndims==2) then
         ! 3D fields
         allocate(data3_g(mesh_fesom%nod2D,nlmax))
         allocate(myData3(nlmax, myDim_nod2D+eDim_nod2D))
         ! read global state
         if (writepe) then
            call check(nf90_inq_varid(fidbgcday(sim), TRIM(sfields(ifield)%variable)//'_mm', varid))
            call check(nf90_get_var(fidbgcday(sim),varid,data3_g,                &
                                    start=(/                1,     1, ReadAtDay /), &
                                    count=(/ mesh_fesom%nod2D, nlmax,      1    /) ))
         endif
         ! broadcast state
         call broadcast_nod(myData3,REAL(TRANSPOSE(data3_g),8))
         ! add to state vector
         do n=1, myDim_nod2D
            do nz=1, nlmax
               s = (n-1) * (nlmax) + nz
               forc_p(offset(ifield)+s) = myData3(nz,n)
            enddo ! nz
         enddo ! n
         deallocate(data3_g,myData3)
      else
         ! surface fields
         allocate(data2_g(mesh_fesom%nod2D))
         allocate(myData2(myDim_nod2D+eDim_nod2D))
         if (writepe) then
            call check(nf90_inq_varid(fidbgcday(sim), TRIM(sfields(ifield)%variable)//'_mm', varid))
            call check(nf90_get_var(fidbgcday(sim),varid,data2_g,          &
                                    start=(/                1, ReadAtDay /), &
                                    count=(/ mesh_fesom%nod2D,         1 /) ))
         endif
         ! broadcast state
         call broadcast_nod(myData2,REAL(data2_g,8))
         ! add to state vector
         do n=1, myDim_nod2D
               forc_p(offset(ifield)+n) = myData2(n)
         enddo ! n
         deallocate(data2_g,myData2)
      endif
   enddo ! i, ifield
   
   deallocate(observedfields)
   
END SUBROUTINE netCDF_getforc


! *********************
! *** PP_DIC_GLODAP ***
! *********************
SUBROUTINE PP_DIC_GLODAP()

   use obs_DIC_glodap_pdafomi, &
       only: init_dim_obs_DIC_glodap, &
       obs_op_DIC_glodap, DIC_glodap_exclude_diff, &
       thisobs, thisobs_PP, thisobs_PP_f
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_DIC_glodap(-1, dim_obs_f)
   ! offset
   offset_day(iDIC,iday+1)=offset_day(iDIC,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing GLODAP-DIC: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_DIC_glodap(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      if (DIC_glodap_exclude_diff > 0.0) then
         allocate(oforc_f(thisobs%dim_obs_f))
         call obs_op_DIC_glodap(dim_state_p, thisobs%dim_obs_f, forc_p, oforc_f)
         WHERE (ABS(oforc_f-thisobs%obs_f)>DIC_glodap_exclude_diff) &
            thisobs_PP_f%isInnoOmit = 1.0
      endif
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iDIC),offset_day(iDIC,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f, oforc_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iDIC),offset_day(iDIC,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth)
   
END SUBROUTINE PP_DIC_GLODAP


! *********************
! *** PP_ALK_GLODAP ***
! *********************
SUBROUTINE PP_Alk_GLODAP()

   use obs_Alk_glodap_pdafomi, &
       only: init_dim_obs_Alk_glodap, &
       obs_op_Alk_glodap, Alk_glodap_exclude_diff, &
       thisobs, thisobs_PP, thisobs_PP_f
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_Alk_glodap(-1, dim_obs_f)
   ! offset
   offset_day(iAlk,iday+1)=offset_day(iAlk,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing GLODAP-Alk: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_Alk_glodap(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      if (Alk_glodap_exclude_diff > 0.0) then
            allocate(oforc_f(thisobs%dim_obs_f))
            call obs_op_Alk_glodap(dim_state_p, thisobs%dim_obs_f, forc_p, oforc_f)
            WHERE (ABS(oforc_f-thisobs%obs_f)>Alk_glodap_exclude_diff) &
               thisobs_PP_f%isInnoOmit = 1.0
      endif
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iAlk),offset_day(iAlk,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f, oforc_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iAlk),offset_day(iAlk,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth)
   
END SUBROUTINE PP_Alk_GLODAP


! *********************
! *** PP_O2_COMFORT ***
! *********************
SUBROUTINE PP_O2_COMFORT()

   use obs_O2_COMF_pdafomi, &
       only: init_dim_obs_O2_COMF, &
       obs_op_O2_COMF, &
       thisobs, thisobs_PP, thisobs_PP_f, &
       o2_comf_exclude_diff
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_O2_COMF(-1, dim_obs_f)
   ! offset
   offset_day(iO2comf,iday+1)=offset_day(iO2comf,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing COMFORT-O2: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_O2_COMF(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f  = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%numrep   , thisobs_PP_f%numrep   , stats)
      
      ! InnoOmit from global observed forecast
      IF (assim_o_o2_comf) THEN
         allocate(oforc_f(thisobs%dim_obs_f))
         call obs_op_O2_COMF(dim_state_p, thisobs%dim_obs_f, forc_p, oforc_f)
         WHERE (ABS(oforc_f-thisobs%obs_f)>o2_comf_exclude_diff) &
            thisobs_PP_f%isInnoOmit = 1.0
      ENDIF
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iO2comf),offset_day(iO2comf,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f)
      IF (assim_o_o2_comf) deallocate(oforc_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iO2comf),offset_day(iO2comf,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth, &
              thisobs_PP%numrep)
   
END SUBROUTINE PP_O2_COMFORT


! *********************
! *** PP_pCO2_SOCAT ***
! *********************
SUBROUTINE PP_PCO2_SOCAT()

   use obs_pCO2_SOCAT_pdafomi, &
       only: init_dim_obs_pCO2_SOCAT, &
       obs_op_pCO2_SOCAT, pCO2_SOCAT_exclude_diff, &
       thisobs, thisobs_PP, thisobs_PP_f
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_pCO2_SOCAT(-1, dim_obs_f)
   ! offset
   offset_day(ipCO2s,iday+1)=offset_day(ipCO2s,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing SOCAT-pCO2: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_pCO2_SOCAT(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%numrep   , thisobs_PP_f%numrep   , stats)
      
      if (pCO2_SOCAT_exclude_diff > 0.0) then
            allocate(oforc_f(thisobs%dim_obs_f))
            call obs_op_pCO2_SOCAT(dim_state_p, thisobs%dim_obs_f, forc_p, oforc_f)
            WHERE (ABS(oforc_f-thisobs%obs_f)>pCO2_SOCAT_exclude_diff) &
               thisobs_PP_f%isInnoOmit = 1.0
      endif
      
      ! write output
      if (writepe) call write_ncfile(fidobs(ipCO2s),offset_day(ipCO2s,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f, oforc_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(ipCO2s),offset_day(ipCO2s,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat, &
              thisobs_PP%numrep)
   
END SUBROUTINE PP_PCO2_SOCAT


! **********************
! *** PP_DIN_COMFORT ***
! **********************
SUBROUTINE PP_DIN_COMFORT()

   use obs_n_comf_pdafomi, &
       only: init_dim_obs_n_comf, &
       obs_op_n_comf, &
       thisobs, thisobs_PP, thisobs_PP_f
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_n_comf(-1, dim_obs_f)
   ! offset
   offset_day(iNcomf,iday+1)=offset_day(iNcomf,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing COMFORT-DIN: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_n_comf(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iNcomf),offset_day(iNcomf,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iNcomf),offset_day(iNcomf,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth)
   
END SUBROUTINE PP_DIN_COMFORT


! **********************
! *** PP_DIN_ARGO    ***
! **********************
SUBROUTINE PP_DIN_ARGO()

   use obs_n_argo_pdafomi, &
       only: init_dim_obs_n_argo, &
       obs_op_n_argo, &
       thisobs, thisobs_PP, thisobs_PP_f
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_n_argo(-1, dim_obs_f)
   ! offset
   offset_day(iNargo,iday+1)=offset_day(iNargo,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing ARGO-DIN: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_n_argo(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iNargo),offset_day(iNargo,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iNargo),offset_day(iNargo,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth, &
              thisobs_PP%numrep)
   
END SUBROUTINE PP_DIN_ARGO


! **********************
! *** PP_O2_ARGO     ***
! **********************
SUBROUTINE PP_O2_ARGO()

   use obs_o2_argo_pdafomi, &
       only: init_dim_obs_o2_argo, &
       obs_op_o2_argo, &
       thisobs, thisobs_PP, thisobs_PP_f
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
  USE PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   ! init observation data
   call init_dim_obs_o2_argo(-1, dim_obs_f)
   ! offset
   offset_day(iOargo,iday+1)=offset_day(iOargo,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing ARGO-OXY: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_o2_argo(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iOargo),offset_day(iOargo,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iOargo),offset_day(iOargo,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth, &
              thisobs_PP%numrep)
   
END SUBROUTINE PP_O2_ARGO


! **********************
! *** PP_DIN_MERGED     ***
! **********************
SUBROUTINE PP_DIN_MERGED()

   use obs_n_merged_pdafomi, &
       only: init_dim_obs_n_merged, &
       obs_op_n_merged, &
       thisobs, thisobs_PP, thisobs_PP_f, &
       mean_n_p
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
   use PDAFomi, &
       only: PDAFomi_set_debug_flag
       
   implicit none
   REAL :: excl_relative_upper, excl_relative_lower ! relative upper and lower limits for exclusion
   
   if (n_merged_excl_relativePP > 0.0) then
      excl_relative_lower = 1.0 - n_merged_excl_relativePP
      excl_relative_upper = 1.0 / excl_relative_lower
   endif
   
   ! provide forecast data
   allocate(mean_n_p(dim_fields(id%DIN)))
   mean_n_p = forc_p(offset(id%DIN)+1 : offset(id%DIN)+dim_fields(id%DIN))
   ! init observation data
   call init_dim_obs_n_merged(-1, dim_obs_f)
   ! offset
   offset_day(iNmerged,iday+1)=offset_day(iNmerged,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing MERGED-DIN: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_n_merged(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      ! InnoOmit from global observed forecast
      allocate(oforc_f(thisobs%dim_obs_f))
      call obs_op_n_merged(dim_state_p, thisobs%dim_obs_f, forc_p, oforc_f)
      ! absolute crit:
      if (n_merged_excl_absolutePP > 0.0) then
         WHERE (ABS(oforc_f-thisobs%obs_f)>n_merged_excl_absolutePP) &
               thisobs_PP_f%isInnoOmit = 1.0
      endif
      ! relative crit:
      if (n_merged_excl_relativePP > 0.0) then
         WHERE (thisobs%obs_f>(excl_relative_upper*oforc_f)) &
               thisobs_PP_f%isInnoOmit = 1.0
         WHERE (thisobs%obs_f<(excl_relative_lower*oforc_f)) &
               thisobs_PP_f%isInnoOmit = 1.0
      endif
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iNmerged),offset_day(iNmerged,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f,oforc_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iNmerged),offset_day(iNmerged,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth, &
              thisobs_PP%numrep,mean_n_p)
   
END SUBROUTINE PP_DIN_MERGED


! **********************
! *** PP_O2_MERGED     ***
! **********************
SUBROUTINE PP_O2_MERGED()

   use obs_o2_merged_pdafomi, &
       only: init_dim_obs_o2_merged, &
       obs_op_o2_merged, &
       thisobs, thisobs_PP, thisobs_PP_f, &
       mean_o2_p
   use PDAFomi_obs_f, &
       only: PDAFomi_gather_obs_f_flex
   use PDAFomi_obs_l, &
       only: PDAFomi_deallocate_obs
   use PDAFomi, &
       ONLY: PDAFomi_set_debug_flag
       
   implicit none
   REAL :: excl_relative_upper, excl_relative_lower ! relative upper and lower limits for exclusion
   
   if (o2_merged_excl_relativePP > 0.0) then
      excl_relative_lower = 1.0 - o2_merged_excl_relativePP
      excl_relative_upper = 1.0 / excl_relative_lower
   endif
   
   ! provide forecast data
   allocate(mean_O2_p(dim_fields(id%O2)))
   mean_O2_p = forc_p(offset(id%O2)+1 : offset(id%O2)+dim_fields(id%O2))
   ! init observation data
   call init_dim_obs_o2_merged(-1, dim_obs_f)
   ! offset
   offset_day(iOmerged,iday+1)=offset_day(iOmerged,iday)+thisobs%dim_obs_f
   
   IF (writepe) &
        WRITE (*, '(a, 5x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Postprocessing MERGED-OXY: Number of observations is', thisobs%dim_obs_f
   
   haveobs: if (thisobs%dim_obs_f>0) then
   
      ! global observed state
      allocate(ostate_f   (2,thisobs%dim_obs_f))
      do isim=1,2
         call obs_op_o2_merged(dim_state_p, thisobs%dim_obs_f, state_p   (isim,:), ostate_f   (isim,:))
      enddo
      
      ! compute model-observation difference and improvement
      allocate(diff_f    (2,thisobs%dim_obs_f))
      allocate(impr_f      (thisobs%dim_obs_f))
      allocate(diffAF_f    (thisobs%dim_obs_f))
      
      do isim=1,2        
         diff_f   (isim,:) = thisobs%obs_f - ostate_f   (isim,:)
      enddo
      impr_f    = abs(diff_f   (FRRN,:)) - abs(diff_f   (ASML,:))
      diffAF_f    = ostate_f   (ASML,:) - ostate_f   (FRRN,:)
      
      ! gather further observation info
      allocate(thisobs_PP_f%depth     (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%isExclObs (thisobs%dim_obs_f), source = 0.0)
      allocate(thisobs_PP_f%nod1_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod2_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nod3_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%elem_g    (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lon       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%lat       (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%nz        (thisobs%dim_obs_f))
      allocate(thisobs_PP_f%numrep    (thisobs%dim_obs_f), source = fill_value)
      allocate(thisobs_PP_f%isInnoOmit(thisobs%dim_obs_f), source = 0.0)
      
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%depth    , thisobs_PP_f%depth    , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%isExclObs, thisobs_PP_f%isExclObs, stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod1_g   , thisobs_PP_f%nod1_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod2_g   , thisobs_PP_f%nod2_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nod3_g   , thisobs_PP_f%nod3_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%elem_g   , thisobs_PP_f%elem_g   , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lon      , thisobs_PP_f%lon      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%lat      , thisobs_PP_f%lat      , stats)
      CALL PDAFomi_gather_obs_f_flex(thisobs%dim_obs_p, thisobs_PP%nz       , thisobs_PP_f%nz       , stats)
      
      ! InnoOmit from global observed forecast
      allocate(oforc_f(thisobs%dim_obs_f))
      call obs_op_O2_merged(dim_state_p, thisobs%dim_obs_f, forc_p, oforc_f)
      ! absolute crit:
      if (o2_merged_excl_absolutePP > 0.0) then
         WHERE (ABS(oforc_f-thisobs%obs_f)>o2_merged_excl_absolutePP) &
               thisobs_PP_f%isInnoOmit = 1.0
      endif
      ! relative crit:
      if (o2_merged_excl_relativePP > 0.0) then
         WHERE (thisobs%obs_f>(excl_relative_upper*oforc_f)) &
               thisobs_PP_f%isInnoOmit = 1.0
         WHERE (thisobs%obs_f<(excl_relative_lower*oforc_f)) &
               thisobs_PP_f%isInnoOmit = 1.0
      endif
      
      ! write output
      if (writepe) call write_ncfile(fidobs(iOmerged),offset_day(iOmerged,iday),thisobs,thisobs_PP_f, &
                                     ostate_f,   diff_f,   impr_f,   diffAF_f)

      ! clean up after have_obs
      deallocate(ostate_f,diff_f,impr_f,diffAF_f,oforc_f)
      deallocate(thisobs_PP_f%isExclObs,thisobs_PP_f%nod1_g, &
                 thisobs_PP_f%nod2_g,thisobs_PP_f%nod3_g, &
                 thisobs_PP_f%elem_g,thisobs_PP_f%lon, &
                 thisobs_PP_f%lat,thisobs_PP_f%nz,thisobs_PP_f%depth, &
                 thisobs_PP_f%numrep,thisobs_PP_f%isInnoOmit)
   else haveobs
      ! no observations: update output file
      if (writepe) call nowrite_ncfile(fidobs(iOmerged),offset_day(iOmerged,iday))
   endif haveobs

   ! clean up after init_dim_obs
   call PDAFomi_deallocate_obs(thisobs)
   deallocate(thisobs_PP%isExclObs,thisobs_PP%nod1_g, &
              thisobs_PP%nod2_g,thisobs_PP%nod3_g, &
              thisobs_PP%elem_g,thisobs_PP%lon, &
              thisobs_PP%lat,thisobs_PP%nz,thisobs_PP%depth, &
              thisobs_PP%numrep,mean_o2_p)
   
END SUBROUTINE PP_O2_MERGED


! *********************
! *** write_ncfile ***
! *********************
! to write output
SUBROUTINE write_ncfile(fid,self_offset_day,self_thisobs,self_thisobsPP, &
                        self_ostate,self_diff,self_impr,self_diffAF)

   USE recom_config,   ONLY: SecondsPerDay
   USE PDAFomi,        ONLY: obs_f          ! Type variable for thisobs
   USE mod_assim_pdaf, ONLY: obs_PP         ! Type variable for thisobs_PP

   integer, intent(in) :: fid               ! file ID
   integer, intent(in) :: self_offset_day   ! write position for observations
   
   type(obs_f)  , intent(in)  :: self_thisobs    ! observation info
   type(obs_PP) , intent(in)  :: self_thisobsPP  ! further observation info
   
   real, intent(in) :: self_ostate(:,:)     ! observed state
   real, intent(in) :: self_diff(:,:)       ! model observation differences
   real, intent(in) :: self_impr(:)         ! improvement
   real, intent(in) :: self_diffAF(:)       ! difference ASML minus FREE
   
   integer, allocatable :: obsindex(:)      ! observation index
   integer, allocatable :: DayInSec(:)      ! day, in seconds since start of year, to specify date
   integer              :: varid
   
   ! set observation index
   allocate(obsindex(self_thisobs%dim_obs_f))
   do i=1,self_thisobs%dim_obs_f
      obsindex(i) = self_offset_day +i -1
   enddo ! i
   
   ! date and index
   allocate(DayInSec(self_thisobs%dim_obs_f))
   DayInSec = (daynew-1)*SecondsPerDay
   
   call check( nf90_inq_varid(fid, "day",  varid))
   call check( nf90_put_var  (fid, varid, (daynew-1)*SecondsPerDay, &
                              start=(/ daynew /)))
   
   call check( nf90_inq_varid(fid, "offset_day", varid))
   call check( nf90_put_var  (fid, varid, self_offset_day, &
                              start=(/ daynew /)))
   
   call check( nf90_inq_varid(fid, "numobs_day", varid))
   call check( nf90_put_var  (fid, varid, self_thisobs%dim_obs_f, start=(/ daynew /)))
   
   call check( nf90_inq_varid(fid, "index", varid))
   call check( nf90_put_var  (fid, varid, obsindex, &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
   
   call check( nf90_inq_varid(fid, "time", varid))
   call check( nf90_put_var  (fid, varid, DayInSec, &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
   
   ! data
   ! observation
   call check( nf90_inq_varid(fid, "OBS", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_thisobs%obs_f,4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
   
   ! observed state same day
   ! nsteps_per_day-1 forecast steps leading up to analysis, and analysis step
   call check( nf90_inq_varid(fid, "ASML", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_ostate(ASML,:),4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "FREE", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_ostate(FRRN,:),4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "diffOASML", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_diff(ASML,:),4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "diffOFREE", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_diff(FRRN,:),4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "diffASMLFREE", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_diffAF,4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "improvement", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_impr,4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
   
   ! observation metadata
   call check( nf90_inq_varid(fid, "depth", varid))
   call check( nf90_put_var  (fid, varid, self_thisobsPP%depth, &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "NOD1_G", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%nod1_g), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "NOD2_G", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%nod2_g), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "NOD3_G", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%nod3_g), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "ELEM_G", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%elem_g), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))

   call check( nf90_inq_varid(fid, "NZ1", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%nz), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))

   call check( nf90_inq_varid(fid, "numrep", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%numrep), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))

   call check( nf90_inq_varid(fid, "lon", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_thisobsPP%lon,4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "lat", varid))
   call check( nf90_put_var  (fid, varid, REAL(self_thisobsPP%lat,4), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "isExclObs", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%isExclObs,to_nf90_byte), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))
                              
   call check( nf90_inq_varid(fid, "isInnoOmit", varid))
   call check( nf90_put_var  (fid, varid, INT(self_thisobsPP%isInnoOmit,to_nf90_byte), &
                              start=(/ self_offset_day /), count=(/ self_thisobs%dim_obs_f /)))

   deallocate(obsindex,DayInSec)
   
END SUBROUTINE write_ncfile


! *********************
! *** nowrite_ncfile ***
! *********************
! to write output
SUBROUTINE nowrite_ncfile(fid,self_offset_day)

   USE recom_config, ONLY: SecondsPerDay

   integer, intent(in) :: fid               ! file ID
   integer, intent(in) :: self_offset_day   ! write position for observations
   
   integer :: varid
   integer :: self_numobs_day = 0
   
   WRITE (*, '(a, 5x, a, 1x, i, 1x, a, 1x, i)') 'FESOM-PDAF', &
        '--- Writing to netCDF at', self_offset_day, 'number of observations', 0
   
   call check( nf90_inq_varid(fid, "day",  varid))
   call check( nf90_put_var  (fid, varid, (daynew-1)*SecondsPerDay, start=(/ daynew /)))
   
   call check( nf90_inq_varid(fid, "offset_day", varid))
   call check( nf90_put_var  (fid, varid, self_offset_day, start=(/ daynew /)))
   
   call check( nf90_inq_varid(fid, "numobs_day", varid))
   call check( nf90_put_var  (fid, varid, self_numobs_day, start=(/ daynew /)))

END SUBROUTINE


! *********************
! *** create_ncfile ***
! *********************
! to write output
SUBROUTINE create_ncfile(obstype,fid)
   
   use netcdf
   
   character(len=20), intent(in) :: obstype
   integer, intent(out) :: fid
   character(len=200) :: filename
   character(40) :: att_text
   integer :: dimID_iter, dimID_day, varid, varid_time, varid_day
   
   ! create file
   filename = trim(DAoutput_path)//'obscompare_'//trim(obstype)//'.'//cyearnew//'.nc'
   call check(nf90_create(filename, NF90_CLOBBER, fid))
   
   ! define dimensions
   call check(NF90_DEF_DIM(fid, 'day'  , ndaysim       , dimID_day ))
   call check(NF90_DEF_DIM(fid, 'index', NF90_UNLIMITED, dimID_iter))
   
   ! define coordinates
   call check(nf90_def_var(fid, 'day'  , NF90_INT, dimID_day , varid_day))
   call check(nf90_def_var(fid, 'index', NF90_INT, dimID_iter, varid))
   
   call check(nf90_def_var(fid, 'offset_day', NF90_INT, dimID_day , varid))
   call check(nf90_def_var(fid, 'numobs_day', NF90_INT, dimID_day , varid))
   
   ! define variables
   call check(nf90_def_var(fid, "time",         NF90_INT,     dimID_iter, varid_time))
   
   call check(nf90_def_var(fid, "OBS",          NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "isExclObs",    NF90_BYTE,    dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(-1,to_nf90_byte)))
   call check(nf90_put_att(fid, varid, 'Description', '1=Excluded; 0=Assimilated; -1=FillValue'))
   call check(nf90_def_var(fid, "isInnoOmit",   NF90_BYTE,    dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(-1,to_nf90_byte)))
   call check(nf90_put_att(fid, varid, 'Description', '1=Excluded; 0=Assimilated; -1=FillValue'))
   call check(nf90_def_var(fid, "ASML",         NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "FREE",         NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "diffOASML",    NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "diffOFREE",    NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "diffASMLFREE", NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "improvement",  NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "NOD1_G",       NF90_INT,     dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(fill_value)))
   call check(nf90_def_var(fid, "NOD2_G",       NF90_INT,     dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(fill_value)))
   call check(nf90_def_var(fid, "NOD3_G",       NF90_INT,     dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(fill_value)))
   call check(nf90_def_var(fid, "ELEM_G",       NF90_INT,     dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(fill_value)))
   call check(nf90_def_var(fid, "NZ1",          NF90_INT,     dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(fill_value)))
   call check(nf90_def_var(fid, "depth",        NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "lon",          NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "lat",          NF90_FLOAT,   dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', REAL(fill_value,4)))
   call check(nf90_def_var(fid, "numrep",       NF90_INT,     dimID_iter, varid))
   call check(nf90_put_att(fid, varid, '_FillValue', INT(fill_value)))
   
   ! standardized attributes time variables
   write(att_text, '(a14,I4.4,a1,I2.2,a1,I2.2,a6)') 'seconds since ', yearnew, '-', 1, '-', 1, ' 0:0:0'
   
   call check( nf90_put_att(fid, varid_time, 'long_name', 'time'))
   call check( nf90_put_att(fid, varid_time, 'standard_name', 'time'))
   call check( nf90_put_att(fid, varid_time, 'units', trim(att_text)))
   call check( nf90_put_att(fid, varid_time, 'axis', 'T'))
   call check( nf90_put_att(fid, varid_time, 'stored_direction', 'increasing'))
   
   call check( nf90_put_att(fid, varid_day, 'long_name', 'time'))
   call check( nf90_put_att(fid, varid_day, 'standard_name', 'time'))
   call check( nf90_put_att(fid, varid_day, 'units', trim(att_text)))
   call check( nf90_put_att(fid, varid_day, 'axis', 'T'))
   call check( nf90_put_att(fid, varid_day, 'stored_direction', 'increasing'))
   
   ! finish
   call check(nf90_enddef(fid))

END SUBROUTINE

END MODULE
