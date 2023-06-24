!$Id: init_dim_l_pdaf.F90 2468 2021-02-26 08:02:43Z lnerger $
!>  Routine to set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!!
!! The routine is called during analysis step
!! in the loop over all local analysis domains.
!! It has to set the dimension of local model 
!! state on the current analysis domain.
!!
!! The routine is called by each filter process.
!!
!! __Revision history:__
!! 2005-09 - Lars Nerger - Initial code
!! 2022-03 - Frauke B    - Adapted for FESOM 2.1
!!

SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)

  USE mod_assim_pdaf, &                    ! Variables for assimilation
       ONLY: id_lstate_in_pstate, &        ! Indices of local state vector in PE-local global state vector
             offset, &                     ! PE-local offsets of fields in state vector
             id, &                         ! Field IDs in state vector
             coords_l, &                   ! Coordinates of local analysis domain
             mesh_fesom, &                 
             nfields, bgcmin, bgcmax
                                          ! mesh_fesom % coord_nod2D, & ! vertex coordinates in radian measure
                                          ! mesh_fesom % nlevels, &     ! number of levels at (below) elem     considering bottom topography
                                          ! mesh_fesom % nlevels_nod2D  ! number of levels at (below) vertices considering bottom topography
                                          ! mesh_fesom % nl             ! number of levels not considering bottom topography
    USE PDAFomi, &
        ONLY: PDAFomi_set_debug_flag
    USE mod_parallel_pdaf, &
        ONLY: mype_filter
    USE PDAF_mod_filter, &
        ONLY: state
    USE mod_nc_out_variables, &
        ONLY: sfields
       
  USE g_rotate_grid, &
       ONLY: r2g                           ! Transform from the mesh (rotated) coordinates 
                                           ! to geographical coordinates  
       ! glon, glat        :: [radian] geographical coordinates
       ! rlon, rlat        :: [radian] mesh (rotated) coordinates
       
  USE mod_assim_pdaf, &
       ONLY: debug_id_nod2
  USE g_parsup, &
       ONLY: myList_nod2D

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     ! Current time step
  INTEGER, INTENT(in)  :: domain_p ! Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    ! Local state dimension

! *** Local variables ***
  INTEGER :: nlay                              ! Number of layers for current domain
  INTEGER :: dim_fields_l(nfields)             ! Field dimensions for current domain
  INTEGER :: offset_l(nfields)                 ! Field offsets for current domain
  INTEGER :: i, b                              ! Counters
  
  ! integer :: myDebug_id(1)
  ! myDebug_id = FINDLOC(myList_nod2D, value=debug_id_nod2)
  
!~   if (mype_filter==0) write(*,*) 'init_dim_l_pdaf: domain_p', domain_p

! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  nlay = mesh_fesom%nlevels_nod2D(domain_p)-1
  
  dim_fields_l (id%SSH)    = 1
  dim_fields_l (id%u)      = nlay
  dim_fields_l (id%v)      = nlay
  dim_fields_l (id%w)      = 0 ! nlay
  dim_fields_l (id%temp)   = nlay
  dim_fields_l (id%salt)   = nlay
  dim_fields_l (id%a_ice)  = 0
  dim_fields_l (id%MLD1)   = 0

  DO b=bgcmin, bgcmax
    ! not updated:
    IF ( .not. (sfields(b)% updated)) THEN
      dim_fields_l(b) = 0
    ELSE
      ! surface fields:
      IF (sfields(b)% ndims == 1)   dim_fields_l(b)=1
      ! 3D fields:
      IF (sfields(b)% ndims == 2)   dim_fields_l(b)=nlay
    ENDIF
  ENDDO

  offset_l(1) = 0
  DO i = 2,nfields
     offset_l(i) = offset_l(i-1) + dim_fields_l(i-1)
  END DO
  
  dim_l = sum(dim_fields_l)

! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Get location of current water column (basis point)
  CALL r2g(coords_l(1), coords_l(2), mesh_fesom%coord_nod2D(1, domain_p), mesh_fesom%coord_nod2D(2, domain_p))
  
!~   IF (domain_p==myDebug_id(1)) THEN
!~   WRITE(*,*) 'OMI-debug (F): init_dim_l_pdaf: coords_l(1): ', coords_l(1), ' coords_l(2): ', coords_l(2)
!~   ENDIF

!~   IF ((mype_filter==55) .AND. (domain_p==669)) THEN
!~         open (2, file = 'dim_fields_l.dat')
!~         write(2,*) dim_fields_l
!~         close(2)
!~         open (3, file = 'offset_l.dat')
!~         write(3,*) offset_l
!~         close(3)
!~   END IF


! ****************************************************
! *** Initialize array of indices for local domain ***
! ****************************************************

  ! Allocate arrays
  IF (ALLOCATED(id_lstate_in_pstate)) DEALLOCATE(id_lstate_in_pstate)
  
  ALLOCATE(id_lstate_in_pstate(dim_l))

  ! *** indices for full state vector ***

  ! SSH
  id_lstate_in_pstate (offset_l(id%ssh)+1) &
        = offset(id%ssh) + domain_p
  ! U
  id_lstate_in_pstate (offset_l(id%u)+1 : offset_l(id%u)+dim_fields_l(id%u)) &
        = offset(id%u) &
        + (domain_p-1)*(mesh_fesom%nl-1) &
        + (/(i, i=1,dim_fields_l(id%u))/)
  ! V
  id_lstate_in_pstate (offset_l(id%v)+1 : offset_l(id%v)+dim_fields_l(id%v)) &
        = offset(id%v) &
        + (domain_p-1)*(mesh_fesom%nl-1) &
        + (/(i, i=1,dim_fields_l(id%v))/)
        
  ! W
  ! id_lstate_in_pstate (offset_l(id%w)+1 : offset_l(id%w+1)) &
  !      = offset(id%w) &
  !      + (domain_p-1)*(mesh_fesom%nl) &
  !      + (/(i, i=1,dim_fields_l(id%w))/)
  
  ! Temp
  id_lstate_in_pstate (offset_l(id%temp)+1 : offset_l(id%temp)+dim_fields_l(id%temp))&
         = offset(id%temp) &
         + (domain_p-1)*(mesh_fesom%nl-1) &
         + (/(i, i=1,dim_fields_l(id%temp))/)
  ! Salt
  id_lstate_in_pstate (offset_l(id%salt)+1 : offset_l(id%salt)+dim_fields_l(id%salt)) &
        = offset(id%salt) &
        + (domain_p-1)*(mesh_fesom%nl-1) &
        + (/(i, i=1,dim_fields_l(id%salt))/)

END SUBROUTINE init_dim_l_pdaf
