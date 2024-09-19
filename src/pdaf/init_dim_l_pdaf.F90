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

SUBROUTINE init_dim_l_pdaf(step, nsweeped_domain_p, dim_l)

  USE mod_assim_pdaf, &                    ! Variables for assimilation
       ONLY: id_lstate_in_pstate, &        ! Indices of local state vector in PE-local global state vector
             offset, &                     ! PE-local offsets of fields in state vector
             id, &                         ! Field IDs in state vector
             coords_l, &                   ! Coordinates of local analysis domain
             mesh_fesom, nlmax, &                 
             nfields, bgcmin, bgcmax, phymin, phymax, &
                                          ! mesh_fesom % coord_nod2D, & ! vertex coordinates in radian measure
                                          ! mesh_fesom % nlevels, &     ! number of levels at (below) elem     considering bottom topography
                                          ! mesh_fesom % nlevels_nod2D  ! number of levels at (below) vertices considering bottom topography
                                          ! mesh_fesom % nl             ! number of levels not considering bottom topography
             dim_fields_l, offset_l, &    ! Domain local (water column at one node) field dimensions
             isweep
    USE PDAFomi, &
        ONLY: PDAFomi_set_debug_flag
    USE mod_parallel_pdaf, &
        ONLY: mype_filter, abort_parallel
    USE PDAF_mod_filter, &
        ONLY: state
    USE mod_nc_out_variables, &
        ONLY: sfields
    USE g_parsup, &
      ONLY: myDim_nod2D, &      ! Process-local number of vertices
            myDim_elem2D        ! Process-local number of elements 
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
  INTEGER, INTENT(in)  :: step              ! Current time step
  INTEGER, INTENT(in)  :: nsweeped_domain_p ! Current local analysis domain, containing repititive sweeps
  INTEGER, INTENT(out) :: dim_l             ! Local state dimension

! *** Local variables ***
  INTEGER :: nlay                                ! Number of layers for current domain
  INTEGER :: i, b, p                             ! Counters
  INTEGER :: domain_p                            ! Local analysis domain accounting for multiple sweeps
  
  ! integer :: myDebug_id(1)
  ! myDebug_id = FINDLOC(myList_nod2D, value=debug_id_nod2)
  
! ********************************************************
! ***  Account for multi sweeps in local analysis loop ***
! ********************************************************

  if (nsweeped_domain_p <= myDim_nod2D) then
  
     domain_p = nsweeped_domain_p
     ! Set index of sweep
     isweep = 1
     
  else
  
     domain_p = nsweeped_domain_p - myDim_nod2D
     ! Set index of sweep
     isweep = 2
  end if

! ****************************************
! *** Initialize local state dimension ***
! ****************************************
  
  nlay = mesh_fesom%nlevels_nod2D(domain_p)-1
  
  IF (nlay > nlmax) THEN
  WRITE(*,*) 'FESOM-PDAF ', 'init_dim_l_pdaf ', 'domain_p ', domain_p, ' nlay exceeds layer bounds!'
  CALL abort_parallel()
  ENDIF
  
!~   dim_fields_l (id%SSH)    = 1
!~   dim_fields_l (id%u)      = nlay
!~   dim_fields_l (id%v)      = nlay
!~   dim_fields_l (id%w)      = 0 ! nlay
!~   dim_fields_l (id%temp)   = nlay
!~   dim_fields_l (id%salt)   = nlay
!~   dim_fields_l (id%a_ice)  = 0
!~   dim_fields_l (id%MLD1)   = 0
  
  ! Physics:
  DO p=phymin, phymax
    ! not updated:
    IF ( .not. (sfields(p)% updated)) THEN
      dim_fields_l(p) = 0
    ELSE
      ! surface fields:
      IF (sfields(p)% ndims == 1)   dim_fields_l(p)=1
      ! 3D fields:
      IF (sfields(p)% ndims == 2)   dim_fields_l(p)=nlay
    ENDIF
  ENDDO

  ! BGC:
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
  if (sfields(id%ssh)%updated) then
  id_lstate_in_pstate (offset_l(id%ssh)+1) &
        = offset(id%ssh) + domain_p
  endif
  
  ! U
  if (sfields(id%u)%updated) then
  id_lstate_in_pstate (offset_l(id%u)+1 : offset_l(id%u)+dim_fields_l(id%u)) &
        = offset(id%u) &
        + (domain_p-1)*(nlmax) &
        + (/(i, i=1,dim_fields_l(id%u))/)
  endif
  
  ! V
  if (sfields(id%v)%updated) then
  id_lstate_in_pstate (offset_l(id%v)+1 : offset_l(id%v)+dim_fields_l(id%v)) &
        = offset(id%v) &
        + (domain_p-1)*(nlmax) &
        + (/(i, i=1,dim_fields_l(id%v))/)
  endif
        
  ! W
  ! id_lstate_in_pstate (offset_l(id%w)+1 : offset_l(id%w+1)) &
  !      = offset(id%w) &
  !      + (domain_p-1)*(nlmax) &
  !      + (/(i, i=1,dim_fields_l(id%w))/)
  
  ! Temp
  if (sfields(id%temp)%updated) then
  id_lstate_in_pstate (offset_l(id%temp)+1 : offset_l(id%temp)+dim_fields_l(id%temp))&
         = offset(id%temp) &
         + (domain_p-1)*(nlmax) &
         + (/(i, i=1,dim_fields_l(id%temp))/)
  endif
  
  ! Salt
  if (sfields(id%salt)%updated) then
  id_lstate_in_pstate (offset_l(id%salt)+1 : offset_l(id%salt)+dim_fields_l(id%salt)) &
        = offset(id%salt) &
        + (domain_p-1)*(nlmax) &
        + (/(i, i=1,dim_fields_l(id%salt))/)
  endif
        
  ! BGC:
  DO b=bgcmin, bgcmax
  
    ! only updated fields:
    IF ((sfields(b)% updated)) THEN
      
      ! surface fields:
      IF (sfields(b)% ndims == 1)   THEN
            id_lstate_in_pstate (offset_l(b)+1) &
                   = offset(b) + domain_p
      ENDIF

      ! 3D fields:
      IF (sfields(b)% ndims == 2)   THEN
            id_lstate_in_pstate (offset_l(b)+1 : offset_l(b)+dim_fields_l(b))&
                   = offset(b) &
                   + (domain_p-1)*(nlmax) &
                   + (/(i, i=1,dim_fields_l(b))/)
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE init_dim_l_pdaf
