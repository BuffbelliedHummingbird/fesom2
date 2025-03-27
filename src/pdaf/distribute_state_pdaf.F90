! PDAF in AWI-CM2 / Fesom 2.0

SUBROUTINE distribute_state_pdaf(dim_p, state_p)

! DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! During the forecast phase of the filter this
! subroutine is called from PDAF\_get\_state
! supplying a model state which has to be evolved. 
! The routine has to initialize the fields of the 
! model (typically available through a module) from 
! the state vector of PDAF. With parallelization, 
! MPI communication might be required to 
! initialize all subdomains on the model PEs.
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! !REVISION HISTORY:
! 2004-10 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2020-03 - Frauke B     - Added velocities, removed sea-ice (FESOM2.0)

  !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_submodel, mype_world, task_id, mype_model, npes_model, &
       COMM_model, MPIerr, mype_filter, writepe
  USE mod_assim_pdaf, &
       ONLY: offset, loc_radius,this_is_pdaf_restart, mesh_fesom, nlmax, &
       dim_fields, id, istep_asml, step_null, start_from_ENS_spinup, &
       topography_p
  USE mod_nc_out_variables, &
       ONLY: sfields, nfields_tr3D, ids_tr3D
  USE g_PARSUP, &
       ONLY: myDim_nod2D, myDim_elem2D, &
             eDim_nod2D, eDim_elem2D 
  USE o_arrays, &
       ONLY: eta_n, uv, wvel, tr_arr, Tsurf, Ssurf, Unode
  USE i_arrays, &
       ONLY: a_ice
  USE g_clock, &
       ONLY: daynew,timenew
  USE g_comm_auto                        ! contains: interface exchange_nod()

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
! ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_dist_state)

! Local variables
  INTEGER :: i, k, b, s, istate, ifesom   ! Counter
  INTEGER :: node                         ! Node index
  
  REAL, ALLOCATABLE :: U_node_upd(:,:,:) ! Velocity update on nodes
  REAL, ALLOCATABLE :: U_elem_upd(:,:,:) ! Velocity update on elements
  
  LOGICAL, save :: first_call = .true.
  

! Debugging:
  LOGICAL            :: debugmode
  LOGICAL            :: write_debug
  INTEGER            :: fileID_debug
  CHARACTER(len=5)   :: mype_string
  CHARACTER(len=3)   :: day_string
  CHARACTER(len=5)   :: tim_string
  
  ! Set debug output
  debugmode    = .false.
  IF (.not. debugmode) THEN
     write_debug = .false.
  ELSE
     IF (mype_world>0) THEN
        write_debug = .false.
     ELSE
        write_debug = .true.
     ENDIF
  ENDIF
    
  IF (write_debug) THEN
         ! print state vector
         WRITE(day_string, '(i3.3)') daynew
         WRITE(tim_string, '(i5.5)') int(timenew)
         fileID_debug=20
         open(unit=fileID_debug, file='distribute_state_pdaf_'//day_string//'_'//tim_string//'.txt', status='unknown')
  ENDIF

! **********************
! *** Initialization ***
! **********************
! no need to distibute the ensemble in case of restart, just skip this routine:

  IF (this_is_pdaf_restart .AND. (istep_asml==step_null)) THEN
    if (mype_world==0) WRITE(*,*) 'FESOM-PDAF This is a restart: Skipping distribute_state_pdaf at initial step'
    
  ELSEIF (start_from_ENS_spinup .AND. (istep_asml==step_null)) THEN
    if (mype_world==0) WRITE(*,*) 'FESOM-PDAF Model starts from perturbed ensemble: Skipping distribute_state_pdaf at initial step'
    
  ELSE
    if (mype_submodel==0) write (*,*) 'FESOM-PDAF distribute_state_pdaf, task: ', task_id
    
  ! ensure to distribute fields with valid topography at first call
  IF (first_call) THEN
      IF (writepe) WRITE (*,'(a, 8x,a)') 'FESOM-PDAF', 'Distribute_state: set topography'
      state_p = state_p * topography_p
      first_call = .FALSE.
  END IF

! *******************************************
! *** Initialize model fields from state  ***
! *** Each model PE knows its sub-state   ***
!********************************************

  ! *** Dimensions of FESOM arrays:
  ! * eta_n          (myDim_nod2D + eDim_nod2D)            ! Dynamic topography
  ! * UV(1,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity u
  ! * UV(2,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity v
  ! * wvel           (nl,   myDim_nod2D + eDim_nod2D)      ! Velocity w
  ! * tr_arr(:,:,1)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Temperature
  ! * tr_arr(:,:,2)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Salinity
  ! * unode(1,:,:)   (1, nl-1, myDim_nod2D + eDim_nod2D)   ! Velocity u interpolated on nodes
  ! * unode(2,:,:)   (1, nl-1, myDim_nod2D + eDim_nod2D)   ! Velocity v interpolated on nodes
  ! * a_ice          (myDim_nod2D + eDim_nod2D)            ! Sea-ice concentration
  ! ***
  
  ! SSH (1)
  DO i = 1, myDim_nod2D
     s = i + offset(id% SSH)
     if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%SSH)%variable, s, eta_n(i), state_p(s)
     eta_n(i) = state_p(s)
  END DO
  
  ! u (2) and v (3) velocities
  ! 1. calculate update on nodes, i.e. analysis state (state_p) minus not-yet-updated model state (Unode)
  allocate(U_node_upd(2, mesh_fesom%nl-1, myDim_nod2D+eDim_nod2D))
  U_node_upd = 0.0
  
  DO i = 1, myDim_nod2D
   DO k = 1, nlmax
      ! u
      s = (i-1) * (nlmax) + k + offset(id% u)
      U_node_upd(1, k, i) = state_p(s) - Unode(1, k, i)
      if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%u)%variable, s, Unode(1, k, i), state_p(s)
      ! v
      s = (i-1) * (nlmax) + k + offset(id% v)
      U_node_upd(2, k, i) = state_p(s) - Unode(2, k, i)
      if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%v)%variable, s, Unode(2, k, i), state_p(s)
   END DO
  END DO
  
  ! 2. interpolate update from nodes to elements
  allocate(U_elem_upd(2, mesh_fesom%nl-1, myDim_elem2D+eDim_elem2D))
  U_elem_upd = 0.0

  call compute_vel_elems(U_node_upd,U_elem_upd)
  
  ! 3. add update to model velocity on elements (UV)
  UV = UV + U_elem_upd
  
  ! 4. adjust diagnostic model velocity on nodes (Unode)
  call compute_vel_nodes(mesh_fesom)


  ! Element-wise version not interpolated onto nodes. Removed in favor of the above:
  ! DO i = 1, myDim_elem2D
  !  DO k = 1, nlmax
  !      UV(1,k,i) = state_p((i-1)*(nlmax) + k + offset(2)) ! u
  !      UV(2,k,i) = state_p((i-1)*(nlmax) + k + offset(3)) ! v
  !  END DO
  ! END DO


  ! w (4) velocity: not updated and thus no need to distribute.
  ! DO i = 1, myDim_nod2D
  !  DO k = 1, nlmax
  !      wvel(k,i) = state_p((i-1)*nlmax + k + offset(id% w)) ! w
  !  END DO
  ! END DO
  
  ! Temp (5) and salt (6)
  ! are included in tracer field loop.
  
  ! Sea-ice concentration is needed in PDAF to not assimilate SST at
  ! sea-ice locations.
  ! But sea-ice itself is not assimilated, thus sea-ice update is
  ! not distributed to the model.
  
! *********************************
! *** Initialize external nodes ***
! *********************************

  call exchange_nod(eta_n)            ! SSH
  call exchange_elem(UV(:,:,:))       ! u and v (element-wise)

! *********************************
! *** Model 3D tracers          ***
! *********************************
  ! tracer field loop
  DO b = 1, nfields_tr3D
  
     istate = ids_tr3D(b)                 ! index of field in state vector
     ifesom = sfields(istate)%trnumfesom  ! index of field in model tracer array
     
     DO i = 1, myDim_nod2D
        DO k = 1, nlmax
        
           ! indeces to flatten 3D arrays
           s = (i-1) * (nlmax) + k
           ! put state values into model tracer array
           tr_arr(k, i,  ifesom) = state_p(s + offset(istate))
           
        ENDDO
     ENDDO
     ! initialize external nodes
     call exchange_nod(tr_arr(:,:,ifesom))
  ENDDO

  ! clean up:
  if (write_debug) close(fileID_debug)
  deallocate(U_node_upd,U_elem_upd)
 
  ENDIF
END SUBROUTINE distribute_state_pdaf
