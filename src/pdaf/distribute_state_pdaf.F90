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
       COMM_model, MPIerr, mype_filter
  USE mod_assim_pdaf, &
       ONLY: offset, loc_radius,this_is_pdaf_restart, mesh_fesom, &
       dim_fields, id, istep_asml, step_null
  USE g_PARSUP, &
       ONLY: myDim_nod2D, myDim_elem2D, &
             eDim_nod2D, eDim_elem2D 
  USE o_arrays, &
       ONLY: eta_n, uv, wvel, tr_arr, Tsurf, Ssurf, Unode
  USE i_arrays, &
       ONLY: a_ice
  USE g_clock, &
       ONLY: daynew
  USE g_comm_auto                        ! contains: interface exchange_nod()

  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  
! ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! PE-local state vector

! CALLING SEQUENCE:
! Called by: PDAF_get_state   (as U_dist_state)

! Local variables
  INTEGER :: i, k         ! Counter
  INTEGER :: node         ! Node index
  
  REAL, ALLOCATABLE :: U_node_upd(:,:,:) ! Velocity update on nodes
  REAL, ALLOCATABLE :: U_elem_upd(:,:,:) ! Velocity update on elements

! **********************
! *** Initialization ***
! **********************
! no need to distibute the ensemble in case of restart, just skip this routine:

  IF (this_is_pdaf_restart .AND. (istep_asml==0)) THEN
    if (mype_submodel==0) WRITE(*,*) 'FESOM-PDAF This is a restart: Skipping distribute_state_pdaf'
  ELSE
    if (mype_submodel==0) write (*,*) 'FESOM-PDAF distribute_state_pdaf, task: ', task_id

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
     eta_n(i) = state_p(i + offset(id% SSH))
  END DO
  
  ! u (2) and v (3) velocities
  ! 1. calculate update on nodes, i.e. analysis state (state_p) minus not-yet-updated model state (Unode)
  allocate(U_node_upd(2, mesh_fesom%nl-1, myDim_nod2D+eDim_nod2D))
  U_node_upd = 0.0
  
  DO i = 1, myDim_nod2D
   DO k = 1, mesh_fesom%nl-1
      U_node_upd(1, k, i) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% u)) - Unode(1, k, i) ! u
      U_node_upd(2, k, i) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% v)) - Unode(2, k, i) ! v
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
  !  DO k = 1, mesh_fesom%nl-1
  !      UV(1,k,i) = state_p((i-1)*(mesh_fesom%nl-1) + k + offset(2)) ! u
  !      UV(2,k,i) = state_p((i-1)*(mesh_fesom%nl-1) + k + offset(3)) ! v
  !  END DO
  ! END DO


  ! w (4) velocity: not updated and thus no need to distribute.
  ! DO i = 1, myDim_nod2D
  !  DO k = 1, mesh_fesom%nl
  !      wvel(k,i) = state_p((i-1)*mesh_fesom%nl + k + offset(id% w)) ! w
  !  END DO
  ! END DO
  
  ! Temp (5) and salt (6)
  
  DO i = 1, myDim_nod2D
   DO k = 1, mesh_fesom%nl-1
      tr_arr(k, i, 1) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% temp)) ! T
      tr_arr(k, i, 2) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% salt)) ! S
   END DO
  END DO
  
  ! Sea-ice concentration is needed in PDAF to not assimilate SST at
  ! sea-ice locations.
  ! But sea-ice itself is not assimilated, thus sea-ice update is
  ! not distributed to the model.
  
! *********************************
! *** Initialize external nodes ***
! *********************************

  call exchange_nod(eta_n)            ! SSH
  call exchange_elem(UV(:,:,:))       ! u and v (element-wise)
  ! call exchange_nod(wvel)           ! No need as vertical velocities are not distributed.
  call exchange_nod(tr_arr(:,:,1))    ! Temp
  call exchange_nod(tr_arr(:,:,2))    ! Salt
  
! *********************************
! *** Biogeochemistry           ***
! *********************************

!~ ! 3D fields:
!~   DO i = 1, myDim_nod2D
!~    DO k = 1, mesh_fesom%nl-1
!~       tr_arr(k, i,  8) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyChl ))
!~       tr_arr(k, i, 17) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaChl ))
!~       tr_arr(k, i,  4) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DIC    ))
!~       tr_arr(k, i, 14) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DOC    ))
!~       tr_arr(k, i,  5) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Alk    ))
!~       tr_arr(k, i,  3) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DIN    ))
!~       tr_arr(k, i, 13) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DON    ))
!~       tr_arr(k, i, 24) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% O2     ))
!~       tr_arr(k, i,  6) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyN   ))
!~       tr_arr(k, i,  7) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaC   ))
!~       tr_arr(k, i, 15) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PAR    ))
!~       tr_arr(k, i, 18) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DetC   ))
!~ !      PAR3D (k, i)     = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PAR    )) ! diagnostic field (not distributed to the model. See int_recom/recom_sms.F90)
!~       tr_arr(k, i, 12) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% HetC   ))
!~       tr_arr(k, i, 10) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DetC   ))
!~ !      diags3D(k, i, 1) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% NPPn   )) ! diagnostic field (not distributed to the model. See int_recom/recom_sms.F90)
!~ !      diags3D(k, i, 2) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% NPPd   )) ! diagnostic field (not distributed to the model. See int_recom/recom_sms.F90)
!~       tr_arr(k, i, 22) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyCalc))
!~            ! small, fast-growing zooplankton
!~            tr_arr(k, i, 12) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo1C)) ! intracellular conc of carbon in zooplankton 1
!~            tr_arr(k, i, 11)state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo1N))    ! intracellular conc of nitrogen in zooplankton 1
!~            ! macrozooplankton/antarctic krill:
!~            tr_arr(k, i, 26) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo2C)) ! intracellular conc of carbon in zooplankton 2
!~            tr_arr(k, i, 27) = state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo2N)) ! intracellular conc of nitrogen in zooplankton 2
!~    END DO
!~   END DO
  
!~ ! 2D fields:
!~ !  DO i = 1, myDim_nod2D
!~ !     GloPCO2surf(i) = state_p(i + offset(id% pCO2s)) ! diagnostic field (not distributed to the model. See int_recom/recom/gasx.F90)
!~ !     GloCO2flux(i)  = state_p(i + offset(id% CO2f )) ! diagnostic field (not distributed to the model. See int_recom/recom/gasx.F90)
!~ !     alphaCO2 and PistonVel are also diagnostic fields.
!~ !  END DO
 
!~ ! Initialize external nodes:
!~   call exchange_nod( tr_arr(:,:, 8) )
!~   call exchange_nod( tr_arr(:,:,17) )
!~   call exchange_nod( tr_arr(:,:, 4) )
!~   call exchange_nod( tr_arr(:,:,14) )
!~   call exchange_nod( tr_arr(:,:, 5) )
!~   call exchange_nod( tr_arr(:,:, 3) )
!~   call exchange_nod( tr_arr(:,:,13) )
!~   call exchange_nod( tr_arr(:,:,24) )
!~   call exchange_nod( tr_arr(:,:, 6) )
!~   call exchange_nod( tr_arr(:,:, 7) )
!~   call exchange_nod( tr_arr(:,:,15) )
!~   call exchange_nod( tr_arr(:,:,16) )
!~   call exchange_nod( tr_arr(:,:,18) )
!~ !  call exchange_nod( PAR3D          )
!~   call exchange_nod( tr_arr(:,:,12) )
!~   call exchange_nod( tr_arr(:,:,10) )
!~ !  call exchange_nod( diags3D(:,:,1) )
!~ !  call exchange_nod( diags3D(:,:,2) )
!~   call exchange_nod( tr_arr(:,:,22) )
!~ !  call exchange_nod( GloPCO2surf )
!~ !  call exchange_nod( GloCO2flux  )

  ! clean up:
  deallocate(U_node_upd,U_elem_upd)
 
  ENDIF
END SUBROUTINE distribute_state_pdaf
