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
       ONLY: offset, loc_radius,this_is_pdaf_restart, mesh_fesom, nlmax, &
       dim_fields, id, istep_asml, step_null, start_from_ENS_spinup
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
  
!~   DO i = 1, myDim_nod2D
!~    DO k = 1, nlmax
!~       ! T
!~       s = (i-1) * (nlmax) + k + offset(id% temp)
!~       tr_arr(k, i, 1) = state_p(s)
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%temp)%variable, s, tr_arr(k, i, 1), state_p(s)
!~       ! S
!~       s = (i-1) * (nlmax) + k + offset(id% salt)
!~       tr_arr(k, i, 2) = state_p(s)
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id%salt)%variable, s, tr_arr(k, i, 2), state_p(s)
!~    END DO
!~   END DO
  
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
!~   call exchange_nod(tr_arr(:,:,1))    ! Temp
!~   call exchange_nod(tr_arr(:,:,2))    ! Salt


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
                            

! *********************************
! *** Biogeochemistry           ***
! *********************************

! 3D fields:
!~   DO i = 1, myDim_nod2D
!~    DO k = 1, nlmax
   
!~       if (write_debug) then
!~       ! small phytoplankton:
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% PhyChl )%variable,   (i-1) * (nlmax) + k + offset(id% PhyChl ),   tr_arr(k, i,  8),   state_p((i-1) * (nlmax) + k + offset(id% PhyChl ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% PhyN   )%variable,   (i-1) * (nlmax) + k + offset(id% PhyN   ),   tr_arr(k, i,  6),   state_p((i-1) * (nlmax) + k + offset(id% PhyN   ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% PhyC   )%variable,   (i-1) * (nlmax) + k + offset(id% PhyC   ),   tr_arr(k, i,  7),   state_p((i-1) * (nlmax) + k + offset(id% PhyC   ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% PhyCalc)%variable,   (i-1) * (nlmax) + k + offset(id% PhyCalc),   tr_arr(k, i, 22),   state_p((i-1) * (nlmax) + k + offset(id% PhyCalc))
!~       ! diatoms:                                                                                                                                                      
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DiaChl )%variable,   (i-1) * (nlmax) + k + offset(id% DiaChl ),   tr_arr(k, i, 17),   state_p((i-1) * (nlmax) + k + offset(id% DiaChl ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DiaN   )%variable,   (i-1) * (nlmax) + k + offset(id% DiaN   ),   tr_arr(k, i, 15),   state_p((i-1) * (nlmax) + k + offset(id% DiaN   ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DiaC   )%variable,   (i-1) * (nlmax) + k + offset(id% DiaC   ),   tr_arr(k, i, 16),   state_p((i-1) * (nlmax) + k + offset(id% DiaC   ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DiaSi  )%variable,   (i-1) * (nlmax) + k + offset(id% DiaSi  ),   tr_arr(k, i, 18),   state_p((i-1) * (nlmax) + k + offset(id% DiaSi  ))
!~       ! small, fast-growing zooplankton                                                                                                                           
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Zo1C   )%variable,   (i-1) * (nlmax) + k + offset(id% Zo1C)   ,   tr_arr(k, i, 12),   state_p((i-1) * (nlmax) + k + offset(id% Zo1C))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Zo1N   )%variable,   (i-1) * (nlmax) + k + offset(id% Zo1N)   ,   tr_arr(k, i, 11),   state_p((i-1) * (nlmax) + k + offset(id% Zo1N))
!~       ! macrozooplankton/antarctic krill:                                                                                                                                
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Zo2C   )%variable,   (i-1) * (nlmax) + k + offset(id% Zo2C)   ,   tr_arr(k, i, 26),   state_p((i-1) * (nlmax) + k + offset(id% Zo2C))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Zo2N   )%variable,   (i-1) * (nlmax) + k + offset(id% Zo2N)   ,   tr_arr(k, i, 25),   state_p((i-1) * (nlmax) + k + offset(id% Zo2N))
!~       ! dissolved tracer pools:                                                                                                                                           
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DIC    )%variable,   (i-1) * (nlmax) + k + offset(id% DIC    ),   tr_arr(k, i,  4),   state_p((i-1) * (nlmax) + k + offset(id% DIC    ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DOC    )%variable,   (i-1) * (nlmax) + k + offset(id% DOC    ),   tr_arr(k, i, 14),   state_p((i-1) * (nlmax) + k + offset(id% DOC    ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Alk    )%variable,   (i-1) * (nlmax) + k + offset(id% Alk    ),   tr_arr(k, i,  5),   state_p((i-1) * (nlmax) + k + offset(id% Alk    ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DIN    )%variable,   (i-1) * (nlmax) + k + offset(id% DIN    ),   tr_arr(k, i,  3),   state_p((i-1) * (nlmax) + k + offset(id% DIN    ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DON    )%variable,   (i-1) * (nlmax) + k + offset(id% DON    ),   tr_arr(k, i, 13),   state_p((i-1) * (nlmax) + k + offset(id% DON    ))
!~       write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% O2     )%variable,   (i-1) * (nlmax) + k + offset(id% O2     ),   tr_arr(k, i, 24),   state_p((i-1) * (nlmax) + k + offset(id% O2     ))
!~       ! detritus:                                                                                                                                                        
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DetC   )%variable,   (i-1) * (nlmax) + k + offset(id% DetC   ),   tr_arr(k, i, 10),   state_p((i-1) * (nlmax) + k + offset(id% DetC    ))
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DetCalc)%variable,   (i-1) * (nlmax) + k + offset(id% DetCalc),   tr_arr(k, i, 23),   state_p((i-1) * (nlmax) + k + offset(id% DetCalc ))
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DetSi  )%variable,   (i-1) * (nlmax) + k + offset(id% DetSi  ),   tr_arr(k, i, 19),   state_p((i-1) * (nlmax) + k + offset(id% DetSi   ))
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% DetN   )%variable,   (i-1) * (nlmax) + k + offset(id% DetN   ),   tr_arr(k, i,  9),   state_p((i-1) * (nlmax) + k + offset(id% DetN    ))
      
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Det2C   )%variable,   (i-1) * (nlmax) + k + offset(id% Det2C   ),   tr_arr(k, i, 28),   state_p((i-1) * (nlmax) + k + offset(id% Det2C    ))
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Det2Calc)%variable,   (i-1) * (nlmax) + k + offset(id% Det2Calc),   tr_arr(k, i, 30),   state_p((i-1) * (nlmax) + k + offset(id% Det2Calc ))
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Det2Si  )%variable,   (i-1) * (nlmax) + k + offset(id% Det2Si  ),   tr_arr(k, i, 29),   state_p((i-1) * (nlmax) + k + offset(id% Det2Si   ))
!~       if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6,G15.6)') sfields(id% Det2N   )%variable,   (i-1) * (nlmax) + k + offset(id% Det2N   ),   tr_arr(k, i, 27),   state_p((i-1) * (nlmax) + k + offset(id% Det2N    ))
!~       endif ! write_debug
   
!~       ! small phytoplankton:
!~       tr_arr(k, i,  8) = state_p((i-1) * (nlmax) + k + offset(id% PhyChl ))
!~       tr_arr(k, i,  6) = state_p((i-1) * (nlmax) + k + offset(id% PhyN   ))
!~       tr_arr(k, i,  7) = state_p((i-1) * (nlmax) + k + offset(id% PhyC   ))
!~       tr_arr(k, i, 22) = state_p((i-1) * (nlmax) + k + offset(id% PhyCalc))
!~       ! diatoms:
!~       tr_arr(k, i, 17) = state_p((i-1) * (nlmax) + k + offset(id% DiaChl ))
!~       tr_arr(k, i, 15) = state_p((i-1) * (nlmax) + k + offset(id% DiaN   ))
!~       tr_arr(k, i, 16) = state_p((i-1) * (nlmax) + k + offset(id% DiaC   ))
!~       tr_arr(k, i, 18) = state_p((i-1) * (nlmax) + k + offset(id% DiaSi  ))
!~       ! small, fast-growing zooplankton
!~       tr_arr(k, i, 12) = state_p((i-1) * (nlmax) + k + offset(id% Zo1C)) ! intracellular conc of carbon in zooplankton 1
!~       tr_arr(k, i, 11) = state_p((i-1) * (nlmax) + k + offset(id% Zo1N)) ! intracellular conc of nitrogen in zooplankton 1
!~       ! macrozooplankton/antarctic krill:
!~       tr_arr(k, i, 26) = state_p((i-1) * (nlmax) + k + offset(id% Zo2C)) ! intracellular conc of carbon in zooplankton 2
!~       tr_arr(k, i, 25) = state_p((i-1) * (nlmax) + k + offset(id% Zo2N)) ! intracellular conc of nitrogen in zooplankton 2
!~       ! dissolved tracer pools:
!~       tr_arr(k, i,  4) = state_p((i-1) * (nlmax) + k + offset(id% DIC    ))
!~       tr_arr(k, i, 14) = state_p((i-1) * (nlmax) + k + offset(id% DOC    ))
!~       tr_arr(k, i,  5) = state_p((i-1) * (nlmax) + k + offset(id% Alk    ))
!~       tr_arr(k, i,  3) = state_p((i-1) * (nlmax) + k + offset(id% DIN    ))
!~       tr_arr(k, i, 13) = state_p((i-1) * (nlmax) + k + offset(id% DON    ))
!~       tr_arr(k, i, 24) = state_p((i-1) * (nlmax) + k + offset(id% O2     ))
!~       ! detritus:
!~       tr_arr(k, i, 10) = state_p((i-1) * (nlmax) + k + offset(id% DetC   ))
!~       tr_arr(k, i, 23) = state_p((i-1) * (nlmax) + k + offset(id% DetCalc))
!~       tr_arr(k, i, 19) = state_p((i-1) * (nlmax) + k + offset(id% DetSi  ))
!~       tr_arr(k, i,  9) = state_p((i-1) * (nlmax) + k + offset(id% DetN   ))
      
!~       tr_arr(k, i, 28) = state_p((i-1) * (nlmax) + k + offset(id% Det2C   ))
!~       tr_arr(k, i, 30) = state_p((i-1) * (nlmax) + k + offset(id% Det2Calc))
!~       tr_arr(k, i, 29) = state_p((i-1) * (nlmax) + k + offset(id% Det2Si  ))
!~       tr_arr(k, i, 27) = state_p((i-1) * (nlmax) + k + offset(id% Det2N   ))
      
      ! diagnostic fields:
!     PAR3D (k, i)     = state_p((i-1) * (nlmax) + k + offset(id% PAR    )) ! diagnostic field (not distributed to the model. See int_recom/recom_sms.F90)
!     diags3D(k, i, 1) = state_p((i-1) * (nlmax) + k + offset(id% NPPn   )) ! diagnostic field (not distributed to the model. See int_recom/recom_sms.F90)
!     diags3D(k, i, 2) = state_p((i-1) * (nlmax) + k + offset(id% NPPd   )) ! diagnostic field (not distributed to the model. See int_recom/recom_sms.F90)      
      
!~    END DO
!~   END DO   
           
! 2D fields:
!  DO i = 1, myDim_nod2D
!     GloPCO2surf(i) = state_p(i + offset(id% pCO2s)) ! diagnostic field (not distributed to the model. See int_recom/recom/gasx.F90)
!     GloCO2flux(i)  = state_p(i + offset(id% CO2f )) ! diagnostic field (not distributed to the model. See int_recom/recom/gasx.F90)
!     alphaCO2 and PistonVel are also diagnostic fields.
!  END DO
 
! Initialize external nodes:
!~   call exchange_nod( tr_arr(:,:, 8) ) ! small phytoplankton
!~   call exchange_nod( tr_arr(:,:, 6) )
!~   call exchange_nod( tr_arr(:,:, 7) )
!~   call exchange_nod( tr_arr(:,:,22) )
!~   call exchange_nod( tr_arr(:,:,17) ) ! diatoms
!~   call exchange_nod( tr_arr(:,:,15) )
!~   call exchange_nod( tr_arr(:,:,16) )
!~   call exchange_nod( tr_arr(:,:,18) )
!~   call exchange_nod( tr_arr(:,:,12) ) ! zooplankton 1
!~   call exchange_nod( tr_arr(:,:,11) )
!~   call exchange_nod( tr_arr(:,:,26) ) ! zooplankton 2
!~   call exchange_nod( tr_arr(:,:,25) )
!~   call exchange_nod( tr_arr(:,:, 4) ) ! dissolved tracers
!~   call exchange_nod( tr_arr(:,:,14) )
!~   call exchange_nod( tr_arr(:,:, 5) )
!~   call exchange_nod( tr_arr(:,:, 3) )
!~   call exchange_nod( tr_arr(:,:,13) )
!~   call exchange_nod( tr_arr(:,:,24) )
!~   call exchange_nod( tr_arr(:,:,10) ) ! detritus
!~   call exchange_nod( tr_arr(:,:,23) )
!~   call exchange_nod( tr_arr(:,:,19) )
!~   call exchange_nod( tr_arr(:,:, 9) )
!~   call exchange_nod( tr_arr(:,:,28) )
!~   call exchange_nod( tr_arr(:,:,30) )
!~   call exchange_nod( tr_arr(:,:,29) )
!~   call exchange_nod( tr_arr(:,:,27) )
!  call exchange_nod( PAR3D          )
!  call exchange_nod( GloPCO2surf )
!  call exchange_nod( GloCO2flux  )
!  call exchange_nod( diags3D(:,:,1) )
!  call exchange_nod( diags3D(:,:,2) )

  ! clean up:
  if (write_debug) close(fileID_debug)
  deallocate(U_node_upd,U_elem_upd)
 
  ENDIF
END SUBROUTINE distribute_state_pdaf
