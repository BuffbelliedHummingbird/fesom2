SUBROUTINE collect_state_pdaf(dim_p, state_p)

! DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEEK/EnKF/SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! This subroutine is called during the forecast 
! phase from PDAF\_put\_state\_X after the 
! propagation of each ensemble member. 
! The supplied state vector has to be initialized
! from the model fields (typically via a module). 
! With parallelization, MPI communication might be 
! required to initialize state vectors for all 
! subdomains on the model PEs. 
!
! The routine is executed by each process that is
! participating in the model integrations.
!
! REVISION HISTORY:
! 2004-11 - Lars Nerger  - Initial code
! 2019-11 - Longjiang Mu - Initial commit for AWI-CM3
! 2022-03 - Frauke B     - Removed sea-ice, added velocities
!
! !USES:
  USE mod_parallel_pdaf, &
       ONLY: mype_world
  USE mod_assim_pdaf, &
       ONLY: offset, mesh_fesom, id, dim_fields
  USE g_PARSUP, &
       ONLY: mydim_nod2d, myDim_elem2D
  USE o_arrays, &
       ONLY: eta_n, uv, wvel, tr_arr, unode, MLD1
  USE REcoM_GloVar, &
       ONLY: GloPCO2surf, GloCO2flux, Diags3D, PAR3D, export, PistonVelocity, alphaCO2
  USE i_arrays, &
       ONLY: a_ice
  USE mod_parallel_pdaf, &
       ONLY: task_id, mype_model
  USE g_clock, &
       ONLY: daynew


  IMPLICIT NONE
  
! ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! Local state vector

! CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
! Called by: init_ens_PDAF

! Local variables
  INTEGER :: i, k, b         ! Counter
  
! Debugging:
  CHARACTER(len=5)   :: mype_string
  CHARACTER(len=3)   :: day_string

! *************************************************
! *** Initialize state vector from model fields ***
! *************************************************

   ! *** Dimensions of FESOM arrays:
   ! * eta_n          (myDim_nod2D + eDim_nod2D)            ! Sea Surface Height
   ! * UV(1,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity u
   ! * UV(2,:,:)      (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity v
   ! * wvel           (nl,   myDim_nod2D + eDim_nod2D)      ! Velocity w
   ! * tr_arr(:,:,1)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Temperature
   ! * tr_arr(:,:,2)  (nl-1, myDim_nod2D + eDim_nod2D, 1)   ! Salinity
   ! * unode(1,:,:)   (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity u interpolated on nodes
   ! * unode(2,:,:)   (1, nl-1, myDim_elem2D + eDim_elem2D) ! Velocity v interpolated on nodes
   ! * a_ice          (myDim_nod2D + eDim_nod2D)            ! Sea-ice concentration
   ! ***

   ! SSH (1)
   DO i = 1, myDim_nod2D
      state_p(i + offset(id% SSH)) = eta_n(i)
   END DO

   ! u (2) and v (3) velocities (interpolated on nodes)
   
   DO i = 1, myDim_nod2D
      DO k =1, mesh_fesom%nl-1
         state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% u)) = Unode(1, k, i) ! u
         state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% v)) = Unode(2, k, i) ! v
      END DO
   END DO
   
   ! w velocity (4)
   DO i = 1, myDim_nod2D
      DO k = 1, mesh_fesom%nl
         state_p((i-1) * mesh_fesom%nl + k + offset(id% w)) = wvel(k, i)
      END DO
   END DO
   
   ! temp (5) and salt (6)
   DO i = 1, myDim_nod2D
      DO k = 1, mesh_fesom%nl-1
         state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% temp)) = tr_arr(k, i, 1) ! T
         state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% salt)) = tr_arr(k, i, 2) ! S
      END DO
   END DO
   
   ! sea-ice concentration (7) and MLD1 (8)
   DO i = 1, myDim_nod2D
      state_p(i + offset(id% a_ice)) = a_ice(i)
      state_p(i + offset(id% MLD1 )) = MLD1(i)
   END DO
   
   ! biogeochemical tracers
   DO i = 1, myDim_nod2D
        ! 3D-fields
        DO k = 1, mesh_fesom%nl-1
           ! nanophytoplankton/calcifiers:
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyChl)) = tr_arr(k, i,  8) ! small phytoplankton chlorophyll
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyN))   = tr_arr(k, i,  6) ! intracellular conc of nitrogen in small phytoplankton
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyC))   = tr_arr(k, i,  7) ! intracellular conc of carbon in small phytoplankton
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyCalc))= tr_arr(k, i, 22) ! small phytoplankton CaCO2
           ! diatoms:
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaChl)) = tr_arr(k, i, 17) ! diatom chlorophyll
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaN))   = tr_arr(k, i, 15) ! intracellular conc of nitrogen in diatoms
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaC))   = tr_arr(k, i, 16) ! intracellular conc of carbon in diatoms
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaSi))  = tr_arr(k, i, 18) ! intracellular conc of Si in diatoms
           ! small, fast-growing zooplankton
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo1C))   = tr_arr(k, i, 12) ! intracellular conc of carbon in zooplankton 1
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo1N))   = tr_arr(k, i, 11) ! intracellular conc of nitrogen in zooplankton 1
           ! macrozooplankton/antarctic krill:
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo2C))   = tr_arr(k, i, 26) ! intracellular conc of carbon in zooplankton 2
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Zo2N))   = tr_arr(k, i, 27) ! intracellular conc of nitrogen in zooplankton 2
           ! dissolved tracers
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DIC))    = tr_arr(k, i,  4) ! dissolved inorganic carbon
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DOC))    = tr_arr(k, i, 14) ! dissolved organic carbon
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Alk))    = tr_arr(k, i,  5) ! alkalinity
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DIN))    = tr_arr(k, i,  3) ! dissolved inorganic nitrogen
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DON))    = tr_arr(k, i, 13) ! dissolved organic nitrogen
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% O2))     = tr_arr(k, i, 24) ! oxygen
           ! others:
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PAR))    = PAR3D(k, i)      ! photosynthetically active radiation
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DetC))   = tr_arr(k, i, 10) ! detritus carbon
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% NPPn))   = diags3D(k, i, 1) ! net primary production small phytoplankton
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% NPPd))   = diags3D(k, i, 2) ! net primary production diatoms
           
           
!~                tr_arr(:,:,25) = tiny                   ! tracer 25 = Zoo2N
!~     tr_arr(:,:,26) = tiny * Redfield        ! tracer 26 = Zoo2C
!~     tr_arr(:,:,27) = tiny                   ! tracer 27 = DetZ2N                              
!~     tr_arr(:,:,28) = tiny                   ! tracer 28 = DetZ2C                                    
!~     tr_arr(:,:,29) = tiny                   ! tracer 29 = DetZ2Si                            
!~     tr_arr(:,:,30) = tiny                   ! tracer 30 = DetZ2Calc 
        ENDDO
        ! 2D-fields
        state_p(i + offset(id% pCO2s ))    = GloPCO2surf(i)    ! surface ocean partial pressure CO2
        state_p(i + offset(id% CO2f ))     = GloCO2flux(i)     ! CO2 flux (from atmosphere into ocean)
        state_p(i + offset(id% export))    = export(i)         ! Export through particle sinking
        state_p(i + offset(id% alphaCO2))  = alphaCO2(i)       ! Solubility of CO2
        state_p(i + offset(id% PistonVel)) = PistonVelocity(i) ! Air-sea gas transfer velocity

   ENDDO
   
!~    ! derived fields:
!~    state_p(offset(id% TChl   )+1 : offset(id% TChl   )+dim_fields(id% TChl   )) = & ! Total chlorophyll
!~    state_p(offset(id% PhyChl )+1 : offset(id% PhyChl )+dim_fields(id% PhyChl )) + &
!~    state_p(offset(id% DiaChl )+1 : offset(id% DiaChl )+dim_fields(id% DiaChl ))
   
!~    state_p(offset(id% TDN ) +1 : offset(id% TDN )+dim_fields(id% TDN )) = & ! Total dissolved nitrogen
!~    state_p(offset(id% DIN ) +1 : offset(id% DIN )+dim_fields(id% DIN )) + &
!~    state_p(offset(id% DON ) +1 : offset(id% DON )+dim_fields(id% DON ))
   
!~    state_p(offset(id% TOC  ) +1 : offset(id% TOC  )+dim_fields(id% TOC  )) = & ! Total organic carbon
!~    state_p(offset(id% PhyC ) +1 : offset(id% PhyC )+dim_fields(id% PhyC )) + &
!~    state_p(offset(id% DiaC ) +1 : offset(id% DiaC )+dim_fields(id% DiaC )) + &
!~    state_p(offset(id% HetC ) +1 : offset(id% HetC )+dim_fields(id% HetC )) + &
!~    state_p(offset(id% DetC ) +1 : offset(id% DetC )+dim_fields(id% DetC )) + &
!~    state_p(offset(id% DOC  ) +1 : offset(id% DOC  )+dim_fields(id% DOC  ))
   
   
  
END SUBROUTINE collect_state_pdaf
