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
       ONLY: offset, mesh_fesom, id
  USE g_PARSUP, &
       ONLY: mydim_nod2d, myDim_elem2D
  USE o_arrays, &
       ONLY: eta_n, uv, wvel, tr_arr, unode, MLD1
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
           state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% PhyChl)) = tr_arr(k, i, 8) ! small phytoplankton chlorophyll
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DiaChl)) = tr_arr(k, i, 17) ! diatom chlorophyll
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DIC)) = tr_arr(k, i, 4) ! DIC
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DOC)) = tr_arr(k, i, 14) ! DOC
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% Alk)) = tr_arr(k, i, 5) ! Alk
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DIN)) = tr_arr(k, i, 3) ! DIN
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% DOC)) = tr_arr(k, i, 13) ! DON
           ! state_p((i-1) * (mesh_fesom%nl-1) + k + offset(id% O2)) = tr_arr(k, i, 24) ! O2
        ENDDO
        ! 2D-fields
        ! state_p(i + offset(id% pCO2s )) = GloPCO2surf(i) ! surface ocean partial pressure CO2
        ! state_p(i + offset(id% CO2f ))  = GloCO2flux(i)  ! CO2 flux (from atmosphere into ocean)
   ENDDO
   
   ! chlorophyll
   ! tr_arr(:,:,8)   --> small phytoplankton chlorophyll
   ! tr_arr(:,:,17)  --> diatom chlorophyll
   ! tr_arr(:,:,3)   --> DIN; dissolved inorganic nitrogen
   ! tr_arr(:,:,13)  --> DON; dissolved organic nitrogen
   ! tr_arr(:,:,5)   --> alkalinity
   ! tr_arr(:,:,4)   --> DIC
   ! tr_arr(:,:,14)  --> DOC
   ! tr_arr(:,:,24)  --> O2
   ! GloPCO2surf(myDim) --> pCO2 surface
   ! GloCO2flux (myDim) --> CO2 flux
   ! --> PAR
   ! --> sum: DIN + DON (?)
   ! --> total organic carbon: PhyC + DiaC + DetC + DOC + HetC
   ! --> "silicate?": DiaSi + DetSi + Si (silicic acid) ???
   
   
  
END SUBROUTINE collect_state_pdaf
