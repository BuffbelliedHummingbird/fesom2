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
       ONLY: offset, mesh_fesom, nlmax, id, dim_fields, dim_state_p
  USE mod_nc_out_variables, &
       ONLY: sfields, ids_tr3D, nfields_tr3D
  USE g_PARSUP, &
       ONLY: mydim_nod2d, myDim_elem2D
  USE o_arrays, &
       ONLY: eta_n, uv, wvel, tr_arr, unode, MLD1, MLD2
  USE REcoM_GloVar, &
       ONLY: GloPCO2surf, GloCO2flux, Diags3D, PAR3D, export, PistonVelocity, alphaCO2
  USE i_arrays, &
       ONLY: a_ice
  USE mod_parallel_pdaf, &
       ONLY: task_id, mype_model
  USE g_clock, &
       ONLY: daynew, timenew
  USE PDAF_mod_filter, &
       ONLY: assim_flag


  IMPLICIT NONE
  
! ARGUMENTS:
  INTEGER, INTENT(in) :: dim_p           ! PE-local state dimension
  REAL, INTENT(inout) :: state_p(dim_p)  ! Local state vector

! CALLING SEQUENCE:
! Called by: PDAF_put_state_X   (as U_coll_state)
! Called by: init_ens_PDAF

! Local variables
  INTEGER :: i, k, b, s, istate, ifesom   ! Counters
  
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
        IF (assim_flag==0) THEN
          write_debug = .false.
        ELSE
          write_debug = .true.
        ENDIF
     ENDIF
  ENDIF
    
  IF (write_debug) THEN
         ! print state vector
         WRITE(day_string, '(i3.3)') daynew
         WRITE(tim_string, '(i5.5)') int(timenew)
         fileID_debug=10
         open(unit=fileID_debug, file='collect_state_pdaf_'//day_string//'_'//tim_string//'.txt', status='unknown')
  ENDIF

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
      s = i + offset(id% SSH)
      state_p(s) = eta_n(i)
      IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%SSH)%variable, s, eta_n(i)
   END DO

   ! u (2) and v (3) velocities (interpolated on nodes)
   
   DO i = 1, myDim_nod2D
      DO k =1, nlmax
         ! u
         s = (i-1) * (nlmax) + k + offset(id% u)
         state_p(s) = Unode(1, k, i)
         IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%u)%variable, s, Unode(1, k, i)
         ! v
         s = (i-1) * (nlmax) + k + offset(id% v)
         state_p(s) = Unode(2, k, i)
         IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%v)%variable, s, Unode(2, k, i)
      END DO
   END DO
   
   ! w velocity (4)
   DO i = 1, myDim_nod2D
      DO k = 1, nlmax
         s = (i-1) * nlmax + k + offset(id% w)
         state_p(s) = wvel(k, i)
         IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%w)%variable, s, wvel(k, i)
      END DO
   END DO
   
!~    ! temp (5) and salt (6)
!~    DO i = 1, myDim_nod2D
!~       DO k = 1, nlmax
!~          ! T
!~          s = (i-1) * (nlmax) + k + offset(id% temp)
!~          state_p(s) = tr_arr(k, i, 1) ! T
!~          IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%temp)%variable, s, tr_arr(k, i, 1)
!~          ! S
!~          s = (i-1) * (nlmax) + k + offset(id% salt)
!~          state_p(s) = tr_arr(k, i, 2) ! S
!~          IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%salt)%variable, s, tr_arr(k, i, 2)
!~       END DO
!~    END DO
   
   ! sea-ice concentration (7) and MLD1 (8 and 9)
   DO i = 1, myDim_nod2D
      ! a_ice
      s = i + offset(id% a_ice)
      state_p(s) = a_ice(i)
      IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%a_ice)%variable, s, a_ice(i)
      ! MLD1
      s = i + offset(id% MLD1)
      state_p(s) = MLD1(i)
      IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%MLD1)%variable, s, MLD1(i)
      ! MLD2
      s = i + offset(id% MLD2)
      state_p(s) = MLD2(i)
      IF (write_debug) WRITE(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%MLD2)%variable, s, MLD2(i)
   END DO
   
   ! ________________
   ! model 3D tracers
   ! ________________
   DO i = 1, myDim_nod2D
        DO k = 1, nlmax
        
           ! indeces to flatten 3D arrays
           s = (i-1) * (nlmax) + k
           
           ! tracer field loop
           DO b = 1, nfields_tr3D
           
              istate = ids_tr3D(b)                 ! index of field in state vector
              ifesom = sfields(istate)%trnumfesom  ! index of field in model tracer array
           
              state_p(s + offset(istate)) = tr_arr(k, i, ifesom)
           
           ENDDO
           
!~            ! nanophytoplankton/calcifiers:
!~            state_p(s + offset(id% PhyChl)) = tr_arr(k, i,  8) ! small phytoplankton chlorophyll
!~            state_p(s + offset(id% PhyN))   = tr_arr(k, i,  6) ! intracellular conc of nitrogen in small phytoplankton
!~            state_p(s + offset(id% PhyC))   = tr_arr(k, i,  7) ! intracellular conc of carbon in small phytoplankton
!~            state_p(s + offset(id% PhyCalc))= tr_arr(k, i, 22) ! small phytoplankton CaCO2
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PhyChl )%variable, s + offset(id% PhyChl) , tr_arr(k, i,  8)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PhyN   )%variable, s + offset(id% PhyN)   , tr_arr(k, i,  6)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PhyC   )%variable, s + offset(id% PhyC)   , tr_arr(k, i,  7)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PhyCalc)%variable, s + offset(id% PhyCalc), tr_arr(k, i, 22)
           
!~            ! diatoms:
!~            state_p(s + offset(id% DiaChl)) = tr_arr(k, i, 17) ! diatom chlorophyll
!~            state_p(s + offset(id% DiaN))   = tr_arr(k, i, 15) ! intracellular conc of nitrogen in diatoms
!~            state_p(s + offset(id% DiaC))   = tr_arr(k, i, 16) ! intracellular conc of carbon in diatoms
!~            state_p(s + offset(id% DiaSi))  = tr_arr(k, i, 18) ! intracellular conc of Si in diatoms
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DiaChl )%variable, s + offset(id% DiaChl), tr_arr(k, i, 17)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DiaN   )%variable, s + offset(id% DiaN)  , tr_arr(k, i, 15)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DiaC   )%variable, s + offset(id% DiaC)  , tr_arr(k, i, 16)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DiaSi  )%variable, s + offset(id% DiaSi) , tr_arr(k, i, 18)
           
!~            ! small, fast-growing zooplankton
!~            state_p(s + offset(id% Zo1C))   = tr_arr(k, i, 12) ! intracellular conc of carbon in zooplankton 1
!~            state_p(s + offset(id% Zo1N))   = tr_arr(k, i, 11) ! intracellular conc of nitrogen in zooplankton 1
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Zo1C )%variable, s + offset(id% Zo1C), tr_arr(k, i, 12)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Zo1N )%variable, s + offset(id% Zo1N), tr_arr(k, i, 11)
           
!~            ! macrozooplankton/antarctic krill:
!~            state_p(s + offset(id% Zo2C))   = tr_arr(k, i, 26) ! intracellular conc of carbon in zooplankton 2
!~            state_p(s + offset(id% Zo2N))   = tr_arr(k, i, 25) ! intracellular conc of nitrogen in zooplankton 2
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Zo2C )%variable, s + offset(id% Zo2C), tr_arr(k, i, 26)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Zo2N )%variable, s + offset(id% Zo2N), tr_arr(k, i, 25)
           
!~            ! dissolved tracers
!~            state_p(s + offset(id% DIC))    = tr_arr(k, i,  4) ! dissolved inorganic carbon
!~            state_p(s + offset(id% DOC))    = tr_arr(k, i, 14) ! dissolved organic carbon
!~            state_p(s + offset(id% Alk))    = tr_arr(k, i,  5) ! alkalinity
!~            state_p(s + offset(id% DIN))    = tr_arr(k, i,  3) ! dissolved inorganic nitrogen
!~            state_p(s + offset(id% DON))    = tr_arr(k, i, 13) ! dissolved organic nitrogen
!~            state_p(s + offset(id% O2))     = tr_arr(k, i, 24) ! oxygen
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DIC )%variable, s + offset(id% DIC), tr_arr(k, i,  4)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DOC )%variable, s + offset(id% DOC), tr_arr(k, i, 14)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Alk )%variable, s + offset(id% Alk), tr_arr(k, i,  5)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DIN )%variable, s + offset(id% DIN), tr_arr(k, i,  3)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DON )%variable, s + offset(id% DON), tr_arr(k, i, 13)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%O2  )%variable, s + offset(id% O2) , tr_arr(k, i, 24)
          
!~            ! detritus
!~            state_p(s + offset(id% DetC))   = tr_arr(k, i, 10) ! detritus carbon
!~            state_p(s + offset(id% DetCalc))= tr_arr(k, i, 23) ! detritus calcite
!~            state_p(s + offset(id% DetSi))  = tr_arr(k, i, 19) ! detritus silicate
!~            state_p(s + offset(id% DetN))   = tr_arr(k, i,  9) ! detritus nitrogen
           
!~            state_p(s + offset(id% Det2C))   = tr_arr(k, i, 28)
!~            state_p(s + offset(id% Det2Calc))= tr_arr(k, i, 30)
!~            state_p(s + offset(id% Det2Si))  = tr_arr(k, i, 29)
!~            state_p(s + offset(id% Det2N))   = tr_arr(k, i, 27)
           
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DetC    )%variable, s + offset(id% DetC   ), tr_arr(k, i, 10)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DetCalc )%variable, s + offset(id% DetCalc), tr_arr(k, i, 23)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DetSi   )%variable, s + offset(id% DetSi  ), tr_arr(k, i, 19)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%DetN    )%variable, s + offset(id% DetN   ), tr_arr(k, i,  9)
           
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Det2C    )%variable, s + offset(id% Det2C   ), tr_arr(k, i, 28)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Det2Calc )%variable, s + offset(id% Det2Calc), tr_arr(k, i, 30)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Det2Si   )%variable, s + offset(id% Det2Si  ), tr_arr(k, i, 29)
!~            if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%Det2N    )%variable, s + offset(id% Det2N   ), tr_arr(k, i, 27)
          
           ! others biogeochem 3D fields
           state_p(s + offset(id% PAR))    = PAR3D(k, i)      ! photosynthetically active radiation
           state_p(s + offset(id% NPPn))   = diags3D(k, i, 1) ! net primary production small phytoplankton
           state_p(s + offset(id% NPPd))   = diags3D(k, i, 2) ! net primary production diatoms
           
           if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PAR  )%variable, s + offset(id% PAR) , PAR3D(k, i)     
           if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%NPPn )%variable, s + offset(id% NPPn), diags3D(k, i, 1)
           if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%NPPd )%variable, s + offset(id% NPPd), diags3D(k, i, 2)

        ENDDO
   ENDDO
   
   ! biogeochem 2D fields
   DO i = 1, myDim_nod2D
   
        state_p(i + offset(id% pCO2s ))    = GloPCO2surf(i)    ! surface ocean partial pressure CO2
        state_p(i + offset(id% CO2f ))     = GloCO2flux(i)     ! CO2 flux (from atmosphere into ocean)
        state_p(i + offset(id% export))    = export(i)         ! Export through particle sinking
        state_p(i + offset(id% alphaCO2))  = alphaCO2(i)       ! Solubility of CO2
        state_p(i + offset(id% PistonVel)) = PistonVelocity(i) ! Air-sea gas transfer velocity

        if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%pCO2s    )%variable, i + offset(id% pCO2s )   , GloPCO2surf(i)    
        if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%CO2f     )%variable, i + offset(id% CO2f )    , GloCO2flux(i)     
        if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%export   )%variable, i + offset(id% export)   , export(i)         
        if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%alphaCO2 )%variable, i + offset(id% alphaCO2) , alphaCO2(i)       
        if (write_debug) write(fileID_debug, '(a10,1x,i8,1x,G15.6)') sfields(id%PistonVel)%variable, i + offset(id% PistonVel), PistonVelocity(i) 

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
   
   ! Close debug-file
   IF (write_debug) close(fileID_debug)
  
END SUBROUTINE collect_state_pdaf
