subroutine recom_sinking_new(tr_num,mesh)

  use REcoM_declarations
  use REcoM_LocVar
  use REcoM_GloVar
  use recom_config
  use REcoM_ciso

  use g_clock
  use o_PARAM
  use g_PARSUP
  use g_rotate_grid
  use g_config
  use mod_MESH
  use i_arrays          ! a_ice, m_ice 
  use o_param           ! num_tracers
  use i_param
  use o_arrays
  use g_forcing_arrays  ! press_air
  use g_comm_auto
  use i_therm_param
  use g_comm
  use g_support
#ifdef use_PDAF
  use mod_carbon_fluxes_diags
  use mod_assim_pdaf, only: nlmax
  use mod_parallel_pdaf, only: mype_filter, abort_parallel
#endif
  
  
  implicit none

  type(t_mesh), intent(in) , target :: mesh

  Integer                           :: node, nz,  id, nzmin, nzmax, n,  tr_num, k, nlevels_nod2D_minimum
  Real(kind=WP)                     :: vd_flux(mesh%nl)
  Real(kind=8)                      :: dz_trr(mesh%nl), aux
  Real(kind=8)                      :: wLoc,wM,wPs
  Real(kind=8)                      :: Rjp,Rj,Rjm

  Real(kind=8)                      :: cfl, d0, d1, thetaP, thetaM, psiP, psiM
  Real(kind=8)                      :: onesixth	= 	1.d0/6.d0
  Real(kind=8)                      :: c1, c2
  Real(kind=8)                      :: Vsink, tv                        ! positive downward
  Real(kind=8),dimension(mesh%nl)   :: Wvel_flux

#include "../associate_mesh.h"

! Constant sinking velocities from namelist

    Vsink=0.0_WP

    if (tracer_id(tr_num)==1007 .or.    &  !idetn
        tracer_id(tr_num)==1008 .or.    &  !idetc
        tracer_id(tr_num)==1017 .or.    &  !idetsi
        tracer_id(tr_num)==1021 ) then     !idetcal
	   
            Vsink = VDet       ! [m/day; positive downward]
          
    elseif(tracer_id(tr_num)==1004 .or. &  !iphyn
           tracer_id(tr_num)==1005 .or. &  !iphyc
           tracer_id(tr_num)==1020 .or. &  !iphycal
           tracer_id(tr_num)==1006 ) then  !ipchl

            Vsink = VPhy       ! [m/day; positive downward]

    elseif(tracer_id(tr_num)==1013 .or. &  !idian
           tracer_id(tr_num)==1014 .or. &  !idiac
           tracer_id(tr_num)==1016 .or. &  !idiasi
           tracer_id(tr_num)==1015 ) then  !idchl

            Vsink = VDia       ! [m/day; positive downward]
            
    elseif(tracer_id(tr_num)==1025 .or. &  !idetz2n
           tracer_id(tr_num)==1026 .or. &  !idetz2c
           tracer_id(tr_num)==1027 .or. &  !idetz2si
           tracer_id(tr_num)==1028 ) then  !idetz2calc 
                
            Vsink = VDet_zoo2  ! [m/day; positive downward]
    end if
    

if (Vsink .gt. 0.1) then ! No sinking if Vsink < 0.1 m/day

   vd_flux = 0.0d0

   do n = 1,myDim_nod2D
      if (ulevels_nod2D(n)>1) cycle 
      nzmin = ulevels_nod2D(n)
      nzmax = nlevels_nod2D(n)-1

      ! distance between tracer points, surface and bottom dz_trr is half 
      ! the layerthickness
      dz_trr                = 0.0d0
      dz_trr(nzmin+1:nzmax) = abs(Z_3d_n(nzmin:nzmax-1,n)-Z_3d_n(nzmin+1:nzmax,n))
      dz_trr(nzmin)         = hnode(nzmin,n)/2.0d0
      dz_trr(nzmax+1)       = hnode(nzmax,n)/2.0d0

      Wvel_flux(nzmin:nzmax+1)= 0.d0  ! Vertical velocity for BCG tracers
                                                                     
      do nz=nzmin,nzmax+1
      
          ! Constant sinking for plankton and second detritus
          
            if(tracer_id(tr_num)==1025 .or. &  !idetz2n
               tracer_id(tr_num)==1026 .or. &  !idetz2c
               tracer_id(tr_num)==1027 .or. &  !idetz2si
               tracer_id(tr_num)==1028 ) then  !idetz2calc
               
               Wvel_flux(nz) = -VSink/SecondsPerDay ! [m/sec; negative downward]
            
            elseif(tracer_id(tr_num)==1004 .or. &  !iphyn
                   tracer_id(tr_num)==1005 .or. &  !iphyc
                   tracer_id(tr_num)==1020 .or. &  !iphycal
                   tracer_id(tr_num)==1006 ) then  !ipchl
          
               Wvel_flux(nz) = -VSink/SecondsPerDay ! [m/sec; negative downward]
          
            elseif(tracer_id(tr_num)==1013 .or. &  !idian
                   tracer_id(tr_num)==1014 .or. &  !idiac
                   tracer_id(tr_num)==1016 .or. &  !idiasi
                   tracer_id(tr_num)==1015 ) then  !idchl
          
               Wvel_flux(nz) = -VSink/SecondsPerDay ! [m/sec; negative downward]
            endif
            
          ! Variable sinking for first detritus
          
          if (tracer_id(tr_num)==1007 .or.    &  !idetn
              tracer_id(tr_num)==1008 .or.    &  !idetc
              tracer_id(tr_num)==1017 .or.    &  !idetsi
              tracer_id(tr_num)==1021 ) then     !idetcal
          
              if (allow_var_sinking) then
                 Wvel_flux(nz) = - (Vsink+Vdet_a * abs(zbar_3d_n(nz,n))) / SecondsPerDay   ! [m/sec; negative downward]
              else
                 Wvel_flux(nz) = -Vsink/SecondsPerDay                                      ! [m/sec; negative downward]
              endif
          endif

      end do

if (1) then ! 3rd Order DST Sceheme with flux limiting. This code comes from old recom

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      vd_flux(nzmin:nzmax+1)= 0.0_WP


      ! loop over shared depth levels
      !do nz=nlevels_nod2D_minimum-1, nzmin+1,-1

      do nz=nzmax, nzmin+1,-1

         Rjp = tr_arr(nz,n,tr_num)              - tr_arr(min(nz+1,nzmax),n,tr_num)
         Rj  = tr_arr(max(nzmin,nz-1),n,tr_num) - tr_arr(nz,n,tr_num) 
         Rjm = tr_arr(max(nzmin,nz-2),n,tr_num) - tr_arr(max(nzmin,nz-1),n,tr_num)

         cfl = abs(Wvel_flux(nz) * dt / dz_trr(nz)) !(Z_n(nz-1)-Z_n(nz)))       ! [m/s] * [s] * [1/m]

         wPs = Wvel_flux(nz) + abs(Wvel_flux(nz)) ! --> Positive vertical velocity
         wM  = Wvel_flux(nz) - abs(Wvel_flux(nz)) ! --> Negative vertical velocity

         d0 = (2.d0 - cfl)*(1.d0 - cfl)*onesixth
         d1 = (1.d0 - cfl*cfl)*onesixth
	
         thetaP = Rjm/(1.d-20+Rj)
         psiP = d0 + d1*thetaP
         psiP = max(0.d0, min(min(1.d0,psiP), &
            (1.d0-cfl)/(1.d-20+cfl)*thetaP))

         thetaM = Rjp/(1.d-20 + Rj)	
         psiM = d0 + d1*thetaM
         psiM = max(0.d0, min(min(1.d0,psiM), &
            (1.d0-cfl)/(1.d-20-cfl)*thetaM))

         tv= (0.5 * wPs * (tr_arr(nz,n,tr_num)              + psiM * Rj)+ &
	      0.5 * wM  * (tr_arr(max(nzmin,nz-1),n,tr_num) + psiP * Rj))
         vd_flux(nz)= - tv*area(nz,n)
         
#ifdef use_PDAF
        IF (nz<=nlmax) THEN
        ! mod_carbonfluxes_diags
        ! sum of tracers here; at begin of tracer loop in oce_ale_tracer.F90, sum is set to zero
        
        ! detritus ("export")
        if (tracer_id(tr_num) == 1008 .or.    &      ! idetc
            tracer_id(tr_num) == 1021 .or.    &      ! idetcal
            tracer_id(tr_num) == 1026 .or.    &      ! idetz2c
            tracer_id(tr_num) == 1028 ) then         ! idetz2calc
            
            t_export(nz,n) = t_export(nz,n) + (-tv)   ! we want positive upwards: inverse sign
        endif
        
        ! alive carbon biomass
        if (tracer_id(tr_num) == 1005 .or.    &   ! iphyc
            tracer_id(tr_num) == 1020 .or.    &   ! iphycal
            tracer_id(tr_num) == 1014 ) then      ! idiac
            
            t_sink_livingmatter(nz,n) = t_sink_livingmatter(nz,n) + (-tv)
        endif
       ENDIF ! (nz<=nlmax)
#endif
         
      end do
end if ! 3rd Order DST Sceheme with flux limiting

if (0) then ! simple upwind

      ! Surface flux
      vd_flux(nzmin)= 0.0_WP

      ! Bottom flux
      vd_flux(nzmax+1)= 0.0_WP

      k=nod_in_elem2D_num(n)
      ! Screening minimum depth in neigbouring nodes around node n
      nlevels_nod2D_minimum=minval(nlevels(nod_in_elem2D(1:k, n))-1)

      do nz=nzmin+1,nzmax !nlevels_nod2D_minimum-1
!         tv = tr_arr(nz,n,tr_num)                                ! simple scheme        - test1
!         tv = 0.5_WP*(tr_arr(nz-1,n,tr_num)+tr_arr(nz,n,tr_num)) ! consider both layers - test2  
!         tv = tv*Wvel_flux(nz) ! Wvel_flux is negative
         tv = - 0.5* & ! - test3
            (tr_arr(nz-1,n,tr_num)*(Wvel_flux(nz)-abs(Wvel_flux(nz))) + &
             tr_arr(nz  ,n,tr_num)*(Wvel_flux(nz)+abs(Wvel_flux(nz))))
         vd_flux(nz)= tv*area(nz,n)

#ifdef use_PDAF
        IF (nz<=nlmax) THEN
        ! mod_carbonfluxes_diags
        ! sum of tracers here; at begin of tracer loop in oce_ale_tracer.F90, sum is set to zero
        
        ! detritus ("export")
        if (tracer_id(tr_num) == 1008 .or.    &      ! idetc
            tracer_id(tr_num) == 1021 .or.    &      ! idetcal
            tracer_id(tr_num) == 1026 .or.    &      ! idetz2c
            tracer_id(tr_num) == 1028 ) then         ! idetz2calc
            
            t_export(nz,n) = t_export(nz,n) + (-tv)  ! we want positive upwards: inverse sign
        endif
        
        ! alive carbon biomass
        if (tracer_id(tr_num) == 1005 .or.    &   ! iphyc
            tracer_id(tr_num) == 1020 .or.    &   ! iphycal
            tracer_id(tr_num) == 1014 ) then      ! idiac
            
            t_sink_livingmatter(nz,n) = t_sink_livingmatter(nz,n) + (-tv)
        endif
       ENDIF ! (nz<=nlmax)
#endif

      end do

end if ! simple upwind

      do nz=nzmin,nzmax

         ! vert_sink is set to zero before recom_sinking is called, so this is just adding to zero?!
         vert_sink(nz,n) = vert_sink(nz,n) + (vd_flux(nz)-vd_flux(nz+1))*dt/areasvol(nz,n)/hnode_new(nz,n) !/(zbar_3d_n(nz,n)-zbar_3d_n(nz+1,n))
         
#ifdef use_PDAF
        IF (nz<=nlmax) THEN
        ! mod_carbonfluxes_diags
        ! sum of tracers here; at begin of tracer loop in oce_ale_tracer.F90, sum is set to zero
        
        ! detritus ("export")
        if (tracer_id(tr_num) == 1008 .or.    &      ! idetc
            tracer_id(tr_num) == 1021 .or.    &      ! idetcal
            tracer_id(tr_num) == 1026 .or.    &      ! idetz2c
            tracer_id(tr_num) == 1028 ) then         ! idetz2calc
            
            if (cfconc) then
               ! concentration
               s_export(nz,n) = s_export(nz,n) + vert_sink(nz,n)/dt
            else
               ! mass
               s_export(nz,n) = s_export(nz,n) + vert_sink(nz,n)/dt*areasvol(nz,n)*hnode_new(nz,n)
            endif
            
        endif
        
        ! alive carbon biomass
        if (tracer_id(tr_num) == 1005 .or.    &   ! iphyc
            tracer_id(tr_num) == 1020 .or.    &   ! iphycal
            tracer_id(tr_num) == 1014 ) then      ! idiac
            
            if (cfconc) then
               ! concentration
               s_sink_livingmatter(nz,n) = s_sink_livingmatter(nz,n) + vert_sink(nz,n)/dt
            else
               ! mass
               s_sink_livingmatter(nz,n) = s_sink_livingmatter(nz,n) + vert_sink(nz,n)/dt*areasvol(nz,n)*hnode_new(nz,n)
            endif
        endif
       ENDIF ! (nz<=nlmax)
#endif
         
      end do
      
      ! sum of detritus tracers export at 200m for state vector (at begin of tracer loop in oce_ale_tracer.F90, sum is zero)
        if (tracer_id(tr_num) == 1008 .or.    &   ! idetc
            tracer_id(tr_num) == 1021 .or.    &   ! idetcal
            tracer_id(tr_num) == 1026 .or.    &   ! idetz2c
            tracer_id(tr_num) == 1028 ) then      ! idetz2calc
            
            IF (area(16,n)>0) then 
              export(n) = export(n) + vd_flux(16)/area(16,n) ! nz=16 at depth of 190m
            endif
        endif
      
   end do ! do n = 1,myDim_nod2D
   
#ifdef use_PDAF
   if (cfdiags_debug) then
      vname_cfdiags = 's_export'
      call debug_vert(vname_cfdiags,s_export)
   endif
#endif

end if ! Vsink .gt. 0.1

end subroutine recom_sinking_new
