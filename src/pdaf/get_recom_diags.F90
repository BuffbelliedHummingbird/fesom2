! initially written by REcoM group, adapted by  O. Gurses 02.03.2020
! computes diagnostic REcoM variables (CO2 flux and pCO2) before timestepping

subroutine get_recom_diags(mesh)

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
  use i_arrays 		! a_ice, m_ice 
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
#endif
  
  
  implicit none
  type(t_mesh), intent(in) , target :: mesh
! ======================================================================================
!! Depth information

!! zbar (depth of layers) and Z (mid depths of layers)
!! zbar is negative 
!! zbar(nl) allocate the array for storing the standard depths
!! Z(nl-1)  mid-depths of cells

!! max. number of levels at node n
!! nzmax = nlevels_nod2D(n)
!! u_ice and v_ice are at nodes
!! u_w, v_w are at nodes (interpolated from elements)
!! u_wind and v_wind are always at nodes
! ======================================================================================

  real(kind=8)               :: SW, Loc_slp
  integer                    :: tr_num
  integer                    :: nz, n, nzmin, nzmax  
  integer                    :: idiags

  real(kind=8)               :: Sali, net, net1, net2
  real (kind=8), allocatable :: Temp(:),  zr(:), PAR(:)
  real(kind=8),  allocatable :: C(:,:)
  character(len=2)           :: tr_num_name
#ifdef use_PDAF
  integer                    :: nlay
#endif
  
#include "../associate_mesh.h"

  allocate(Temp(nl-1), zr(nl-1) , PAR(nl-1))
  allocate(C(nl-1,bgc_num))


  if (.not. use_REcoM) return

! ======================================================================================
!************************* READ SURFACE BOUNDARY FILES *********************************			

  if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> Atm_input'//achar(27)//'[0m'
  
    call Atm_input(mesh)        !<  read surface atmospheric deposition for Fe, N, CO2
    call River_input(mesh)      !<  read riverine input
    call Erosion_input(mesh)    !<  read erosion input
  
  if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> bio_fluxes'//achar(27)//'[0m'
  
    call bio_fluxes(mesh)       !<  alkalinity restoring/ virtual flux is possible

! ======================================================================================
!********************************* LOOP STARTS *****************************************			

  do n=1, myDim_nod2D  ! needs exchange_nod in the end
!     if (ulevels_nod2D(n)>1) cycle 
!            nzmin = ulevels_nod2D(n)

     !!---- Number of vertical layers
     nzmax = nlevels_nod2D(n)-1

     !!---- This is needed for piston velocity 
     Loc_ice_conc = a_ice(n) 

     !!---- Mean sea level pressure 
     Loc_slp = press_air(n)

     !!---- Benthic layers
     LocBenthos(1:benthos_num) = Benthos(n,1:benthos_num)

     !!---- Local conc of [H+]-ions from last time time step. Stored in LocVar
     !!---- used as first guess for H+ conc.in subroutine CO2flux (provided by recom_init)
     Hplus = GloHplus(n)                                  

     !!---- Interpolated wind from atmospheric forcing 
     !!---- temporarily stored in module LocVar
     ULoc = sqrt(u_wind(n)**2+v_wind(n)**2)

     !!---- Atmospheric CO2 in LocVar                                                                        
     LocAtmCO2         = AtmCO2(month)   
     if (ciso) then
        LocAtmCO2_13 = AtmCO2_13(month)
        LocAtmCO2_14 = AtmCO2_14(month)
        r_atm_13 = LocAtmCO2_13(1) / LocAtmCO2(1)
        r_atm_14 = LocAtmCO2_14(1) / LocAtmCO2(1)
     end if

     !!---- Shortwave penetration
     SW = parFrac * shortwave(n)
     SW = SW * (1.d0 - a_ice(n))

     !!---- Temperature in water column
     Temp(1:nzmax) = tr_arr(1:nzmax, n, 1)

     !!---- Surface salinity
     Sali = tr_arr(1,       n, 2)

     !!---- Biogeochemical tracers
     C(1:nzmax,1:bgc_num) = tr_arr(1:nzmax, n, 3:num_tracers)             

     !!---- Depth of the nodes in the water column 
     zr(1:nzmax) = Z_3d_n(1:nzmax, n)                          

     !!---- The PAR in the local water column is initialized
     PAR(1:nzmax) = 0.d0                                        

     !!---- a_ice(row): Ice concentration in the local node
     FeDust = GloFeDust(n) * (1 - a_ice(n)) * dust_sol    
     NDust = GloNDust(n)  * (1 - a_ice(n))

     allocate(Diags3Dloc(nzmax,8))
     Diags3Dloc(:,:) = 0.d0
     
#ifdef use_PDAF
     nlay = min(nzmax,nlmax)
     ! initialize local carbon flux diags
     allocate(cffields(id_s_bio_dic)         %loc(nlay))
     allocate(cffields(id_s_bio_livingmatter)%loc(nlay))
     allocate(cffields(id_s_bio_deadmatter)  %loc(nlay))
     allocate(cffields(id_s_bio_alk)         %loc(nlay))
     cffields(id_s_bio_dic)         %loc = 0.d0
     cffields(id_s_bio_livingmatter)%loc = 0.d0
     cffields(id_s_bio_deadmatter)  %loc = 0.d0
     cffields(id_s_bio_alk)         %loc = 0.d0
#endif

if (recom_debug .and. mype==0) print *, achar(27)//'[36m'//'     --> REcoM_Forcing'//achar(27)//'[0m'


! ======================================================================================
!******************************** RECOM FORCING ****************************************

     call get_REcoM_Forcing(zr, n, nzmax, C, SW, Loc_slp, Temp, Sali, PAR, mesh) ! includes call to recom_sms

     !!---- Local variables that have been changed during the time-step are stored so they can be saved
     Benthos(n,1:benthos_num)     = LocBenthos(1:benthos_num)                      ! Updating Benthos values

     Diags2D(n,1:8)               = LocDiags2D(1:8)                                ! Updating diagnostics
     GloPCO2surf(n)               = pco2surf(1)
     GlodPCO2surf(n)              = dpco2surf(1)

     GloCO2flux(n)                = dflux(1)                              !  [mmol/m2/day]
     GloCO2flux_seaicemask(n)     = co2flux_seaicemask(1)                 !  [mmol/m2/s]
     GloO2flux_seaicemask(n)      = o2flux_seaicemask(1)                  !  [mmol/m2/s]
     if (ciso) then
        GloCO2flux_seaicemask_13(n)     = co2flux_seaicemask_13(1)        !  [mmol/m2/s]
        GloCO2flux_seaicemask_14(n)     = co2flux_seaicemask_14(1)        !  [mmol/m2/s]
     end if

     GloHplus(n)                  = ph(1) ! hplus

     AtmFeInput(n)                = FeDust
     AtmNInput(n)                 = NDust 

     PistonVelocity(n)         = kw660(1) ! CN: write Piston velocity                                
     alphaCO2(n)               = K0(1)  ! CN: write CO2 solubility

     GlodecayBenthos(n, 1:benthos_num) = decayBenthos(1:benthos_num)/SecondsPerDay ! convert from [mmol/m2/d] to [mmol/m2/s]  

     PAR3D(1:nzmax,n)             = PAR(1:nzmax) !     PAR3D(inds(1:nn))   = PAR(1:nn)
   
     do idiags = 1,diags3d_num
       Diags3D(1:nzmax,n,idiags)  = Diags3Dloc(1:nzmax,idiags) ! 1=NPPnano, 2=NPPdia
     end do
     deallocate(Diags3Dloc)
#ifdef use_PDAF
     ! local diagnostics to global field
     
     cffields(id_s_bio_dic         )%instantconc(1:nlay,n) = cffields(id_s_bio_dic)         %loc(1:nlay)          
     cffields(id_s_bio_livingmatter)%instantconc(1:nlay,n) = cffields(id_s_bio_livingmatter)%loc(1:nlay) 
     cffields(id_s_bio_deadmatter  )%instantconc(1:nlay,n) = cffields(id_s_bio_deadmatter)  %loc(1:nlay)   
     cffields(id_s_bio_alk         )%instantconc(1:nlay,n) = cffields(id_s_bio_alk)         %loc(1:nlay)          
     deallocate(cffields(id_s_bio_dic)         %loc )
     deallocate(cffields(id_s_bio_livingmatter)%loc )
     deallocate(cffields(id_s_bio_deadmatter)  %loc )
     deallocate(cffields(id_s_bio_alk)         %loc )
#endif

  end do

! ======================================================================================
!************************** EXCHANGE NODAL INFORMATION *********************************			

  do n=1, benthos_num
    call exchange_nod(Benthos(:,n))
  end do
  
  do n=1, 8
    call exchange_nod(Diags2D(:,n))	
  end do

  call exchange_nod(GloPCO2surf)	
  call exchange_nod(GloCO2flux)	
  call exchange_nod(GloCO2flux_seaicemask)

  do n=1, 4
    call exchange_nod(GlodecayBenthos(:,n))
  end do

  if (ciso) then
    call exchange_nod(GloPCO2surf_13)
    call exchange_nod(GloPCO2surf_14)
    call exchange_nod(GloCO2flux_13)
    call exchange_nod(GloCO2flux_14)
    call exchange_nod(GloCO2flux_seaicemask_13)
    call exchange_nod(GloCO2flux_seaicemask_14)  
  end if
  call exchange_nod(GloO2flux_seaicemask)	
  call exchange_nod(GloHplus)	
  call exchange_nod(AtmFeInput)	
  call exchange_nod(AtmNInput)	
  call exchange_nod(PAR3D)	

end subroutine get_recom_diags
