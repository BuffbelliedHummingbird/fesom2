MODULE mod_perturbation_pdaf
   
   ! USE
   IMPLICIT NONE
   SAVE
   
   LOGICAL :: perturb_params_bio = .false.    ! whether to initialize REcoM with perturbed parameters
   LOGICAL :: perturb_params_phy = .false.    ! whether to initialize FESOM with perturbed parameters
   REAL    :: perturb_scale      = 0.25       ! scaling factor for BGC parameter perturbation
   REAL    :: perturb_scaleD     = 0.25       ! scaling factor for BGC parameter perturbation: for spread in DIC, DIN, Alk, O2
   REAL    :: perturb_scaleMix   = 0.10       ! scaling factor for mixing perturbation
   
   CONTAINS

! *************************
! *** PERTURB LOGNORMAL ***
! *************************

subroutine perturb_lognormal(value, stddev, iseed)

  implicit none

! *** Arguments ***
  real, intent(inout) :: value       ! value to be perturbed
  real, intent(in)    :: stddev      ! Standard deviation of lognormal distribution
  integer, intent(in)      :: iseed(4)    ! Seed for dlarnv

! *** Local variables ***
  real :: sigma2     ! Variance
  real :: logval     ! Logrithmic mean value of input value
  real :: rndval     ! Normal random value

  ! Generate random number
  CALL dlarnv(3, iseed, 1, rndval)

  sigma2 = log(1.0d0 + stddev*stddev)
  logval = log(value) -0.5d0 * sigma2

  value = exp(logval + sqrt(sigma2) * rndval) ! Eq. (A10) from Ciavatta et al. (2016)

end subroutine perturb_lognormal

! ******************************
! *** PERTURB PARAMETERS BIO ***
! ******************************

subroutine do_perturb_param_bio()

        USE mod_parallel_pdaf, &
            ONLY: mype_model, task_id
        USE recom_config, &
            ONLY: alfa, alfa_d, P_cm, P_cm_d, Chl2N_max, Chl2N_max_d, &
            deg_Chl, deg_Chl_d, graz_max, graz_max2, grazEff, grazEff2, &
            VDet, VDet_zoo2, Vdet_a, k_din, k_din_d, res_phy, res_phy_d, &
            rho_N, rho_C1, lossN, lossN_d, lossC, lossC_d, reminN, &
            reminC, calc_prod_ratio, res_het, res_zoo2, &
            biosynth, calc_diss_rate, calc_diss_rate2

        implicit none
        
        INTEGER :: iseed(4)          ! Seed for random number generation to perturb BGC parameters

        ! Selected parameters to disturb chlorophyll and biomass:
        ! alfa       = 0.14   ! Initial slope of P I curve small phytoplankton
        ! alfa_d     = 0.19   ! Initial slope of P I curve Diatoms
        ! P_cm        = 3.0   ! Small phytoplankton maximum rate of phtotosynthesis
        ! P_cm_d      = 3.5   ! Diatom maximum rate of phtotosynthesis
        ! Chl2N_max   = 3.78  ! Small phytoplankton maximum Chlorophyll a to nitrogen ratio
        ! Chl2N_max_d = 4.2   ! Diatom maximum Chlorophyll a to nitrogen ratio
        ! deg_Chl     = 0.1   ! degradation rate constant
        ! deg_Chl_d   = 0.1   ! degradation rate constant
        ! graz_max, graz_max2 ! maximum grazing rates
        ! grazEff, grazEff2   ! grazing effeciency of ZooPlankton
        
        ! Selected parameters to disturb dissolved tracers, most importantly DIC:
        ! VDet, VDet_zoo2, Vdet_a   ! Sinking speed detritus
        ! k_din, k_din_d            ! Nitrate uptake
        ! res_phy, res_phy_d        ! Respiration
        ! res_het, res_zoo2
        ! biosynth, (biosynthSi=0)
        ! calc_diss_rate, calc_diss_rate2 ! Dissolution during sinking
        ! rho_N, rho_C1             ! Remineralization of dissolved organic material
        !                             (DON -> DIN; DOC --> DIC)
        ! lossN, lossN_d            ! Loss terms phytoplankton (PhyN --> DON)
        ! lossC, lossC_d
        ! reminN, reminC            ! Remineralization of detritus
        ! calc_prod_ratio           ! How much of small phytoplankton are calcifiers

        ! iseed is the seed of the random number generator
        ! be aware that elements must be between 0 and 4095,
        ! and ISEED(4) must be odd.

        ! alfa
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'alfa ', alfa, ' NOT PERTURBED'
        iseed(1)=1
        iseed(2)=3
        iseed(3)=31+task_id ! different seed for each ensemble member
        iseed(4)=37
        CALL perturb_lognormal(alfa, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'alfa ', alfa, ' on task ', task_id
        
        ! alfa_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'alfa_d ', alfa_d, ' NOT PERTURBED'
        iseed(1)=2
        iseed(2)=5
        iseed(3)=29+task_id
        iseed(4)=41
        CALL perturb_lognormal(alfa_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'alfa_d ', alfa_d, ' on task ', task_id
       
        ! P_cm
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'P_cm ', P_cm, ' NOT PERTURBED'
        iseed(1)=2
        iseed(2)=4
        iseed(3)=23 + task_id
        iseed(4)=43
        CALL perturb_lognormal(P_cm, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'P_cm ', P_cm, ' on task ', task_id
        
        ! P_cm_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'P_cm_d ', P_cm_d, ' NOT PERTURBED'
        iseed(1)=7
        iseed(2)=4
        iseed(3)=19 + task_id
        iseed(4)=47
        CALL perturb_lognormal(P_cm_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'P_cm_d ', P_cm_d, ' on task ', task_id
        
        ! Chl2N_max
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'Chl2N_max ', Chl2N_max, ' NOT PERTURBED'
        iseed(1)=14
        iseed(2)=15
        iseed(3)=17 + task_id
        iseed(4)=53
        CALL perturb_lognormal(Chl2N_max, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'Chl2N_max ', Chl2N_max, ' on task ', task_id
        
        ! Chl2N_max_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'Chl2N_max_d ', Chl2N_max_d, ' NOT PERTURBED'
        iseed(1)=15
        iseed(2)=16
        iseed(3)=13 + task_id
        iseed(4)=59
        CALL perturb_lognormal(Chl2N_max_d,perturb_scale,iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'Chl2N_max_d ', Chl2N_max_d, ' on task ', task_id
        
        ! deg_Chl
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'deg_Chl ', deg_Chl, ' NOT PERTURBED'
        iseed(1)=8
        iseed(2)=6
        iseed(3)=11 + task_id
        iseed(4)=61
        CALL perturb_lognormal(deg_Chl, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'deg_Chl ', deg_Chl, ' on task ', task_id
        
        ! deg_Chl_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'deg_Chl_d ', deg_Chl_d, ' NOT PERTURBED'
        iseed(1)=104
        iseed(2)=97
        iseed(3)=7 + task_id
        iseed(4)=67
        CALL perturb_lognormal(deg_Chl_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'deg_Chl_d ', deg_Chl_d, ' on task ', task_id
        
        ! graz_max
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'graz_max ', graz_max, ' NOT PERTURBED'
        iseed(1)=33
        iseed(2)=16
        iseed(3)=3 + task_id
        iseed(4)=73
        CALL perturb_lognormal(graz_max, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'graz_max ', graz_max, ' on task ', task_id
        
        ! graz_max2
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'graz_max2 ', graz_max2, ' NOT PERTURBED'
        iseed(1)=3
        iseed(2)=6
        iseed(3)=5 + task_id
        iseed(4)=71
        CALL perturb_lognormal(graz_max2, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'graz_max2 ', graz_max2, ' on task ', task_id

        ! grazEff
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'grazEff ', grazEff, ' NOT PERTURBED'
        iseed(1)=2
        iseed(2)=1
        iseed(3)=6+3*task_id
        iseed(4)=13
        CALL perturb_lognormal(grazEff, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'grazEff ', grazEff, ' on task ', task_id
        
        ! grazEff2
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'grazEff2 ', grazEff2, ' NOT PERTURBED'
        iseed(1)=6
        iseed(2)=3
        iseed(3)=97+101*task_id
        iseed(4)=11
        CALL perturb_lognormal(grazEff2, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'grazEff2 ', grazEff2, ' on task ', task_id
        
        ! VDet
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'VDet ', VDet, ' NOT PERTURBED'
        iseed(1)=68
        iseed(2)=54
        iseed(3)=63+3*task_id
        iseed(4)=21
        CALL perturb_lognormal(VDet, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'VDet ', VDet, ' on task ', task_id
        
        ! VDet_zoo2
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'vdet_zoo2 ', vdet_zoo2, ' NOT PERTURBED'
        iseed(1)=87
        iseed(2)=28
        iseed(3)=93+5*task_id
        iseed(4)=23
        CALL perturb_lognormal(vdet_zoo2, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'vdet_zoo2 ', vdet_zoo2, ' on task ', task_id
        
        ! Vdet_a
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'vdet_a ', vdet_a, ' NOT PERTURBED'
        iseed(1)=45
        iseed(2)=65
        iseed(3)=41+7*task_id
        iseed(4)=27
        CALL perturb_lognormal(vdet_a, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'vdet_a ', vdet_a, ' on task ', task_id
        
        ! k_din
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'k_din ', k_din, ' NOT PERTURBED'
        iseed(1)=26
        iseed(2)=87
        iseed(3)=56+4*task_id
        iseed(4)=29
        CALL perturb_lognormal(k_din, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'k_din ', k_din, ' on task ', task_id
        
        ! k_din_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'k_din_d ', k_din_d, ' NOT PERTURBED'
        iseed(1)=61
        iseed(2)=15
        iseed(3)=14+3*task_id
        iseed(4)=33
        CALL perturb_lognormal(k_din_d, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'k_din_d ', k_din_d, ' on task ', task_id
        
        ! res_phy
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'res_phy ', res_phy, ' NOT PERTURBED'
        iseed(1)=67
        iseed(2)=50
        iseed(3)=31+2*task_id
        iseed(4)=35
        CALL perturb_lognormal(res_phy, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'res_phy ', res_phy, ' on task ', task_id
        
        ! res_phy_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'res_phy_d ', res_phy_d, ' NOT PERTURBED'
        iseed(1)=20
        iseed(2)=50
        iseed(3)=30+9*task_id
        iseed(4)=37
        CALL perturb_lognormal(res_phy_d, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'res_phy_d ', res_phy_d, ' on task ', task_id
        
        ! rho_N
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'rho_N ', rho_N, ' NOT PERTURBED'
        iseed(1)=12
        iseed(2)=24
        iseed(3)=36+11*task_id
        iseed(4)=45
        CALL perturb_lognormal(rho_N, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'rho_N ', rho_N, ' on task ', task_id
        
        ! rho_C1
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'rho_c1 ', rho_c1, ' NOT PERTURBED'
        iseed(1)=9
        iseed(2)=38
        iseed(3)=44+3*task_id
        iseed(4)=47
        CALL perturb_lognormal(rho_c1, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'rho_c1 ', rho_c1, ' on task ', task_id
        
        ! lossN
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'lossN ', lossn, ' NOT PERTURBED'
        iseed(1)=18
        iseed(2)=51
        iseed(3)=54+8*task_id
        iseed(4)=53
        CALL perturb_lognormal(lossn, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'lossN ', lossn, ' on task ', task_id
        
        ! lossN_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'lossN_d ', lossN_d, ' NOT PERTURBED'
        iseed(1)=48
        iseed(2)=73
        iseed(3)=35+5*task_id
        iseed(4)=55
        CALL perturb_lognormal(lossN_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'lossN_d ', lossN_d, ' on task ', task_id
        
        ! lossC
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'lossC ', lossC, ' NOT PERTURBED'
        iseed(1)=94
        iseed(2)=2
        iseed(3)=22+12*task_id
        iseed(4)=63
        CALL perturb_lognormal(lossC, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'lossC ', lossc, ' on task ', task_id
        
        ! lossC_d
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'lossC_d ', lossC_d, ' NOT PERTURBED'
        iseed(1)=41
        iseed(2)=74
        iseed(3)=31+8*task_id
        iseed(4)=67
        CALL perturb_lognormal(lossC_d, perturb_scale, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'lossC_d ', lossc_d, ' on task ', task_id
        
        ! reminN
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'reminN ', reminn, ' NOT PERTURBED'
        iseed(1)=80
        iseed(2)=25
        iseed(3)=97+101*task_id
        iseed(4)=69
        CALL perturb_lognormal(reminn, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'reminN ', reminn, ' on task ', task_id
        
        ! reminC    
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'reminC ', reminC, ' NOT PERTURBED'
        iseed(1)=13
        iseed(2)=15
        iseed(3)=66+17*task_id
        iseed(4)=71
        CALL perturb_lognormal(reminC, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'reminC ', reminc, ' on task ', task_id
        
        ! calc_prod_ratio   
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'calc_prod_ratio ', calc_prod_ratio, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(calc_prod_ratio, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'calc_prod_ratio ', calc_prod_ratio, ' on task ', task_id
        
        ! res_het
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'res_het ', res_het, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(res_het, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'res_het ', res_het, ' on task ', task_id
        
        ! res_zoo2
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'res_zoo2 ', res_zoo2, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(res_zoo2, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'res_zoo2 ', res_zoo2, ' on task ', task_id
        
        ! biosynth
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'biosynth ', biosynth, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(biosynth, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'biosynth ', biosynth, ' on task ', task_id
        
        ! calc_diss_rate
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'calc_diss_rate ', calc_diss_rate, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(calc_diss_rate, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'calc_diss_rate ', calc_diss_rate, ' on task ', task_id
        
        ! calc_diss_rate2
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'calc_diss_rate2 ', calc_diss_rate2, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(calc_diss_rate2, perturb_scaleD, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'calc_diss_rate2 ', calc_diss_rate2, ' on task ', task_id
        
end subroutine do_perturb_param_bio


! ******************************
! *** PERTURB PARAMETERS PHY ***
! ******************************

subroutine do_perturb_param_phy()

        USE mod_parallel_pdaf, &
            ONLY: mype_model, task_id
       USE o_param, &
            ONLY: K_ver

        implicit none
        
        INTEGER :: iseed(4)          ! Seed for random number generation to perturb BGC parameters
        
        ! Selected parameters to disturb mixing:
        ! K_Ver                     ! constant background diffusivity

        ! iseed is the seed of the random number generator
        ! be aware that elements must be between 0 and 4095,
        ! and ISEED(4) must be odd.
        
        ! K_ver
        IF (mype_model==0 .and. task_id==1) WRITE(*,'(a16,es14.4,a15)') 'K_ver', K_ver, ' NOT PERTURBED'
        iseed(1)=11
        iseed(2)=12
        iseed(3)=55+5*task_id
        iseed(4)=73
        CALL perturb_lognormal(K_ver, perturb_scaleMix, iseed)
        IF (mype_model==0) WRITE(*,'(a16,es14.4,a9,i3)') 'K_ver', K_ver, ' on task ', task_id
        
end subroutine do_perturb_param_phy

END MODULE mod_perturbation_pdaf
