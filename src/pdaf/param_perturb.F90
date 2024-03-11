MODULE mod_perturbation_pdaf
   
   ! USE
   IMPLICIT NONE
   SAVE
   
   LOGICAL :: perturb_parameters = .false.  ! whether to initialize with BGC parameter perturbation
   REAL    :: perturb_scale = 0.25       ! scaling factor for BGC parameter perturbation
   
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

END MODULE mod_perturbation_pdaf
