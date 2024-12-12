module parameters
  
  use iso_fortran_env, only: dp => real64
  implicit none

  integer, parameter :: MX = 105, NP = 100, NQ = 101, NEQ = 4   ! 기초변수
  real(dp), parameter :: SML = 1.0_dp * 10.0_dp**(-14)          ! 축옮길 때 사용하는 값
  real(dp), parameter :: PI = 3.14159265358979_dp
  real(dp), parameter :: PI05 = PI * 0.5_dp
  real(dp), parameter :: PI2 = PI * 2.0_dp

end module parameters

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로