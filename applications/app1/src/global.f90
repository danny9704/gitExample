module global
  ! 전역변수
  use iso_fortran_env, only: dp => real64
  use parameters
  implicit none

  integer :: IRAO, IDIF, IRAD          ! 입출력 단위
  real(dp) :: CMAS,C22,OG,KZZ,GM       ! 질량,복원력계수,중력중심,관성반경, 메타센터높이
  real(dp) :: XP(MX),YP(MX),XQ(NQ),YQ(NQ) ! 패널 중심( P), 패널 경계( Q)
  real(dp) :: VN(3,NP)                 ! 좌표계,패널의 면법선벡터
  complex(dp) :: ZFI(4,NP)             ! 복소포텐셜 분포
  complex(dp) :: ZHA(4),ZHB(4)         ! Kochin 함수 계산
  complex(dp) :: ZAB(3,3),ZEXF(3)      ! 부가질량, 감쇠, 파랑외력
  complex(dp) :: ZMTN(3), ZMTNO(3)     ! 중심에서의 Heave, Sway, Roll

  integer :: IPRINT,NPRINT,NB,NT       ! 출력 제어, 패널 수, 총 점개수
  real(dp) :: H0,SIG1,SIG2,OGD,KZZB    ! 선체 형상 관련 값
  integer :: IOFF, IAD, II
  real(dp) :: DTH,SIGMA,RSUB,AMD,A1,A3,AMB,TH
  real(dp) :: DS,SUM,S1,S2,S3,OBM

end module global

!주요 변경 사항 :
!1.double precession을 real(dP)로