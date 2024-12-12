program MAIN
    use iso_fortran_env, only: dp => real64
    use parameters
    use global
    implicit none
  
    integer :: nTPlus
    real(dp) :: kStart, kEnd, dk, akb
    integer :: nK, iK
  
    IPRINT=1
    NPRINT=0
  
    write(6,*) "INPUT NB,H0,SIG1,SIG2: "
    read(5,*) NB,H0,SIG1,SIG2
    write(6,600) NB,H0,SIG1,SIG2
  
    OGD =0.05_dp  
    KZZB=0.35_dp
  
    write(*,*) "NT plus value:"
    read(5, *) nTPlus
    NT=NB + nTPlus
  
    ! 선체 형상 정의 및 기초 파라미터 설정
    call OFFSET(NB,NT,H0,SIG1,SIG2,OGD,KZZB,NPRINT)
  
    open(newunit=IRAO,file="MotionRAO.dat",status="replace")
    open(newunit=IDIF,file="DiffForce.dat",status="replace")
    open(newunit=IRAD,file="RadiationForce.dat",status="replace")
  
    write(6,*) "KB Start, KB End, nK: "
    read(5,*) kStart, kEnd, nK
  
    dk = (kEnd - kStart) / ( nK - 1.0_dp )
    akb = kStart
  
    ! 주파수 변화에 따라 절차실행
    do iK = 1, nK
      call SOLVE (NB,NT,akb)
      call KOCHIN(NB,akb,NPRINT)
      call FORCE (NB,akb,IPRINT)
      call MOTION(akb,IPRINT)
      call TRCOEF(akb,IPRINT)
      akb = akb + dk
    end do
  
  600 format(//14X,48('*')
       &   /19X,'2-D RADIATION AND DIFFRACTION PROBLEMS',
       &   /19X,'    OF A GENERAL-SHAPED 2-D BODY',
       &   /19X,'     BY INTEGRAL-EQUATION METHOD',/14X,48('*'),
       &  //15X,'NUMBER OF PANELS OVER THE WHOLE BODY---(NB)=',I4,
       &   /15X,'HALF-BEAM TO DRAFT RATIO---------H0(=B/2/D)=',F8.4,
       &   /15X,'SECTIONAL AREA RATIO FOR RIGHT-SIG1(=S/B/D)=',F8.4,
       &   /15X,'SECTIONAL AREA RATIO FOR LEFT--SIG2(=S/B/D)=',F8.4/)
  
  end program MAIN
  
!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로