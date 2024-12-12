subroutine EZE1Z(XX,YY,EC,ES)
    ! E(Z) = EC + i*ES 형태의 특수 복소함수를 계산하는 서브루틴.
    ! 주어진 점 (XX, YY)에 대한 파동 문제의 특정 보조 함수 값을 제공.
    use parameters
    implicit none
  
    real(dp), intent(in) :: XX,YY   ! 입력 좌표
    complex(dp), intent(out) :: EC,ES  ! 결과로 반환되는 복소함수의 실수부(EC), 허수부(ES)
    
    ! 지역 변수들:
    real(dp) :: X,Y,R,C,ER,EI,SB,FN,CN,OLD,NEW,EXC,EXS,CC,SS
    integer :: N,J
    complex(dp) :: Z,Z1,ZSUB,ZS
    real(dp), parameter :: GAMMA=0.5772156649015d0  ! 오일러-마스케로니 상수
  
    ! 좌표 변환 및 기본값 계산
    X = XX
    Y = dabs(YY)
    R = dsqrt(X*X+Y*Y)
    C = datan2(Y,X)
  
    ! 조건에 따라 다른 계산 방식 적용
    ! 1) R이 작지 않거나 특정 조건에서 시리즈 전개 사용
    ! 2) 특정 조건에서 다른 해석적 방법(continued fraction 형태?)
    ! 3) 큰 R에서 또 다른 근사 적용
  
    ! 첫 번째 분기: R이 충분히 작거나 중간범위
    if (R.gt.25.0_dp) goto 30
    if (X.gt.0.0_dp .and. R.gt.8.0_dp) goto 20
    if (X.le.0.0_dp .and. Y.gt.10.0_dp) goto 20
  
    ! [CASE 1] 시리즈 확장법
    ER = -GAMMA - dlog(R) + R*dcos(C)
    EI = -C + R*dsin(C)
    SB = -R
    do N=2,100
      FN = real(N,dp)
      CN = C*FN
      SB = -SB*R*(FN-1.0_dp)/FN/FN
      ER = ER - SB*dcos(CN)
      EI = EI - SB*dsin(CN)
      if (N.eq.100)  goto 1
      if (EI.eq.0.0_dp)  goto 10
      if (dabs(SB/EI).le.1.0e-8_dp) goto 10
      cycle
  10: if (dabs(SB/ER).le.1.0e-8_dp) goto 1
    end do
  1: CC = dexp(X)*dcos(Y)
     SS = dexp(X)*dsin(Y)
     ! 복소함수 E(Z) = (CC+ i*0)*(ER + i*EI 변환)
     EC = cmplx(CC*ER - SS*EI, CC*EI + SS*ER, dp)
     ES = EC
     if (YY.lt.0.0_dp) ES = cmplx(dreal(EC), -dimag(EC), dp)
     return
  
  20: ! [CASE 2] 다른 해석적 방법 사용
     Z = cmplx(X,Y,dp)
     Z1 = cmplx(1.0_dp,0.0_dp,dp)
     ZSUB = cmplx(10.0_dp,0.0_dp,dp)
     ZS = Z + ZSUB/(Z1 + ZSUB/Z)
     do J=1,9
       ZSUB = cmplx(real(10-J,dp),0.0_dp,dp)
       ZS = Z + ZSUB/(Z1 + ZSUB/ZS)
     end do
     ZSUB = Z1/ZS
     EC = cmplx(dreal(ZSUB), dimag(ZSUB), dp)
     ES = EC
     if (YY.lt.0.0_dp) ES = cmplx(dreal(EC), -dimag(EC), dp)
     return
  
  30: ! [CASE 3] 큰 R에서의 근사
     OLD = -1.0_dp/R
     EXC = OLD*dcos(C)
     EXS = OLD*dsin(C)
     do N=2,100
       NEW = -OLD/R*real(N-1,dp)
       if(EXS.eq.0.0_dp) goto 31
       if(dabs(NEW/EXS).le.1.0e-8_dp) goto 31
       goto 32
  31:  if(EXC.eq.0.0_dp) goto 32
       if(dabs(NEW/EXC).le.1.0e-8_dp) goto 33
  32:  if(dabs(OLD).lt.dabs(NEW)) goto 33
       OLD=NEW
       EXC=EXC+OLD*dcos(C*real(N,dp))
       EXS=EXS+OLD*dsin(C*real(N,dp))
     end do
  33: EC=cmplx(-EXC,0.0_dp,dp)
     ES=cmplx(0.0_dp,EXS,dp)
     if(dabs(PI-dabs(C)).lt.1.0e-10_dp) ES=cmplx(dreal(ES),-PI*dexp(X),dp)
     if(YY.lt.0.0_dp) ES=cmplx(dreal(ES),-dimag(ES),dp)
     return
  
  end subroutine EZE1Z
  

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로
!3.goto 부문 압축 및 cycle함수 추가