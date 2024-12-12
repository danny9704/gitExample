subroutine SOLVE(NB,NT,AK)
    ! 방정식을 구성한 뒤 ZSWEEP을 통해 해를 구한다.
    ! ZFI에 파동 포텐셜 계수 결과가 저장됨.
    use parameters
    use global
    implicit none
    integer, intent(in) :: NB,NT
    real(dp), intent(in) :: AK
    complex(dp) :: ZSA(MX,NP),ZSB(MX,NEQ),ZAA(NP,NP),ZBB(NP,NEQ)
    complex(dp) :: ZS(NP),ZD(NP)
    real(dp) :: SS(NP),DD(NP)
    integer :: I,J,K,M
    complex(dp) :: Z0,ZI
  
    Z0=cmplx(0.0_dp,0.0_dp,dp)
    ZI=cmplx(0.0_dp,1.0_dp,dp)
    
    ! 초기화
    do I=1,NB
      do J=1,NB
        ZAA(I,J)=Z0
      end do
      do M=1,NEQ
        ZBB(I,M)=Z0
      end do
    end do
  
    do I=1,NT
      do J=1,NB
        ZSA(I,J)=Z0
      end do
      do M=1,NEQ
        ZSB(I,M)=Z0
      end do
      if(I.le.NB) ZSA(I,I)=cmplx(PI,0.0_dp,dp)
    end do
  
    ! 각 점에 대해 SDSUB, SDCAL 호출로 ZSA, ZSB 구성
    do I=1,NT
      call SDSUB(XP(I),YP(I),NB,SS,DD)
      call SDCAL(XP(I),YP(I),AK,NB,ZS,ZD)
  
      do J=1,NB
        ZSA(I,J)=ZSA(I,J)+DD(J)+ZD(J)
      end do
  
      do M=1,3
        do J=1,NB
          ZSB(I,M)=ZSB(I,M)+(SS(J)+ZS(J))*VN(M,J)
        end do
      end do
      ZSB(I,4)=PI2*cexp(-AK*(YP(I)-ZI*XP(I)))
    end do
  
    ! ZAA, ZBB 구성
    do I=1,NB
      do J=1,NB
        do K=1,NT
          ZAA(I,J)=ZAA(I,J)+ZSA(K,I)*ZSA(K,J)
        end do
      end do
      do M=1,NEQ
        do K=1,NT
          ZBB(I,M)=ZBB(I,M)+ZSA(K,I)*ZSB(K,M)
        end do
      end do
    end do
  
    ! 복소행렬 가우스 소거
    call ZSWEEP(NP,NB,ZAA,ZBB,NEQ,SML)
    if(cdabs(ZAA(1,1)).lt.SML) write(6,600)
  
    do M=1,NEQ
      do I=1,NB
        ZFI(M,I)=ZBB(I,M)
      end do
    end do
  
  600 format(//10X,'*** ERROR: ZSWEEP IN SUBROUTINE (SOLVE) WAS ABNORMALLY DONE.',/23X,'PLEASE CHECK!'///)
  
  end subroutine SOLVE
  

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로
!3.cmplx로 실수와 허수를 받아서 복소수 변환