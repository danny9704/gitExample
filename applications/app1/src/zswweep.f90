subroutine ZSWEEP(NDIM,N,ZA,ZB,NEQ,EPS)
    ! 복소 선형시스템 ZA * X = ZB를 가우스 소거로 풀이

    use parameters
    implicit none
  
    integer, intent(in) :: NDIM,N,NEQ
    real(dp), intent(in) :: EPS
    complex(dp), intent(inout) :: ZA(NDIM,NDIM),ZB(NDIM,NEQ)
  
    integer :: I,J,K,IP
    real(dp) :: P
    complex(dp) :: ZW
  
    do K=1,N
      P=0.0_dp
      ! k열에서 가장 절대값 큰 원소 찾기
      do I=K,N
        if(P.ge.cdabs(ZA(I,K))) cycle
        P=cdabs(ZA(I,K))
        IP=I
      end do
  
      if(P.le.EPS) goto 6
  
      ! 행교환
      if(IP.ne.K) then
        do J=K,N
          ZW=ZA(K,J)
          ZA(K,J)=ZA(IP,J)
          ZA(IP,J)=ZW
        end do
        do J=1,NEQ
          ZW=ZB(K,J)
          ZB(K,J)=ZB(IP,J)
          ZB(IP,J)=ZW
        end do
      end if
  
      if(K.eq.N) goto 70
  
     
      do J=K+1,N
        ZA(K,J)=ZA(K,J)/ZA(K,K)
      end do
      do J=1,NEQ
        ZB(K,J)=ZB(K,J)/ZA(K,K)
      end do
  
      do I=1,N
        if(I.eq.K) cycle
        if(K.ne.N) then
          do J=K+1,N
            ZA(I,J)=ZA(I,J)-ZA(I,K)*ZA(K,J)
          end do
        end if
        do J=1,NEQ
          ZB(I,J)=ZB(I,J)-ZA(I,K)*ZB(K,J)
        end do
      end do
  
 
      ZA(1,1)=cmplx(1.0_dp,0.0_dp,dp)
      return
  
  6:  ! 해가 불안정한 경우
      ZA(1,1)=cmplx(dabs(P),0.0_dp,dp)
      return
  
  70 continue
    end do
  
  end subroutine ZSWEEP
  

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로
!3.goto 부문 압축 및 cycle함수 추가