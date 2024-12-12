subroutine FORCE(NB,AK,IPRINT)
    ! 부가질량, 감쇠계수, 파랑외력 등 힘 계수를 계산
    use parameters
    use global
    implicit none
    integer, intent(in) :: NB,IPRINT
    real(dp), intent(in) :: AK
    real(dp) :: A(3,3),B(3,3),BE(3,3),EAMP(3),EPHA(3)
    integer :: I,J,K
    real(dp) :: DX,DY,D,C1,C2,CHK
  
    ! ZAB, ZEXF 초기화
    do I=1,3
      do J=1,3
        ZAB(I,J)=cmplx(0.0_dp,0.0_dp,dp)
      end do
      ZEXF(I)=cmplx(0.0_dp,0.0_dp,dp)
    end do
  
    ! 표면압 일적을 통해 부가질량/감쇠/외력 계산
    do K=1,NB
      DX=XQ(K+1)-XQ(K)
      DY=YQ(K+1)-YQ(K)
      D =dsqrt(DX*DX+DY*DY)
      do I=1,3
        do J=1,3
          ZAB(I,J)=ZAB(I,J)-ZFI(J,K)*VN(I,K)*D
        end do
        ZEXF(I)=ZEXF(I)+ZFI(4,K)*VN(I,K)*D
      end do
    end do
  
    ! 실수부, 허수부 분리
    do I=1,3
      do J=1,3
        A(I,J)= dreal(ZAB(I,J))
        B(I,J)=-dimag(ZAB(I,J))
        BE(I,J)=0.5_dp*(ZHA(I)*dconjg(ZHA(J))+ZHB(I)*dconjg(ZHB(J)))
      end do
      EAMP(I)=cdabs(ZEXF(I))
      EPHA(I)=datan2(dimag(ZEXF(I)),dreal(ZEXF(I)))*180.0_dp/PI 
    end do
  
    write(IDIF, 650) AK, (EAMP(I), I=1,3), (EPHA(I), I=1,3)
    write(IRAD, 650) AK
    do I = 1, 3
      write(IRAD, 650) ( A(I, J), J = 1, 3)
    end do
    do I = 1, 3
      write(IRAD, 650) ( B(I, J), J = 1, 3)
    end do
  
    if(IPRINT.eq.0) return
    write(6,600) NB,AK
    do I=1,3
      C1=B(I,I)
      C2=BE(I,I)
      CHK=dabs(C1-C2)/dabs(C1+C2)*200.0_dp
      write(6,610) I,I,A(I,I),B(I,I),BE(I,I),CHK
    end do
    write(6,615)
    do I=1,3
      do J=1,3
        if(I.eq.J) cycle
        write(6,610) I,J,A(I,J),B(I,J),BE(I,J)
      end do
    end do
    write(6,630)
    do I=1,3
      write(6,640) I,ZEXF(I),ZHA(I),EAMP(I),EPHA(I)
    end do  
  600 format(//5X,'++++++++ ADDED-MASS & DAMPING COEFF. ( NB=',I3,', K*B/2=',F8.4,' )  +++++++',//10X,
       &    'I  J',8X,'ADDED-MASS',6X,'DAMPING',9X,
       &    'ENERGY',8X,'ERROR(%)')
  610 format(8X,'(',I2,',',I2,')',3X,E13.4,3(2X,E13.4))
  615 format(' ')
  630 format(//5X,'+++++ WAVE EXCITING FORCE +++++',
       &  //17X,'PRESSURE INTEGRAL',13X,'HASKIND-NEWMAN',/9X,'J',
       &   2(7X,'REAL',9X,'IMAG',4X),7X,'AMP',5X,'PHASE(DEG)')
  640 format(8X,I2,2E13.4,2X,2E13.4,3X,E11.4,2X,F9.3)
  650 format(99(1pe15.6))
  
  end subroutine FORCE

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로
!3.불필요한 함수 제거