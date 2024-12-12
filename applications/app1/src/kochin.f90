subroutine KOCHIN(NB,AK,IPRINT)
    use parameters
    use global
    implicit none
    integer, intent(in) :: NB,IPRINT
    real(dp), intent(in) :: AK
    real(dp) :: ABAR(3),EPS(3),RTD
    integer :: J,I,M
    real(dp) :: DX,DY,D,CDEL,SDEL
    complex(dp) :: Z0,ZI,ZETA,ZEOLD,ZENEW,ZFFA,ZFGA,ZFFB,ZFGB,ZSYM,ZANT,ZHRA,ZHRB
  
    Z0=cmplx(0.0_dp,0.0_dp,dp)
    ZI=cmplx(0.0_dp,1.0_dp,dp)
    do M=1,4
      ZHA(M)=Z0
      ZHB(M)=Z0
    end do
  
    ZETA=-AK*(YQ(1)-ZI*XQ(1))
    ZEOLD=cexp(ZETA)
    do J=1,NB
      DX=XQ(J+1)-XQ(J)
      DY=YQ(J+1)-YQ(J)
      D =dsqrt(DX*DX+DY*DY)
      CDEL=DX/D
      SDEL=DY/D
      ZFFA =-(SDEL+ZI*CDEL)/AK*(cexp(-AK*(YQ(J+1)-ZI*XQ(J+1)))-ZEOLD)
      ZFGA =-ZI*(cexp(-AK*(YQ(J+1)-ZI*XQ(J+1)))-ZEOLD)
      ZFFB =dconjg(ZFFA)
      ZFGB =dconjg(ZFGA)
      ZEOLD=cexp(-AK*(YQ(J+1)-ZI*XQ(J+1)))
  
      do M=1,3
        ZHA(M)=ZHA(M)+VN(M,J)*ZFFA - ZFI(M,J)*ZFGA
        ZHB(M)=ZHB(M)+VN(M,J)*ZFFB - ZFI(M,J)*ZFGB
      end do
      ZHA(4)=ZHA(4)-ZFI(4,J)*ZFGA
      ZHB(4)=ZHB(4)-ZFI(4,J)*ZFGB
    end do
  
    do I=1,3
      ABAR(I)=AK*cdabs(ZHA(I))
      EPS(I)=datan2(-dreal(ZHA(I)),dimag(ZHA(I)))
    end do
  
    ZSYM=cexp(ZI*EPS(2))*dcos(EPS(2))*ZI
    ZANT=cexp(ZI*EPS(1))*dsin(EPS(1))
    ZHRA=ZSYM-ZANT
    ZHRB=ZSYM+ZANT
  
    if(IPRINT.eq.0) return
    RTD=180.0_dp/PI
    write(6,600) AK,(ABAR(I),EPS(I)*RTD,I=1,3)
    write(6,610) ZHA(4),ZHRA,ZHB(4),ZHRB
  
  600 format(/6X,'******* KOCHIN FUNCTION & ACCURACY CHECK ( K*B/2=',F8.4,' ) *******',//
       & 5X,2(15X,'A-BAR',6X,'EPS(DEG)'),
       & /9X,'SWAY: ',E12.5,2X,F9.3,4X,'HEAVE: ',E12.5,2X,F9.3,
       & /9X,'ROLL: ',E12.5,2X,F9.3)
  610 format(/27X,'DIRECT CALCULATION',7X,'COMPUTED FROM RADIATION',
       & /5X,'DIFFRACTION (+) ',2E13.4,2X,2E13.4,
       & /5X,'DIFFRACTION (-) ',2E13.4,2X,2E13.4)
  
  end subroutine KOCHIN

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로