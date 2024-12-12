subroutine SDCAL(XPI,YPI,AK,NB,ZS,ZD)
    ! 특정 점에 대한 복소 적분값을 계산
    use parameters
    use global
    implicit none
    integer, intent(in) :: NB
    real(dp), intent(in) :: XPI,YPI,AK
    complex(dp), intent(out) :: ZS(NB),ZD(NB)
    integer :: J
    real(dp) :: XX,YY,SGNX,XE,YE,RFL1,RFT1,RFL2,RFT2,DX,DY,D,CDEL,SDEL,SUB
    complex(dp) :: Z0,ZI,ZETA,ZEOLD,ZENEW,ZFC1,ZFC2,ZFS1,ZFS2,ZSUB
    complex(dp) :: EC,ES
  
    ! 초기화
    Z0=cmplx(0.0_dp,0.0_dp,dp)
    ZI=cmplx(0.0_dp,1.0_dp,dp)
    do J=1,NB
      ZS(J)=Z0
      ZD(J)=Z0
    end do
  
    ! 첫번째 점에 대한 초기값
    XX=XPI-XQ(1)
    YY=YPI+YQ(1)
    SGNX=sign(1.0_dp,XX)
    if(dabs(XX).lt.1.0e-10_dp) SGNX=0.0_dp
    XE=-AK*YY
    YE=-AK*dabs(XX)
    call EZE1Z(XE,YE,EC,ES)
    RFL1=0.5_dp*dlog(XX**2+YY**2)
    RFT1=datan2(YY,XX)
    ZFC1= EC - PI*cexp(cmplx(XE,YE,dp))*ZI
    ZFS1=(ES - PI*cexp(cmplx(XE,YE,dp)))*SGNX
  
    ! 패널을 따라 적분
    do J=1,NB
      XX=XPI-XQ(J+1)
      YY=YPI+YQ(J+1)
      SGNX=sign(1.0_dp,XX)
      if(dabs(XX).lt.1.0e-10_dp) SGNX=0.0_dp
      XE=-AK*YY
      YE=-AK*dabs(XX)
      call EZE1Z(XE,YE,EC,ES)
      RFL2=0.5_dp*dlog(XX**2+YY**2)
      RFT2=datan2(YY,XX)
      ZFC2= EC - PI*cexp(cmplx(XE,YE,dp))*ZI
      ZFS2=(ES - PI*cexp(cmplx(XE,YE,dp)))*SGNX
  
      DX=XQ(J+1)-XQ(J)
      DY=YQ(J+1)-YQ(J)
      D =dsqrt(DX*DX+DY*DY)
      CDEL=DX/D
      SDEL=DY/D
      SUB =SDEL*(RFL2-RFL1)+CDEL*(RFT2-RFT1)
      ZSUB=SDEL*(ZFC2-ZFC1)-CDEL*(ZFS2-ZFS1)
      ZS(J)=ZS(J)+2.0_dp/AK*(SUB+ZSUB)
      ZD(J)=ZD(J)-2.0_dp*(ZFS2-ZFS1)
      RFL1=RFL2
      RFT1=RFT2
      ZFC1=ZFC2
      ZFS1=ZFS2
    end do
  
  end subroutine SDCAL

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로