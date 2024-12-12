subroutine TRCOEF(AK,IPRINT)
    ! 방사, 회절 문제로부터 전파된 파랑에 대한 투과, 반사 계수 계산
    use parameters
    use global
    implicit none
    real(dp), intent(in) :: AK
    integer, intent(in) :: IPRINT
    real(dp) :: TT,RR,CDIF,CT,CR,CFRE
    complex(dp) :: ZI,ZTDIF,ZRDIF,ZTFRE,ZRFRE
    real(dp) :: S
    integer I
  
    ZI=cmplx(0.0_dp,1.0_dp,dp)
  
    ZTDIF=1.0_dp+ZI*ZHB(4)
    ZRDIF=ZI*ZHA(4)
    TT   =cdabs(ZTDIF)
    RR   =cdabs(ZRDIF)
    CDIF =TT**2+RR**2
  
    ZTFRE=ZTDIF
    ZRFRE=ZRDIF
    S=1.0_dp
    do I=1,3
      S=-S
      ZTFRE=ZTFRE-ZI*AK*ZMTN(I)*ZHB(I)
      ZRFRE=ZRFRE-ZI*AK*ZMTN(I)*ZHA(I)
    end do
    CT   =cdabs(ZTFRE)
    CR   =cdabs(ZRFRE)
    CFRE =CT**2+CR**2
  
    if(IPRINT.eq.0) return
    write(6,600) AK,TT,RR,CDIF,CT,CR,CFRE
  600 format(//5X,'******* TRANSMISSION & REFLECTION COEFF. ( K*B/2=',F8.4,' ) *******',/29X,'CT',12X,'CR',8X,'CT**2+CR**2',
       &   /10X,'DIFFRACTION',2X,E11.4,3X,E11.4,3X,E12.5,
       &   /10X,'MOTION FREE',2X,E11.4,3X,E11.4,3X,E12.5)
  
  end subroutine TRCOEF  

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로