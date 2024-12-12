subroutine MOTION(AK,IPRINT)
    ! 부가질량, 감쇠, 복원력, 파랑외력을 통해 운동 방정식 계산
    use parameters
    use global
    implicit none
    real(dp), intent(in) :: AK
    integer, intent(in) :: IPRINT
    complex(dp) :: ZAA(3,3),ZBB(3)
    real(dp) :: AMPG(3),PHAG(3)
    complex(dp) :: ZMTNG(3)
    integer I
  
    ! 운동 방정식 계수행렬 구성
    ZAA(1,1)=-AK*(CMAS+ZAB(1,1))
    ZAA(1,2)=-AK* ZAB(1,2)
    ZAA(1,3)=-AK*(ZAB(1,3)+OG*ZAB(1,1))
    ZBB(1  )= ZEXF(1)
  
    ZAA(2,1)=-AK*ZAB(2,1)
    ZAA(2,2)=-AK*(CMAS+ZAB(2,2))+C22
    ZAA(2,3)=-AK*(ZAB(2,3)+OG*ZAB(2,1))
    ZBB(2  )= ZEXF(2)
  
    ZAA(3,1)=-AK*(ZAB(3,1)+OG*ZAB(1,1))
    ZAA(3,2)=-AK*(ZAB(3,2)+OG*ZAB(1,2))
    ZAA(3,3)=-AK*(CMAS*KZZ**2+ZAB(3,3)+OG*ZAB(1,3)
       &         +OG*(ZAB(3,1)+OG*ZAB(1,1)))+CMAS*GM
    ZBB(3  )= ZEXF(3)+OG*ZEXF(1)
  
    call ZSWEEP(3,3,ZAA,ZBB,1,SML)
    if(cdabs(ZAA(1,1)).lt.SML) write(6,600)
    do I=1,3
      ZMTNG(I)=ZBB(I)
    end do
    ZMTNO(1)=ZMTNG(1)+OG*ZMTNG(3)
    ZMTNO(2)=ZMTNG(2)
    ZMTNO(3)=ZMTNG(3)
  
    do I=1,3
      AMPG(I)=cdabs(ZMTNG(I))
      if(I.eq.3) AMPG(I)=AMPG(I)/AK
      PHAG(I)=datan2(dimag(ZMTNG(I)),dreal(ZMTNG(I)))*180.0_dp/PI
    end do
  
    write(IRAO, 620) AK, ( AMPG(I), I=1, 3), ( PHAG(I), I=1, 3)
  
    if(IPRINT.eq.0) return
    write(6,610) AK,(AMPG(I),PHAG(I),I=1,3)
  600 format(///10X,'+++ ERROR: ZSWEEP IN (MOTION) +++'///)
  610 format(//5X,'+++++ MOTIONS ABOUT ''G'' FOR K*B/2=',F7.3,
       &   '+++++',/20X,'AMP.',7X,'PHASE',/9X,'SWAY ',E11.4,
       &   2X,F9.3,' (DEG)' ,/9X, 'HEAVE ',E11.4,2X,F9.3,' (OEG)',
       &   /9X, 'ROLL ',E11.4,2X,F9.3,' (DEG)')
  620 format(99(1pe15.6))
  
  end subroutine MOTION
  
!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로
!3.불필요한 함수 제거