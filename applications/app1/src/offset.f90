subroutine OFFSET(NB,NT,H0,SIG1,SIG2,OGD,KZZB,IPRINT)
    ! 선체 형상 정의 및 계산
    ! KZZB: 관성반경 비율
  
    use parameters
    use global
    implicit none
    integer, intent(in) :: NB,NT,IPRINT
    real(dp), intent(in) :: H0,SIG1,SIG2,OGD,KZZB
    integer :: J,I,II
    real(dp) :: TH,DX,DY,D,AMB,RSUB,AMD,A1,A3,SUM,OBM,DS
  
    IAD=NT-NB
    DTH=PI/real(NB,dp)
  
    SIGMA=SIG1
    RSUB=(H0+1.0_dp)**2+8.0_dp*H0*(1.0_dp-4.0_dp*SIGMA/PI)
    AMD =0.25_dp*(3.0_dp*(H0+1.0_dp)-dsqrt(RSUB))
    A1  =(0.5_dp*(H0-1.0_dp))/AMD
    A3  =(0.5_dp*(H0+1.0_dp))/AMD-1.0_dp
    AMB =AMD/H0
  
    do J=1,NB/2+1
      TH=PI05-DTH*real(J-1,dp)
      XQ(J)=AMB*((1.0_dp+A1)*dsin(TH)-A3*dsin(3.0_dp*TH))
      YQ(J)=AMB*((1.0_dp-A1)*dcos(TH)+A3*dcos(3.0_dp*TH))
    enddo
  
    SIGMA=SIG2
    RSUB=(H0+1.0_dp)**2+8.0_dp*H0*(1.0_dp-4.0_dp*SIGMA/PI)
    AMD=0.25_dp*(3.0_dp*(H0+1.0_dp)-dsqrt(RSUB))
    A1=0.5_dp*(H0-1.0_dp)/AMD
    A3=0.5_dp*(H0+1.0_dp)/AMD-1.0_dp
    AMB=AMD/H0
  
    do J=NB/2+2,NB+1
      TH = PI05-DTH*real(J-1,dp)
      XQ(J)=AMB*((1.0_dp+A1)*dsin(TH)-A3*dsin(3.0_dp*TH))
      YQ(J)=AMB*((1.0_dp-A1)*dcos(TH)+A3*dcos(3.0_dp*TH))
    end do
  
    open(newunit=IOFF, file="Offset.dat", status="replace")
    do J=1,NB+1
      write(IOFF, "(i5,99(1pe15.6))") J, XQ(J), YQ(J)
    end do
    close(IOFF)
  
    ! 패널 중심과 법선 계산
    do I=1, NB
      XP(I)=(XQ(I+1)+XQ(I))/2.0_dp
      YP(I)=(YQ(I+1)+YQ(I))/2.0_dp
      DX=XQ(I+1)-XQ(I)
      DY=YQ(I+1)-YQ(I)
      D =dsqrt(DX*DX+DY*DY)
      VN(1,I)= DY/D
      VN(2,I)=-DX/D
      VN(3,I)=XP(I)*VN(2,I)-YP(I)*VN(1,I)
    end do
  
    if(IAD.ne.0) then
      DS=(XQ(1)-XQ(NB+1))/real(IAD+1,dp)
      do I=1,IAD
        II=NB+I
        XP(II)=XQ(NB+1)+DS*real(I,dp)
        YP(II)=0.0_dp
      end do
    end if
  
    ! 질량, 복원력, GM 계산
    CMAS=(SIG1+SIG2)/H0
    C22 =(XQ(1)-XQ(NB+1))/XQ(1)
    OG  =OGD/H0
    KZZ =KZZB
    SUM=0.0_dp
    do J=1,NB
      S1 =YQ(J+1)-YQ(J)
      S2 =XQ(J  )*(2.0_dp*YQ(J  )+YQ(J+1))
      S3 =XQ(J+1)*(2.0_dp*YQ(J+1)+YQ(J  ))
      SUM=SUM+S1*(S2+S3)
    end do
    OBM=SUM/6.0_dp
    GM =(2.0_dp/3.0_dp-OBM)/CMAS+OG
  
    write(6,600) CMAS,C22,OGD,KZZ,GM
    if(IPRINT.eq.0) return
    write(6,610)
    do J=1,NB+1
      write(6,620) J,XQ(J),YQ(J),XP(J),YP(J)
    end do
  
  600 format(
       &    15X,'NONDIMENSIONAL MASS------- S/(B/2)**2=',F8.5,
       &   /15X,'HEAVE RESTORING FORCE COEFF--AW/(B/2)=',F8.5,
       &   /15X,'CENTER OF GRAVITY----------------OG/D=',F8.5,
       &   /15X,'GYRATIONAL RADIUS-----------KZZ/(B/2)=',F8.5,
       &   /15X,'METACENTRIC HEIGHT-----------GM/(B/2)=',F8.5/)
  610 format(/15X,'***** CHECK OF ORDINATES *****'
       &   /8X,'J',6X,'XQ',8X,'YQ',10X,'XP',8X,'YP')
  620 format(7X,I2,1X,2F10.5,2X,2F10.5)
  
  end subroutine OFFSET
  
  
  !주요 변경 사항 :
  !1.double precession을 real(dP)로 
  !2.d0를 정확도 확보를 위해 _dp로