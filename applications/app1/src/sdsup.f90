subroutine SDSUB(XPI,YPI,NB,SS,DD)
    use global
    use parameters
    implicit none
    integer, intent(in) :: NB
    double precision, intent(in) :: XPI,YPI
    double precision, intent(out) :: SS(NB),DD(NB)
    integer :: J,L
    double precision :: SWA,DWA,DX,DY,D,CDEL,SDEL,XA,XB,YA,YB,SUBA,SUBB,COEF,ABSC
    double precision :: WA1,WA2,WA3
    double precision :: SL
  
    do J=1,NB
      SWA=0.0D0
      DWA=0.0D0
      if(dabs(YPI).lt.1.0D-8) goto 10
        DX=XQ(J+1)-XQ(J)
        DY=YQ(J+1)-YQ(J)
        D=dsqrt(DX*DX+DY*DY)
        CDEL=DX/D
        SDEL=DY/D
        XA=XPI-XQ(J)
        XB=XPI-XQ(J+1)
  
        SL=-1.0D0
        do L=1,2
          SL=-SL
          YA=SL*YPI-YQ(J)
          YB=SL*YPI-YQ(J+1)
          SUBA=XA*CDEL+YA*SDEL
          SUBB=XB*CDEL+YB*SDEL
          COEF=XA*SDEL - YA*CDEL
          ABSC=dabs(COEF)
          WA1=0.5D0*(SUBB*dlog(XB*XB+YB*YB)-SUBA*dlog(XA*XA+YA*YA))
          if(ABSC.lt.1.0D-10) then
             WA2=0.0D0
             WA3=0.0D0
          else
             WA2=ABSC*(datan(SUBB/ABSC)-datan(SUBA/ABSC))
             WA3=WA2/COEF
          endif
          SWA=SWA-(WA1+WA2)*SL
          DWA=DWA+WA3*SL
        end do
  10: SS(J)=SWA
      DD(J)=DWA
    enddo
  
  end subroutine SDSUB

!주요 변경 사항 :
!1.double precession을 real(dP)로 
!2.d0를 정확도 확보를 위해 _dp로