MODULE RossLiBRDF
USE RTTYPE, ONLY : PI

IMPLICIT NONE
LOGICAL :: fRossLiBRDF=.FALSE.

REAL*8 :: fiso,fvol,fgeo
REAL*8, PARAMETER :: PARA_HoB=2.0d0,PARA_BoR=1.0d0
!REAL*8, PARAMETER :: PI=3.141592653589793d0

CONTAINS

SUBROUTINE RossLi_Kernel_vol(COST1,PHI1,COST2,PHI2,BRDML)
IMPLICIT none

REAL*8,INTENT(IN) :: COST1,PHI1,COST2,PHI2
REAL*8,DIMENSION(4,4),INTENT(OUT) :: BRDML

REAL*8 :: SINT1,SINT2,THETA1,THETA2, PHIV

REAL*8 :: thetascat,costhetascat,sinthetascat

SINT1=DSQRT(1.0D0-COST1*COST1)
SINT2=DSQRT(1.0D0-COST2*COST2)
PHIV=PHI2-PHI1+PI
!IF(PHIV<0.0D0)PHIV=PHIV+2.0D0*PI
!IF(PHIV>2.0D0*PI)PHIV=PHIV-2.0D0*PI

costhetascat=abs(COST1*COST2)+SINT1*SINT2*COS(PHIV)
IF(ABS(costhetascat)>1.0D0) STOP 'CHECK COST1 AND COST2 IN RossLi_Kernel_vol'
sinthetascat=DSQRT(1.0D0-costhetascat*costhetascat)
thetascat=ACOS(costhetascat)
BRDML=0.0D0
BRDML(1,1)=((0.5D0*PI-thetascat)*costhetascat+sinthetascat)/(abs(COST1)+abs(COST2))
BRDML(1,1)=BRDML(1,1)-0.25D0*PI
ENDSUBROUTINE RossLi_Kernel_vol

SUBROUTINE RossLi_Kernel_geo(COST1,PHI1,COST2,PHI2,BRDML)
IMPLICIT none

REAL*8,INTENT(IN) :: COST1,PHI1,COST2,PHI2
REAL*8,DIMENSION(4,4),INTENT(OUT) :: BRDML

REAL*8 :: SINT1,SINT2,THETA1,THETA2, PHIV,TANT1,TANT2

REAL*8 :: thetaP1,thetaP2,cosphaseanglep,CapD,rtmp1,rtmp2,costl,sintl,secP1,secP2,Ofunc

SINT1=DSQRT(1.0D0-COST1*COST1)
SINT2=DSQRT(1.0D0-COST2*COST2)
TANT1=SINT1/ABS(COST1)
TANT2=SINT2/ABS(COST2)

thetaP1=atan(PARA_BoR*TANT1)
thetaP2=atan(PARA_BoR*TANT2)

PHIV=PHI2-PHI1+PI
!IF(PHIV<0.0D0)PHIV=PHIV+2.0D0*PI
!IF(PHIV>2.0D0*PI)PHIV=PHIV-2.0D0*PI

cosphaseanglep=cos(thetaP1)*cos(thetaP2)+sin(thetaP1)*sin(thetaP2)*cos(PHIV)
rtmp1=tan(thetaP1)
rtmp2=tan(thetaP2)

CapD=sqrt(rtmp1**2+rtmp2**2-2.0d0*rtmp1*rtmp2*cos(PHIV))
secP1=1.0d0/cos(thetaP1)
secP2=1.0d0/cos(thetaP2)

costl=PARA_HoB*sqrt(CapD**2+(rtmp1*rtmp2*sin(PHIV))**2)/(secP1+secP2)
if(abs(costl)>1.0d0)then
	Ofunc=0.0d0
else
	sintl=sqrt(1.0d0-costl*costl)
	Ofunc=(acos(costl)-sintl*costl)*(secP1+secP2)/PI
!  write(*,*)'debugging Ross-Li geo kernel,thetaP1,thetaP2,CapD,costl='
!  write(*,'(4(f12.6,2x))')thetaP1,thetaP2,CapD,costl
endif
BRDML=0.0d0

BRDML(1,1)=Ofunc-secP1-secP2+0.5d0*(1.0d0+cosphaseanglep)*secP1*secP2

ENDSUBROUTINE RossLi_Kernel_geo

SUBROUTINE pBRDM_RossLi(COST1,PHI1,COST2,PHI2,BRDML)
IMPLICIT none
REAL*8,INTENT(IN) :: COST1,PHI1,COST2,PHI2
REAL*8,DIMENSION(4,4),INTENT(OUT) :: BRDML

REAL*8,DIMENSION(4,4) :: BRDML_vol,BRDML_geo

BRDML=0.0d0
call RossLi_Kernel_vol(COST1,PHI1,COST2,PHI2,BRDML_vol)
call RossLi_Kernel_geo(COST1,PHI1,COST2,PHI2,BRDML_geo)
BRDML(1,1)=fiso+fvol*BRDML_vol(1,1)+fgeo*BRDML_geo(1,1)
BRDML(1,1)=BRDML(1,1)*abs(COST1)

IF(BRDML(1,1)<0.0D0) STOP 'CHECK Ross-Li BRDF routine or coefficients'

END SUBROUTINE pBRDM_RossLi

END MODULE RossLiBRDF

!program main
!use RossLiBRDF
!
!implicit none
!REAL*8,DIMENSION(4,4) :: BRDML_vol,BRDML_geo
!
!integer :: inttmp0,inttmp,inttmp1
!real*8 :: thetav,thetai,PHI1,PHI2,rtmp,rtmp1,rtmp0
!real*8 :: COST1,COST2,albedo_geo,albedo_vol
!PHI1=0.0d0
!
!rtmp0=PI/2.0d0/18.0d0
!
!rtmp=PI/2.0d0/180.0d0
!rtmp1=2.0d0*PI/180.0d0
!write(*,*)'thetai/PI*180.0d0,albedo_vol,albedo_geo'
!
!do inttmp0=0,17
!
!thetai=inttmp0*rtmp0
!COST1=cos(thetai)
!
!albedo_geo=0.0d0
!albedo_vol=0.0d0
!do inttmp=0,179
!do inttmp1=0,179
!	thetav=inttmp*rtmp
!	PHI2=inttmp1*rtmp1
!	COST2=cos(thetav)
!
!	call RossLi_Kernel_vol(COST1,PHI1,COST2,PHI2,BRDML_vol)
!	call RossLi_Kernel_geo(COST1,PHI1,COST2,PHI2,BRDML_geo)
!    albedo_geo=albedo_geo+BRDML_geo(1,1)*COST2*sin(thetav)
!    albedo_vol=albedo_vol+BRDML_vol(1,1)*COST2*sin(thetav)
!enddo
!enddo
!albedo_geo=albedo_geo*rtmp*rtmp1
!albedo_vol=albedo_vol*rtmp*rtmp1
!write(*,'(3(f12.6,2x))')thetai/PI*180.0d0,albedo_vol,albedo_geo
!enddo
!
!thetai=0.0d0*PI/180.0d0
!write(*,*)'thetai=',thetai*180.0d0/PI
!write(*,*)'thetav,K_vol,K_geo'
!COST1=cos(thetai)
!
!PHI2=0.5*PI
!rtmp=PI/2.0d0/18.0d0
!
!do inttmp=17,0,-1
!    thetav=inttmp*rtmp
!    COST2=cos(thetav)
!
!    call RossLi_Kernel_vol(COST1,PHI1,COST2,PHI2,BRDML_vol)
!    call RossLi_Kernel_geo(COST1,PHI1,COST2,PHI2,BRDML_geo)
!    write(*,'(3(f12.6,2x))')-thetav/PI*180.0d0,BRDML_vol(1,1),BRDML_geo(1,1)
!enddo
!
!PHI2=1.5*PI
!do inttmp=1,17
!	thetav=inttmp*rtmp
!	COST2=cos(thetav)
!
!	call RossLi_Kernel_vol(COST1,PHI1,COST2,PHI2,BRDML_vol)
!	call RossLi_Kernel_geo(COST1,PHI1,COST2,PHI2,BRDML_geo)
!	write(*,'(3(f12.6,2x))')thetav/PI*180.0d0,BRDML_vol(1,1),BRDML_geo(1,1)
!enddo
!endprogram main
