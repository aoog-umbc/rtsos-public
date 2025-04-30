
      SUBROUTINE gauleg(x1,x2,x,w,n)
      INTEGER n
      DOUBLE PRECISION x1,x2,x(n),w(n)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=3.d-14)

      INTEGER i,j,m
      DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
      DOUBLE PRECISION fi,fj,fn
      fn=dfloat(n)
	  
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
	    fi=dfloat(i)
        z=cos(3.14159265358979323846d0*(fi-.25d0)/(fn+.5d0))
1       continue
          p1=1.d0
          p2=0.d0
          do 11 j=1,n
		    fj=dfloat(j)
            p3=p2
            p2=p1
            p1=((2.d0*fj-1.d0)*z*p2-(fj-1.d0)*p3)/fj
11        continue
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        w(n+1-i)=w(i)
12    continue
      return
      END


! Pengwang's new math libraries
FUNCTION bisearch(nxx,xx,x)
IMPLICIT NONE
INTEGER, INTENT(IN) :: nxx
DOUBLE PRECISION, DIMENSION(nxx), INTENT(IN) :: xx
DOUBLE PRECISION, INTENT(IN) :: x
INTEGER :: bisearch
INTEGER :: jlo,jmid,jub
LOGICAL :: mono_increase

mono_increase = (xx(nxx) > xx(1))
jlo=0
jub=nxx+1
do while (jub-jlo > 1)
	jmid=(jub+jlo)/2
	if (mono_increase .eqv. (x >= xx(jmid))) then
		jlo=jmid
	else
		jub=jmid
	end if
end do
if (x == xx(1)) then
	bisearch=1
else if (x == xx(nxx)) then
	bisearch=nxx-1
else
	bisearch=jlo
end if
END FUNCTION bisearch

DOUBLE PRECISION FUNCTION Func_UVIP3P(NXX,XX,YY,X,ORDER)
IMPLICIT NONE
INTEGER,intent(in) :: NXX,ORDER
DOUBLE PRECISION,intent(in) :: XX(NXX),YY(NXX)
DOUBLE PRECISION,intent(in) :: X

DOUBLE PRECISION,DIMENSION(1) :: XI,YI
INTEGER :: IU_LO,M_PL,K_LO
integer,external :: bisearch
DOUBLE PRECISION :: RTMP,RTMP1
LOGICAL :: mono_increase,mono_decrease

mono_increase=XX(2)>XX(1)
mono_decrease=.not. mono_increase
IF((mono_increase .AND. (X<=XX(1).OR.X>=XX(NXX))) .OR. &
   (mono_decrease .AND. (X>=XX(1).OR.X<=XX(NXX))))       THEN

if( (mono_increase .and. X<=XX(1)) .or. &
    (mono_decrease .and. X>=XX(1)) ) then
   RTMP=XX(2)-XX(1)
   RTMP1=X-XX(1)
   IF(ABS(RTMP1)<=abs(RTMP))THEN
		Func_UVIP3P=YY(1)+RTMP1*(YY(2)-YY(1))/RTMP
        RETURN
   ELSE
		WRITE(*,*) 'Func_UVIP3P, X=',X,'OUT OF XX RANGE,XX(1),XX(2)=',XX(1),XX(2)
        STOP
   ENDIF
else
	RTMP=XX(NXX-1)-XX(NXX)
	RTMP1=X-XX(NXX)
    IF(ABS(RTMP1)<=abs(RTMP))THEN
		Func_UVIP3P=YY(NXX)+RTMP1*(YY(NXX-1)-YY(NXX))/RTMP
		RETURN
   ELSE
	 WRITE(*,*) 'Func_UVIP3P, X=',X,'OUT OF XX RANGE,XX(NXX-1),XX(NXX)=',XX(NXX-1),XX(NXX)
	 STOP
   ENDIF
endif

RETURN ! FINISH HANDLING OF LINEAR EXTRAPOLATION

ENDIF


M_PL=ORDER+1
XI=X

IU_LO=bisearch(NXX,XX,X)
K_LO=min(max(IU_LO-(M_PL-1)/2,1),NXX+1-M_PL)

CALL UVIP3P(ORDER,M_PL,XX(K_LO),YY(K_LO),1,XI,YI)
Func_UVIP3P=YI(1)

RETURN

END FUNCTION Func_UVIP3P

SUBROUTINE interp2(x1a,x2a,ya,m,n,x1,x2,y,dy)
INTEGER,intent(in):: m,n
DOUBLE PRECISION,intent(in):: x1,x2,x1a(m),x2a(n),ya(m,n)
DOUBLE PRECISION,intent(out):: dy,y
DOUBLE PRECISION,external :: Func_UVIP3P
DOUBLE PRECISION, dimension(:), allocatable:: ymtmp,yntmp

INTEGER j,k,ORDER
ORDER=2
allocate(ymtmp(m),yntmp(n))
do j=1,m
  do k=1,n
	yntmp(k)=ya(j,k)
  enddo
  ymtmp(j)=Func_UVIP3P(n,x2a,yntmp,x2,ORDER)
enddo
y=Func_UVIP3P(m,x1a,ymtmp,x1,ORDER)
deallocate(ymtmp,yntmp)
return
END SUBROUTINE interp2


DOUBLE PRECISION FUNCTION TRAPZOID(NXX,XX,YY)
IMPLICIT NONE
! CALCULATE INTEGRATION(YY,d_XX)
INTEGER,INTENT(IN) :: NXX
DOUBLE PRECISION,INTENT(IN),DIMENSION(1:NXX) :: XX,YY
INTEGER :: IXX

TRAPZOID=0.0D0
DO IXX=1,NXX-1
  TRAPZOID=TRAPZOID+0.5D0*(YY(IXX+1)+YY(IXX))*(XX(IXX+1)-XX(IXX))
ENDDO

END FUNCTION TRAPZOID

SUBROUTINE  UVIP3P(NP,ND,XD_IN,YD_IN,NI,XI, YI)

! Univariate Interpolation (Improved Akima Method)

! Hiroshi Akima
! U.S. Department of Commerce, NTIA/ITS
! Version of 89/07/04

! This subroutine performs univariate interpolation.  It is based
! on the improved A method developed by Hiroshi Akima, 'A method
! of univariate interpolation that has the accuracy of a third-
! degree polynomial,' ACM TOMS, vol. xx, pp. xxx-xxx, 19xx.  (The
! equation numbers referred to in the comments below are those in
! the paper.)

! In this method, the interpolating function is a piecewise
! function composed of a set of polynomials applicable to
! successive intervals of the given data points.  This method
! uses third-degree polynomials as the default, but the user has
! an option to use higher-degree polynomial to reduce undulations
! in resulting curves.

! This method has the accuracy of a third-degree polynomial if
! the degree of the polynomials for the interpolating function is
! set to three.

! The input arguments are
!   NP = degree of the polynomials for the interpolating
!        function,
!   ND = number of input data points
!        (must be equal to 2 or greater),
!   XD = array of dimension ND, containing the abscissas of
!        the input data points
!        (must be in a monotoni! increasing order),
!   YD = array of dimension ND, containing the ordinates of
!        the input data points,
!   NI = number of points for which interpolation is desired
!        (must be equal to 1 or greater),
!   XI = array of dimension NI, containing the abscissas of
!        the desired points.

! The output argument is
!   YI = array of dimension NI, where the ordinates of the
!        desired points are to be stored.

! If an integer value smaller than 3 is given to the NP argument,
! this subroutine assumes NP = 3.

! The XI array elements need not be monotonic, but this
! subroutine interpolates faster if the XI array elements are
! given in a monotoni! order.

! If the XI array element is less than XD(1) or greater than
! XD(ND), this subroutine linearly interpolates the YI value.


! Specification statement
IMPLICIT NONE
INTEGER,INTENT(IN) :: NP,ND,NI
DOUBLE PRECISION,INTENT(IN),DIMENSION(ND)::   XD_IN,YD_IN
DOUBLE PRECISION,INTENT(IN),DIMENSION(NI):: XI
DOUBLE PRECISION,INTENT(OUT),DIMENSION(NI):: YI

DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::   XD,YD

INTEGER :: ID,ID0,ID1,ID2,ID3,IDMD,IDMN,IDMX,II,IINT,IPE,NP0,NPM1,IEPT,IINTPV
DOUBLE PRECISION :: A0,A1,A12,A13,A2,A3,AA0,AA1,B0,B1,DLT,DNM,DX,DY,DY0,DY1,&
					DY2,DY3,EPSLN,PE,RENNM2,RENPM1,SMPEF,SMPEI,SMWTF,SMWTI, &
					SX,SXX,SXY,SY,T0,T1,U,UC,V,VOL,WT,X0,X1,X2,X3,XII,XX,Y0,&
					Y1,Y2,Y3,YP,YP0,YP1
! Error check
   10 IF (ND.LE.1)   GO TO 90
	  IF (NI.LE.0)   GO TO 91

   ALLOCATE(XD(ND),YD(ND))
   IF (XD_IN(2).LT.XD_IN(1)) THEN
		DO ID=1,ND
           XD(ID)=XD_IN(ND-ID+1)
           YD(ID)=YD_IN(ND-ID+1)
		ENDDO
   ELSE
		XD(1:ND)=XD_IN(1:ND)
		YD(1:ND)=YD_IN(1:ND)
   ENDIF


	  DO 11  ID=2,ND
		IF (XD(ID).LE.XD(ID-1)) GO TO 92
   11 CONTINUE
! Branches off special cases.
	  IF (ND.LE.4)   GO TO 50
! General case  --  Five data points of more
! Calculates some local variables.
   20 NP0=MAX(3,NP)
	  NPM1=NP0-1
	  RENPM1=NPM1
	  RENNM2=NP0*(NP0-2)
! Main calculation for the general case
! First (outermost) DO-loop with respect to the desired points
   30 DO 39  II=1,NI
		IF (II.EQ.1)      IINTPV=-1
		XII=XI(II)
! Locates the interval that includes the desired point by binary
! search.
		IF (XII.LE.XD(1))  THEN
		  IINT=0
		ELSE IF (XII.LT.XD(ND))  THEN
		  IDMN=1
		  IDMX=ND
		  IDMD=(IDMN+IDMX)/2
   31     IF (XII.GE.XD(IDMD))  THEN
			IDMN=IDMD
		  ELSE
			IDMX=IDMD
		  END IF
		  IDMD=(IDMN+IDMX)/2
		  IF (IDMD.GT.IDMN)    GO TO 31
		  IINT=IDMD
		ELSE
		  IINT=ND
		END IF
! End of locating the interval of interest
! Interpolation or extrapolation in one of the three subcases
		IF (IINT.LE.0)  THEN
! Subcase 1  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the first data
!                point or less.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
		  IF (IINT.NE.IINTPV)  THEN
			IINTPV=IINT
			X0=XD(1)
			X1=XD(2)-X0
			X2=XD(3)-X0
			X3=XD(4)-X0
			Y0=YD(1)
			Y1=YD(2)-Y0
			Y2=YD(3)-Y0
			Y3=YD(4)-Y0
			DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
			A1=(((X2*X3)**2)*(X3-X2)*Y1 &
			  +((X3*X1)**2)*(X1-X3)*Y2 &
			  +((X1*X2)**2)*(X2-X1)*Y3)/DLT
		  END IF
! Evaluates the YI value.
		  YI(II)=Y0+A1*(XII-X0)
! End of Subcase 1
		ELSE IF (IINT.GE.ND)  THEN
! Subcase 2  --  Linear extrapolation when the abscissa of the
!                desired point is equal to that of the last data
!                point or greater.
! Estimates the first derivative when the interval is not the
! same as the one for the previous desired point.  --
! cf. Equation (8)
		  IF (IINT.NE.IINTPV)  THEN
			IINTPV=IINT
			X0=XD(ND)
			X1=XD(ND-1)-X0
			X2=XD(ND-2)-X0
			X3=XD(ND-3)-X0
			Y0=YD(ND)
			Y1=YD(ND-1)-Y0
			Y2=YD(ND-2)-Y0
			Y3=YD(ND-3)-Y0
			DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
			A1=(((X2*X3)**2)*(X3-X2)*Y1 &
			  +((X3*X1)**2)*(X1-X3)*Y2 &
			  +((X1*X2)**2)*(X2-X1)*Y3)/DLT
		  END IF
! Evaluates the YI value.
		  YI(II)=Y0+A1*(XII-X0)
! End of Subcase 2
		ELSE
! Subcase 3  --  Interpolation when the abscissa of the desired
!                point is  between those of the first and last
!                data points.
! Calculates the coefficients of the third-degree polynomial (for
! NP.LE.3) or the factors for the higher-degree polynomials (for
! NP.GT.3), when the interval is not the same as the one for the
! previous desired point.
		  IF (IINT.NE.IINTPV)  THEN
			IINTPV=IINT
! The second DO-loop with respect to the two endpoints of the
! interval
			DO 37  IEPT=1,2
! Calculates the estimate of the first derivative at an endpoint.
! Initial setting for calculation
			  ID0=IINT+IEPT-1
			  X0=XD(ID0)
			  Y0=YD(ID0)
			  SMPEF=0.0
			  SMWTF=0.0
			  SMPEI=0.0
			  SMWTI=0.0
! The third (innermost) DO-loop with respect to the four primary
! estimate of the first derivative
			  DO 36  IPE=1,4
! Selects point numbers of four consecutive data points for
! calculating the primary estimate of the first derivative.
				IF (IPE.EQ.1)  THEN
				  ID1=ID0-3
				  ID2=ID0-2
				  ID3=ID0-1
				ELSE IF (IPE.EQ.2)  THEN
				  ID1=ID0+1
				ELSE IF (IPE.EQ.3)  THEN
				  ID2=ID0+2
				ELSE
				  ID3=ID0+3
				END IF
! Checks if any point number falls outside the legitimate range
! (between 1 and ND).  Skips calculation of the primary estimate
! if any does.
				IF (ID1.LT.1.OR.ID2.LT.1.OR.ID3.LT.1.OR. &
				   ID1.GT.ND.OR.ID2.GT.ND.OR.ID3.GT.ND) &
					GO TO 36
! Calculates the primary estimate of the first derivative  --
! cf. Equation (8)
				X1=XD(ID1)-X0
				X2=XD(ID2)-X0
				X3=XD(ID3)-X0
				Y1=YD(ID1)-Y0
				Y2=YD(ID2)-Y0
				Y3=YD(ID3)-Y0
				DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
				PE=(((X2*X3)**2)*(X3-X2)*Y1 &
				  +((X3*X1)**2)*(X1-X3)*Y2 &
				  +((X1*X2)**2)*(X2-X1)*Y3)/DLT
! Calculates the volatility factor, VOL, and distance factor,
! SXX, for the primary estimate.  --  cf. Equations (9) and (11)
				SX=X1+X2+X3
				SY=Y1+Y2+Y3
				SXX=X1*X1+X2*X2+X3*X3
				SXY=X1*Y1+X2*Y2+X3*Y3
				DNM=4.0*SXX-SX*SX
				B0=(SXX*SY-SX*SXY)/DNM
				B1=(4.0*SXY-SX*SY)/DNM
				DY0=-B0
				DY1=Y1-(B0+B1*X1)
				DY2=Y2-(B0+B1*X2)
				DY3=Y3-(B0+B1*X3)
				VOL=DY0*DY0+DY1*DY1+DY2*DY2+DY3*DY3
! Calculates the EPSLN value, which is used to decide whether or
! not the volatility factor, VOL, is essentially zero.
				EPSLN=(YD(ID0)**2+YD(ID1)**2  &
					 +YD(ID2)**2+YD(ID3)**2)*1.0E-12
! Accumulates the weighted primary estimates.  --
! cf. Equations (13) and (14)
				IF (VOL.GT.EPSLN)  THEN
! - For finite weight.
				  WT=1.0/(VOL*SXX)
				  SMPEF=SMPEF+PE*WT
				  SMWTF=SMWTF+WT
				ELSE
! - For infinite weight.
				  SMPEI=SMPEI+PE
				  SMWTI=SMWTI+1.0
				END IF
   36         CONTINUE
! End of the third DO-loop
! Calculates the final estimate of the first derivative.  --
! cf. Equation (14)
			  IF (SMWTI.LT.0.5)  THEN
! - When no infinite weights exist.
				YP=SMPEF/SMWTF
			  ELSE
! - When infinite weights exist.
				YP=SMPEI/SMWTI
			  END IF
			  IF (IEPT.EQ.1)  THEN
				YP0=YP
			  ELSE
				YP1=YP
			  END IF
! End of the calculation of the estimate of the first derivative
! at an endpoint
   37       CONTINUE
! End of the second DO-loop
			IF (NP0.LE.3)  THEN
! Calculates the coefficients of the third-degree polynomial
! (when NP.LE.3).  --  cf. Equation (4)
			  DX=XD(IINT+1)-XD(IINT)
			  DY=YD(IINT+1)-YD(IINT)
			  A0=YD(IINT)
			  A1=YP0
			  YP1=YP1-YP0
			  YP0=YP0-DY/DX
			  A2=-(3.0*YP0+YP1)/DX
			  A3= (2.0*YP0+YP1)/(DX*DX)
			ELSE
! Calculates the factors for the higher-degree polynomials
! (when NP.GT.3).  --  cf. Equation (20)
			  DX=XD(IINT+1)-XD(IINT)
			  DY=YD(IINT+1)-YD(IINT)
			  T0=YP0*DX-DY
			  T1=YP1*DX-DY
			  AA0= (T0+RENPM1*T1)/RENNM2
			  AA1=-(RENPM1*T0+T1)/RENNM2
			END IF
		  END IF
! End of the calculation of the coefficients of the third-degree
! polynomial (when NP.LE.3) or the factors for the higher-degree
! polynomials (when NP.GT.3), when the interval is not the same
! as the one for the previous desired point.
! Evaluates the YI value.
		  IF (NP0.LE.3)  THEN
! - With a third-degree polynomial.  --  cf. Equation (3)
			XX=XII-XD(IINT)
			YI(II)=A0+XX*(A1+XX*(A2+XX*A3))
		  ELSE
! - With a higher-degree polynomial.  --  cf. Equation (19)
			U=(XII-XD(IINT))/DX
			UC=1.0-U
			V=AA0*((U**NP0)-U)+AA1*((UC**NP0)-UC)
			YI(II)=YD(IINT)+DY*U+V
		  END IF
! End of Subcase 3
		END IF
   39 CONTINUE
! End of the first DO-loop
! End of general case
	  DEALLOCATE(XD,YD)
	  RETURN
! Special cases  --  Four data points or less
! Preliminary processing for the special cases
   50 X0=XD(1)
	  Y0=YD(1)
	  X1=XD(2)-X0
	  Y1=YD(2)-Y0
	  IF (ND.EQ.2)   GO TO 60
	  X2=XD(3)-X0
	  Y2=YD(3)-Y0
	  IF (ND.EQ.3)   GO TO 70
	  X3=XD(4)-X0
	  Y3=YD(4)-Y0
	  GO TO 80
! Special Case 1  --  Two data points
! (Linear interpolation and extrapolation)
   60 A1=Y1/X1
	  DO 61  II=1,NI
		YI(II)=Y0+A1*(XI(II)-X0)
   61 CONTINUE
! End of Special Case 1
		DEALLOCATE(XD,YD)
	  RETURN
! Special Case 2  --  Three data points
! (Quadrati! interpolation and linear extrapolation)
   70 DLT=X1*X2*(X2-X1)
	  A1=(X2*X2*Y1-X1*X1*Y2)/DLT
	  A2=(X1*Y2-X2*Y1)/DLT
	  A12=2.0*A2*X2+A1
	  DO 71  II=1,NI
		XX=XI(II)-X0
		IF (XX.LE.0.0)  THEN
		  YI(II)=Y0+A1*XX
		ELSE IF (XX.LT.X2) THEN
		  YI(II)=Y0+XX*(A1+XX*A2)
		ELSE
		  YI(II)=Y0+Y2+A12*(XX-X2)
		END IF
   71 CONTINUE
! End of Special Case 2
	  DEALLOCATE(XD,YD)
	  RETURN
! Special Case 3  --  Four data points
! (Cubi! interpolation and linear extrapolation)
   80 DLT=X1*X2*X3*(X2-X1)*(X3-X2)*(X3-X1)
	  A1=(((X2*X3)**2)*(X3-X2)*Y1 &
		+((X3*X1)**2)*(X1-X3)*Y2 &
		+((X1*X2)**2)*(X2-X1)*Y3)/DLT
	  A2=(X2*X3*(X2*X2-X3*X3)*Y1 &
		+X3*X1*(X3*X3-X1*X1)*Y2 &
		+X1*X2*(X1*X1-X2*X2)*Y3)/DLT
	  A3=(X2*X3*(X3-X2)*Y1 &
		+X3*X1*(X1-X3)*Y2 &
		+X1*X2*(X2-X1)*Y3)/DLT
	  A13=(3.0*A3*X3+2.0*A2)*X3+A1
	  DO 81  II=1,NI
		XX=XI(II)-X0
		IF (XX.LE.0.0)  THEN
		  YI(II)=Y0+A1*XX
		ELSE IF (XX.LT.X3) THEN
		  YI(II)=Y0+XX*(A1+XX*(A2+XX*A3))
		ELSE
		  YI(II)=Y0+Y3+A13*(XX-X3)
		END IF
   81 CONTINUE
! End of Special Case 3
	  DEALLOCATE(XD,YD)
	  RETURN
! Error exit
   90 WRITE (*,99090) ND
	  GO TO 99
   91 WRITE (*,99091) NI
	  GO TO 99
   92 WRITE (*,99092)ND,XD(1),XD(2),ID,XD(ID-1),XD(ID)
   99 WRITE (*,99099)
	  DEALLOCATE(XD,YD)
	  RETURN
! Format statements for error messages
99090 FORMAT (1X/ ' ***   Insufficient data points.',7X,'ND =',I3)
99091 FORMAT (1X/ ' ***   No desired points.',7X,'NI =',I3)
99092 FORMAT (1X/ ' ***   Two data points identical or out of ','sequence.'/ &
	   7X,'ND,XD(1),XD(2),ID, XD(ID-1), XD(ID) =',I5,2F10.3,I5,2F10.3)
99099 FORMAT (' Error detected in the UVIP3P subroutine'/)
	  ENDSUBROUTINE  UVIP3P



SUBROUTINE bisection(func, x1, x2, root, tol)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INout) :: x1, x2, tol
  DOUBLE PRECISION, INTENT(OUT) :: root
  INTEGER, PARAMETER :: max_iter=100
  INTEGER :: iter
  DOUBLE PRECISION :: fx1, fx2, froot, xmid

  INTERFACE
	  DOUBLE PRECISION FUNCTION func(x)
		  DOUBLE PRECISION, INTENT(IN) :: x
      END FUNCTION func
  END INTERFACE

  fx1 = func(x1)
  fx2 = func(x2)

  IF (fx1 * fx2 >= 0.0) THEN
	STOP "Root not bracketed"
	RETURN
  END IF

  DO iter = 1, max_iter
	xmid = (x1 + x2) / 2.0
	froot = func(xmid)

	IF (ABS(froot) < tol) THEN
	  root = xmid
	  RETURN
	END IF

	IF (fx1 * froot < 0.0) THEN
	  x2 = xmid
	  fx2 = froot
	ELSE
	  x1 = xmid
	  fx1 = froot
	END IF
  END DO

  PRINT*, "Max iterations reached"
END SUBROUTINE bisection

SUBROUTINE Err_MSG(string)
CHARACTER(LEN=*), INTENT(IN) :: string
write (*,*) 'Error Message: ',string
STOP
END SUBROUTINE Err_MSG

!SUBROUTINE swap(a,b)
!Double Precision, INTENT(INOUT) :: a,b
!Double Precision :: dum
!dum=a
!a=b
!b=dum
!END SUBROUTINE swap

!	subroutine myspline(n,x,y,nn,xn,yn)
!
!!c---------spline fit to derive the yy value at point xx
!
!!c  Inputs:
!!c     n:    the length of x and y
!!c     x(n): the x values which x(1) < x(2) ... < x(n)
!!c     y(n): the y value which correspondent to x(n)
!!c     nn:  the length of vector xx and yy
!!c     xn:  the x value at which y value is wanted
!!c
!
!!c  Outputs:
!
!!c     yn: the wanted y value from the fitting
!!c  Internal variables:
!
!!c     yp1: the derivative of y over x at x(1), for natural bc, yp1=1.e31
!!c     ypn: the derivative of y over x at x(n), for natural bc, ypn=1.e31
!!c     y2(n): the second derivatives
!
!      integer n,nn,ny2,i
!      parameter (ny2=5000)
!      DOUBLE PRECISION x(*),y(*),xn(*),yn(*),y2(ny2),xx,yy,yp1,ypn
!
!!c--------the sorting which makes sure x(1)<x(2)<...<x(n)-------
!
!!        call sort2(n,x,y)
!
!!c--------start spline------------
!
!        yp1=1.e31
!        ypn=1.e31
!
!!    yp1=(y(2)-y(1))/(x(2)-x(1))
!!	ypn=(y(n)-y(n-1))/(x(n)-x(n-1))
!
!	call spline2(x,y,n,yp1,ypn,y2)
!	do i=1,nn
!        xx = xn(i)
!        call splint2(x,y,y2,n,xx,yy)
!        yn(i) = yy
!	enddo
!
!	return
!
!   	end

	



!      SUBROUTINE spline2(x,y,n,yp1,ypn,y2)
!
!      INTEGER n,NMAX, ny2
!      PARAMETER (NMAX=5000)
!      parameter (ny2=5000)
!      DOUBLE PRECISION yp1,ypn,x(n),y(n),y2(ny2)
!      INTEGER i,k
!      DOUBLE PRECISION p,qn,sig,un,u(NMAX)
!
!      if (yp1.gt..99e30) then
!        y2(1)=0.0d0
!        u(1)=0.0d0
!      else
!        y2(1)=-0.50d0
!        u(1)=(3.0d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
!      endif
!
!      do i=2,n-1
!
!        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
!        p=sig*y2(i-1)+2.0d0
!        y2(i)=(sig-1.0d0)/p
!
!        u(i)=(6.0d0*((y(i+1)-y(i))/(x(i+1)      &
!    -x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*    &
!     u(i-1))/p
!
!      enddo
!
!      if (ypn.gt..99e30) then
!        qn=0.0d0
!        un=0.0d0
!      else
!        qn=0.50d0
!        un=(3.0d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
!      endif
!
!      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
!
!      do k=n-1,1,-1
!        y2(k)=y2(k)*y2(k+1)+u(k)
!      enddo
!
!      return
!
!      END
!
!
!
!      SUBROUTINE splint2(xa,ya,y2a,n,x,y)
!
!      INTEGER n
!      DOUBLE PRECISION x,y,xa(n),y2a(n),ya(n)
!      INTEGER k,khi,klo
!      DOUBLE PRECISION a,b,h
!
!      klo=1
!      khi=n
!
!100     if (khi-klo.gt.1) then
!
!        k=(khi+klo)/2
!
!        if(xa(k).gt.x)then
!          khi=k
!        else
!          klo=k
!        endif
!
!      goto 100
!
!      endif
!
!      h=xa(khi)-xa(klo)
!      if (h.eq.0.0d0) stop 'bad xa input in splint'
!
!      a=(xa(khi)-x)/h
!      b=(x-xa(klo))/h
!
!       y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+ &
!        (b**3-b)*y2a(khi))*(h*h)/6.0d0
!
!      return
!
!      END


!      SUBROUTINE hunt(xx,n,x,jlo)
!      INTEGER jlo,n
!      DOUBLE PRECISION x,xx(n)
!      INTEGER inc,jhi,jm
!      LOGICAL ascnd
!      ascnd=xx(n).ge.xx(1)
!      if(jlo.le.0.or.jlo.gt.n)then
!        jlo=0
!        jhi=n+1
!        goto 3
!      endif
!      inc=1
!      if(x.ge.xx(jlo).eqv.ascnd)then
!1       jhi=jlo+inc
!        if(jhi.gt.n)then
!          jhi=n+1
!        else if(x.ge.xx(jhi).eqv.ascnd)then
!          jlo=jhi
!          inc=inc+inc
!          goto 1
!        endif
!      else
!        jhi=jlo
!2       jlo=jhi-inc
!        if(jlo.lt.1)then
!          jlo=0
!        else if(x.lt.xx(jlo).eqv.ascnd)then
!          jhi=jlo
!          inc=inc+inc
!          goto 2
!        endif
!      endif
!3     if(jhi-jlo.eq.1)then
!        if(x.eq.xx(n))jlo=n-1
!        if(x.eq.xx(1))jlo=1
!        return
!      endif
!      jm=(jhi+jlo)/2
!      if(x.ge.xx(jm).eqv.ascnd)then
!        jlo=jm
!      else
!        jhi=jm
!      endif
!      goto 3
!      END

!      SUBROUTINE polint(xa,ya,n,x,y,dy)
!      INTEGER n,NMAX
!      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
!      PARAMETER (NMAX=10)
!      INTEGER i,m,ns
!      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
!      ns=1
!      dif=abs(x-xa(1))
!      do i=1,n
!        dift=abs(x-xa(i))
!        if (dift.lt.dif) then
!          ns=i
!          dif=dift
!        endif
!        c(i)=ya(i)
!        d(i)=ya(i)
!      enddo
!      y=ya(ns)
!      ns=ns-1
!      do m=1,n-1
!        do i=1,n-m
!          ho=xa(i)-x
!          hp=xa(i+m)-x
!          w=c(i+1)-d(i)
!          den=ho-hp
!          if(den.eq.0.0d0) then
!              write(*,*)'polint data=',xa,ya,n,x
!              stop 'failure in polint'
!		  endif
!          den=w/den
!          d(i)=hp*den
!          c(i)=ho*den
!        enddo
!        if (2*ns.lt.n-m)then
!          dy=c(ns+1)
!        else
!          dy=d(ns)
!          ns=ns-1
!        endif
!        y=y+dy
!      enddo
!      return
!      END

!	DOUBLE PRECISION FUNCTION rtsec(func,x1,x2,xacc)
!	IMPLICIT NONE
!	DOUBLE PRECISION, INTENT(IN) :: x1,x2,xacc
!	INTERFACE
!		DOUBLE PRECISION FUNCTION func(x)
!		IMPLICIT NONE
!		DOUBLE PRECISION, INTENT(IN) :: x
!		END FUNCTION func
!	END INTERFACE
!	INTEGER, PARAMETER :: MAXIT=30
!	INTEGER :: j
!	DOUBLE PRECISION :: dx,f,fl,xl
!	fl=func(x1)
!	f=func(x2)
!	if (abs(fl) < abs(f)) then
!		rtsec=x1
!		xl=x2
!		call swap(fl,f)
!	else
!		xl=x1
!		rtsec=x2
!	end if
!	do j=1,MAXIT
!		dx=(xl-rtsec)*f/(f-fl)
!		xl=rtsec
!		fl=f
!		rtsec=rtsec+dx
!		f=func(rtsec)
!		if (abs(dx) < xacc .or. f == 0.0d0) RETURN
!	end do
!	call Err_MSG('rtsec: exceed maximum iterations')
!	END FUNCTION rtsec

!	FUNCTION rtsafe(funcd,x1,x2,xacc)
!	IMPLICIT NONE
!	DOUBLE PRECISION, INTENT(IN) :: x1,x2,xacc
!	DOUBLE PRECISION :: rtsafe
!	INTERFACE
!		SUBROUTINE funcd(x,fval,fderiv)
!		IMPLICIT NONE
!		DOUBLE PRECISION, INTENT(IN) :: x
!		DOUBLE PRECISION, INTENT(OUT) :: fval,fderiv
!		END SUBROUTINE funcd
!	END INTERFACE
!	INTEGER, PARAMETER :: MAXIT=100
!	INTEGER :: j
!	DOUBLE PRECISION :: df,dx,dxold,f,fh,fl,temp,xh,xl
!	call funcd(x1,fl,df)
!	call funcd(x2,fh,df)
!	if ((fl > 0.0 .and. fh > 0.0) .or. &
!		(fl < 0.0 .and. fh < 0.0)) &
!	call Err_MSG('root must be bracketed in rtsafe')
!	if (fl == 0.0) then
!		rtsafe=x1
!		RETURN
!	else if (fh == 0.0) then
!		rtsafe=x2
!		RETURN
!	else if (fl < 0.0) then
!		xl=x1
!		xh=x2
!	else
!		xh=x1
!		xl=x2
!	end if
!	rtsafe=0.5d0*(x1+x2)
!	dxold=abs(x2-x1)
!	dx=dxold
!	call funcd(rtsafe,f,df)
!	do j=1,MAXIT
!		if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) > 0.0 .or. &
!			abs(2.0d0*f) > abs(dxold*df) ) then
!			dxold=dx
!			dx=0.5d0*(xh-xl)
!			rtsafe=xl+dx
!			if (xl == rtsafe) RETURN
!		else
!			dxold=dx
!			dx=f/df
!			temp=rtsafe
!			rtsafe=rtsafe-dx
!			if (temp == rtsafe) RETURN
!		end if
!		if (abs(dx) < xacc) RETURN
!		call funcd(rtsafe,f,df)
!		if (f < 0.0) then
!			xl=rtsafe
!		else
!			xh=rtsafe
!		end if
!	end do
!	call Err_MSG('rtsafe: exceeded maximum iterations')
!	END FUNCTION rtsafe

!      DOUBLE PRECISION FUNCTION expbessi1(x)
!      DOUBLE PRECISION bessi1,x
!      DOUBLE PRECISION ax
!      DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
!      SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
!      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0, &
!           0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
!      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1, &
!        -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,-0.2895312d-1, &
!        0.1787654d-1,-0.420059d-2/
!      ax=abs(x)
!      if (abs(x).lt.3.75d0) then
!        y=(x/3.75D0)**2
!        bessi1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
!        expbessi1=dexp(-ax)*bessi1
!      else
!        y=3.75D0/ax
!!        bessi1=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4    &
!!		       +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
!        expbessi1=(1.0d0/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4    &
!		       +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
!
!        if(x.lt.0.0D0)bessi1=-bessi1
!      endif
!      return
!      END


!DOUBLE PRECISION FUNCTION expbessi0(x)
!DOUBLE PRECISION bessi0,x
!DOUBLE PRECISION ax
!DOUBLE PRECISION p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9,y
!SAVE p1,p2,p3,p4,p5,p6,p7,q1,q2,q3,q4,q5,q6,q7,q8,q9
!DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,  &
!	 1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
!DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,      &
!0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,0.2635537d-1,&
!-0.1647633d-1,0.392377d-2/
!ax=abs(x)
!if (abs(x).lt.3.75d0) then
!  y=(x/3.75d0)**2
!  bessi0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
!  expbessi0=dexp(-ax)*bessi0
!else
!  y=3.75d0/ax
!!        bessi0=(exp(ax)/sqrt(ax))*(q1+y*(q2+y*(q3+y*(q4   &
!!		        +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
!  expbessi0=(1.0d0/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4   &
!		  +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
!
!endif
!return
!END

