!-----------------------------------------------------------------
!用于计算第一类零阶、一阶贝塞尔函数；第二类零阶、一阶贝塞尔函数
!第二类变型零阶、一阶贝塞尔函数
!程序编制：孙葳  2011.12.23
!-----------------------------------------------------------------
MODULE BESSEL_MOD

IMPLICIT NONE

PRIVATE BESSELI0,BESSELI1
PUBLIC DBSJ0,DBSJ1,DBSY0,DBSY1,DBSK0,DBSK1

CONTAINS


!-------------------------------------------------------------------
!计算第一类零阶BESSEL函数BESSELJ0(X)
!-------------------------------------------------------------------
FUNCTION DBSJ0(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5
REAL(8)::b0,b1,b2,b3,b4,b5
REAL(8)::E0,E1,E2,E3,E4
REAL(8)::F0,F1,F2,F3,F4
real(8)::T,Y,z,temp1,temp2,TEMP
REAL(8),parameter::	PI=3.14159265358979

a0=5.7568490574D10;a1=-1.3362590354d10;a2=6.516196407d8;a3=-1.121442418d7;a4=77392.33017;a5=-184.9052456
b0=5.7568490411d10;b1=1.029532985d9;b2=9.494680718d6;b3=59272.64853;b4=267.8532712;b5=1.0

E0=1.0;E1=-0.1098628627D-2;E2=0.2734510407D-4;E3=-0.2073370639D-5;E4=0.2093887211D-6
F0=-0.1562499995D-01;F1=0.1430488765D-03;F2=-0.6911147651D-05;F3=0.7621095161D-06;F4=-0.934935152D-07

T=ABS(X)
TEMP1=0.;TEMP2=0.
X_RESULT=0.0

if(T<8.0) THEN
   Y=X**2
   temp1=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)
   temp2=b0+b1*y+b2*(y**2)+b3*(y**3)+b4*(y**4)+b5*(y**5)
   X_RESULT=temp1/temp2
else if(T>=8.0) then
   z=8.0/T
   Y=Z**2
   TEMP=T-PI/4.0
   TEMP1=E0+E1*y+E2*(y**2)+E3*(y**3)+E4*(y**4)
   TEMP2=F0+F1*y+F2*(y**2)+F3*(y**3)+F4*(y**4)
   X_RESULT=SQRT(2.0/(PI*T))*(TEMP1*COS(TEMP)-Z*TEMP2*SIN(TEMP))
END IF

END FUNCTION


!-------------------------------------------------------------------
!计算第一类一阶BESSEL函数BESSELJ1(X)
!-------------------------------------------------------------------
FUNCTION DBSJ1(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5
REAL(8)::b0,b1,b2,b3,b4,b5
REAL(8)::E0,E1,E2,E3,E4
REAL(8)::F0,F1,F2,F3,F4
real(8)::T,Y,z,temp1,temp2,TEMP
REAL(8),parameter::	PI=3.14159265358979

a0=72362614232.0;a1=-7895059235.0;a2=242396853.1;a3=-2972611.439;a4=15704.4826;a5=-30.16036606
b0=144725228443.0;b1=2300535178.0;b2=18583304.74;b3=99447.43394;b4=376.9991397;b5=1.0

E0=1.0;E1=0.183105D-02;E2=-0.3516396496D-04;E3=0.2457520174D-05;E4=-0.240337019D-06
F0=0.4687499995D-01;F1=-0.2002690873D-03;F2=0.8449199096D-05;F3=-0.88228987D-06;F4=0.105787412D-06

T=ABS(X)
TEMP1=0.;TEMP2=0.
X_RESULT=0.0

if(T<8.0) THEN
   Y=X**2
   temp1=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)
   temp2=b0+b1*y+b2*(y**2)+b3*(y**3)+b4*(y**4)+b5*(y**5)
   X_RESULT=X*temp1/temp2
else if(T>=8.0) then
   z=8.0/T
   Y=Z**2
   TEMP=T-3.0*PI/4.0
   TEMP1=E0+E1*y+E2*(y**2)+E3*(y**3)+E4*(y**4)
   TEMP2=F0+F1*y+F2*(y**2)+F3*(y**3)+F4*(y**4)
   X_RESULT=SQRT(2.0/(PI*T))*(TEMP1*COS(TEMP)-Z*TEMP2*SIN(TEMP))
   IF(X<0.0)X_RESULT=-X_RESULT
END IF

END FUNCTION


!-------------------------------------------------------------------
!计算第二类零阶BESSEL函数BESSELY0(X)
!-------------------------------------------------------------------
FUNCTION DBSY0(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5
REAL(8)::b0,b1,b2,b3,b4,b5
REAL(8)::E0,E1,E2,E3,E4
REAL(8)::F0,F1,F2,F3,F4
real(8)::Y,z,temp1,temp2,TEMP
! REAL(8)::BESSELJ0
REAL(8),parameter::	PI=3.14159265358979

a0=-2.957821389D+09;a1=7.062834065D+09;a2=-5.123598036D+08;a3=1.087988129D+07;a4=-8.632792757D+04;a5=2.284622733D+02
b0=4.0076544269D+10;b1=7.452499648D+08;b2=7.189466438D+06;b3=4.74472647D+04;b4=2.261030244D+02;b5=1.0

E0=1.0;E1=-0.1098628627D-02;E2=0.2734510407D-04;E3=-0.2073370639D-05;E4=0.2093887211D-06
F0=-0.1562499995D-01;F1=0.1430488765D-03;F2=-0.6911147651D-05;F3=0.7621095161D-06;F4=-0.934935152D-07

! T=ABS(X)
TEMP1=0.;TEMP2=0.
X_RESULT=0.0

if(X<8.0) THEN
   Y=X**2
   temp1=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)
   temp2=b0+b1*y+b2*(y**2)+b3*(y**3)+b4*(y**4)+b5*(y**5)
   X_RESULT=temp1/temp2+2.0/PI*DBSJ0(X)*LOG(X)
else if(X>=8.0) then
   z=8.0/X
   Y=Z**2
   TEMP=X-PI/4.0
   TEMP1=E0+E1*y+E2*(y**2)+E3*(y**3)+E4*(y**4)
   TEMP2=F0+F1*y+F2*(y**2)+F3*(y**3)+F4*(y**4)
   X_RESULT=SQRT(2.0/(PI*X))*(TEMP1*SIN(TEMP)+Z*TEMP2*COS(TEMP))
END IF

END FUNCTION


!-------------------------------------------------------------------
!计算第二类一阶BESSEL函数BESSELY1(X)
!-------------------------------------------------------------------
FUNCTION DBSY1(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5
REAL(8)::b0,b1,b2,b3,b4,b5,b6
REAL(8)::E0,E1,E2,E3,E4
REAL(8)::F0,F1,F2,F3,F4
real(8)::Y,z,temp1,temp2,TEMP
! REAL(8)::BESSELJ1
REAL(8),parameter::	PI=3.14159265358979

a0=-4.900604943D+12;a1=1.27527439D+12;a2=-5.153438139D+10;a3=7.349264551D+08;a4=-4.237922726D+06;a5=8.511937935D+03
b0=2.49958057D+13;b1=4.244419664D+11;b2=3.733650367D+09;b3=2.245904002D+07;b4=1.02042605D+05;b5=3.549632885D+02;b6=1.0

E0=1.0;E1=0.183105D-02;E2=-0.3516396496D-04;E3=0.2457520174D-05;E4=-0.240337019D-06
F0=0.4687499995D-01;F1=-0.2002690873D-03;F2=0.8449199096D-05;F3=-0.88228987D-06;F4=0.105787412D-06

! T=ABS(X)
TEMP1=0.;TEMP2=0.
X_RESULT=0.0

if(X<8.0) THEN
   Y=X**2
   temp1=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)
   temp2=b0+b1*y+b2*(y**2)+b3*(y**3)+b4*(y**4)+b5*(y**5)+b6*(Y**6)
   X_RESULT=X*temp1/temp2+2.0/PI*(DBSJ1(X)*LOG(X)-1.0/X)
else if(X>=8.0) then
   z=8.0/X
   Y=Z**2
   TEMP=X-3.0*PI/4.0
   TEMP1=E0+E1*y+E2*(y**2)+E3*(y**3)+E4*(y**4)
   TEMP2=F0+F1*y+F2*(y**2)+F3*(y**3)+F4*(y**4)
   X_RESULT=SQRT(2.0/(PI*X))*(TEMP1*SIN(TEMP)+Z*TEMP2*COS(TEMP))
END IF

END FUNCTION


!-------------------------------------------------------------------
!计算变型第一类零阶BESSEL函数BESSELI0(X)
!-------------------------------------------------------------------
FUNCTION BESSELI0(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5,a6
REAL(8)::E0,E1,E2,E3,E4,E5,E6,E7,E8
real(8)::T,Y,TEMP

a0=1.0;a1=3.5156229;a2=3.0899424;a3=1.2067492;a4=0.2659732;a5=0.0360768;a6=0.0045813

E0=0.39894228;E1=0.01328592;E2=0.00225319;E3=-0.00157565;E4=0.00916281
E5=-0.02057706;E6=0.02635537;E7=-0.01647663;E8=0.00392377

T=ABS(X)
TEMP=0.
X_RESULT=0.0

if(T<3.75) THEN
   Y=(X/3.75)**2
   X_RESULT=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)+a6*(y**6)
else if(T>=3.75) then
   Y=3.75/T
   TEMP=E0+E1*y+E2*(y**2)+E3*(y**3)+E4*(y**4)+E5*(y**5)+E6*(y**6)+E7*(Y**7)+E8*(Y**8)
   X_RESULT=TEMP*EXP(T)/SQRT(T)
END IF

END FUNCTION

!-------------------------------------------------------------------
!计算变型第一类一阶BESSEL函数BESSELI1(X)
!-------------------------------------------------------------------
FUNCTION BESSELI1(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5,a6
REAL(8)::E0,E1,E2,E3,E4,E5,E6,E7,E8
real(8)::T,Y,TEMP

a0=0.5;a1=0.87890594;a2=0.51498869;a3=0.15084934;a4=0.02658773;a5=0.00301532;a6=0.00032411

E0=0.39894228;E1=-0.03988024;E2=-0.00362018;E3=0.00163801;E4=-0.01031555
E5=0.02282967;E6=-0.02895312;E7=0.01787654;E8=-0.00420059

T=ABS(X)
TEMP=0.
X_RESULT=0.0

if(T<3.75) THEN
   Y=(X/3.75)**2
   X_RESULT=X*(a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)+a6*(y**6))
else if(T>=3.75) then
   Y=3.75/T
   TEMP=E0+E1*y+E2*(y**2)+E3*(y**3)+E4*(y**4)+E5*(y**5)+E6*(y**6)+E7*(Y**7)+E8*(Y**8)
   X_RESULT=TEMP*EXP(T)/SQRT(T)
   IF(X<0.0)X_RESULT=-X_RESULT
END IF

END FUNCTION


!-------------------------------------------------------------------
!计算变型第二类零阶BESSEL函数BESSELK0(X)
!-------------------------------------------------------------------
FUNCTION DBSK0(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5,A6
REAL(8)::b0,b1,b2,b3,b4,b5,B6

real(8)::Y,TEMP
! REAL(8)::BESSELI0
REAL(8),parameter::	PI=3.14159265358979

a0=-0.57721566;a1=0.4227842;a2=0.23069756;a3=0.0348859;a4=0.00262698;a5=0.0001075;A6=0.0000074
b0=1.25331414;b1=-0.07832358;b2=0.02189568;b3=-0.01062446;b4=0.00587872;b5=-0.0025154;B6=0.00053208


! T=ABS(X)
TEMP=0.
X_RESULT=0.0

if(X<=2.0) THEN
   Y=X**2/4.0
   temp=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)+A6*(Y**6)
   X_RESULT=temp-BESSELI0(X)*LOG(X/2.0)
else if(X>2.0) then
   Y=2.0/X
   TEMP=B0+B1*y+B2*(y**2)+B3*(y**3)+B4*(y**4)+b5*(y**5)+b6*(Y**6)
   X_RESULT=EXP(-X)/SQRT(X)*TEMP
END IF

END FUNCTION


!-------------------------------------------------------------------
!计算变型第二类一阶BESSEL函数BESSELK1(X)
!-------------------------------------------------------------------
FUNCTION DBSK1(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::a0,a1,a2,a3,a4,a5,A6
REAL(8)::b0,b1,b2,b3,b4,b5,B6

real(8)::Y,TEMP
! REAL(8)::BESSELI1
REAL(8),parameter::	PI=3.14159265358979

a0=1.0;a1=0.15443144;a2=-0.67278579;a3=-0.18156897;a4=-0.01919402;a5=-0.00110404;A6=-0.00004686
b0=1.25331414;b1=0.23498619;b2=-0.0365562;b3=0.01504268;b4=-0.00780353;b5=0.00325614;B6=-0.00068245


! T=ABS(X)
TEMP=0.
X_RESULT=0.0

if(X<=2.0) THEN
   Y=X**2/4.0
   temp=a0+a1*y+a2*(y**2)+a3*(y**3)+a4*(y**4)+a5*(y**5)+A6*(Y**6)
   X_RESULT=1.0/X*TEMP+BESSELI1(X)*LOG(X/2.0)
else if(X>2.0) then
   Y=2.0/X
   TEMP=B0+B1*y+B2*(y**2)+B3*(y**3)+B4*(y**4)+b5*(y**5)+b6*(Y**6)
   X_RESULT=EXP(-X)/SQRT(X)*TEMP
END IF

END FUNCTION


!------------------------------------------------------
END MODULE BESSEL_MOD