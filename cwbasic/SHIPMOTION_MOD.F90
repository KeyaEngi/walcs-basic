MODULE SHIPMOTION_MOD


USE INOUTACCESS_MOD !文件路径变量
USE MESH_MOD,ONLY:SL,SB,ST,VOL,HRM
USE CAL_MOD,ONLY:ROllDAMPING
! USE PRINT_MOD,ONLY:ANSMETHOD
USE ENVIRONMENT_MOD,ONLY:ROU,G0,MASSMATRIX,PI


IMPLICIT NONE
PRIVATE
PUBLIC::SHIPMOTION

CONTAINS

!************************************************************************************
!程序功能：求解船舶运动
!          MASSMATRIX(1:6,1:6),HRM(1:6,1:6)为船舶的质量矩阵和静恢复力矩阵
!          F(1:2,1:6,1:9)
!输出参数：motion(1:6,1:2)为船舶六自由度运动的实部和虚部 
!************************************************************************************
subroutine shipmotion(Beta,omeI,F,motion,Fn0,M4x,bv44)
implicit none

real(8),intent(in)::omeI,Beta,FN0
real(8),dimension(1:2,1:6,1:9),intent(in)::F 
real(8),dimension(1:2,1:6),intent(out)::motion
real(8),intent(out),dimension(1:2)::M4x
real(8),intent(out)::bv44

real(8),dimension(1:6,1:7)::A,B
real(8)	::B44,FN1
integer(4)::i,j
! real(8)::a22

!纵向运动
motion=0.0
B44=F(2,4,4)
FN1=FN0/sqrt(g0*sl)
!F(2,3,3)=F(2,3,3)+
do J=1,3
	do I=1,3
! 		if(flag==0)then
			A(i,j)=HRM(2*i-1,2*j-1)-omeI**2*(F(1,2*i-1,2*j-1)+MASSMATRIX(2*i-1,2*j-1))
			A(i,j+3)=-omeI*F(2,2*i-1,2*j-1)
			A(i+3,j)=-A(i,j+3)
			A(i+3,j+3)=A(i,j)
! 		else
! 			A(i,j)=CC(2*i-1,2*j-1)-ome**2*(fitF(2*i-1,2*j-1,1)+MM(2*i-1,2*j-1))
! 			A(i,j+3)=-ome*fitF(2*i-1,2*j-1,2)
! 			A(i+3,j)=-A(i,j+3)
! 			A(i+3,j+3)=A(i,j)
! 		endif
	       A(i,7)=F(1,2*i-1,9)
	       A(i+3,7)=F(2,2*i-1,9)
	enddo
enddo

b=a
call Gauss(A,6,7)

do i=1,3
	motion(1,2*i-1)=A(i,7)
	motion(2,2*i-1)=A(i+3,7)
enddo
!横向运动
if(abs(sin(Beta))>0.001)then
	do J=1,3
		do I=1,3
! 			if(flag==0)then
				A(i,j)=HRM(2*i,2*j)-omeI**2*(F(1,2*i,2*j)+MASSMATRIX(2*i,2*j))
				A(i,j+3)=-omeI*F(2,2*i,2*j)
				A(i+3,j)=-A(i,j+3)
				A(i+3,j+3)=A(i,j)
! 			else
! 				A(i,j)=CC(2*i,2*j)-ome**2*(fitF(2*i,2*j,1)+MM(2*i,2*j))
! 				A(i,j+3)=-ome*fitF(2*i,2*j,2)
! 				A(i+3,j)=-A(i,j+3)
! 				A(i+3,j+3)=A(i,j)
! 			endif
		        A(i,7)=F(1,2*i,9)
		        A(i+3,7)=F(2,2*i,9)
		enddo

	enddo
	B=A
	call Gauss(A,6,7)						
	
! 	a22=F(1,2,2)+MASSMATRIX(2,2)

	if(ROllDAMPING/=0)call RollMotion(omeI,a,b,F,Fn1,HRM(4,4),MASSMATRIX(4,4),B44,M4x,bv44)
	do i=1,3
		motion(1,2*i)=A(i,7)
		motion(2,2*i)=A(i+3,7)
	enddo
endif


end subroutine shipmotion
!----------------------------------------------------------------



!----------------------------------------------------------------
!JC=1--->rma,rmb
!   2--->rma,rmb,ek0
!   3--->Empirical Method of Miller
!   4--->Empirical Method for Barges
!   5
!----------------------------------------------------------------
subroutine RollMotion(omeI,a,b,F,Fn0,C44,M44,B440,M4x,b44)
implicit none

real(8),intent(in)::OMEI,FN0
real(8),intent(in),dimension(1:2,1:6,1:9)::F 
real(8),intent(in)	::	M44,C44,B440
real(8),intent(inout),dimension(1:6,1:7)::A,B
real(8),intent(out),dimension(1:2)::M4x
real(8),intent(out)::b44

real(8)::rma,rmb,rat0,ek0,Bcoef
real(8)::Lbk,hbk,rb !舭龙骨长度，宽度，舭部距中心线距离

real(8)::temp1,temp2,rat,esp
real(8)::dvel,ds,av,Cv,Abk,k1,k2,hx,Cb
CHARACTER(LEN=80)::RDMINPUT
integer::i
! PI=4.0*atan(1.0)
esp=1e-6
M4x=0.0
if(ROllDAMPING/=4)	then    !读入横摇阻尼修正文件
    RDMINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.RDM'
	open(13,file=RDMINPUT,status='old')
	call	IUTMP(13)
ENDIF

temp1=sqrt(A(2,7)**2+A(5,7)**2)

if(ROllDAMPING==1)then

    read(13,*)	rma,rmb
	
	do
		b44=rma+8.0*rmb/(3.0*PI)*omeI*temp1
		B(2,5)=-omeI*b44
		B(5,2)=-B(2,5)
		do i=1,3
			B(i,7)=F(1,2*i,9)
			B(i+3,7)=F(2,2*i,9)
		enddo
		A=B
		call Gauss(A,6,7)
		temp2=sqrt(A(2,7)**2+A(5,7)**2)
		if(abs(temp1-temp2)<esp)exit

		temp1=temp2
	enddo
	m4x(1)=-b44*A(5,7)*omeI+b440*A(5,7)*omeI
	m4x(2)=b44*A(2,7)*omeI-b440*A(2,7)*omeI

elseif(ROllDAMPING==2)then

   	read(13,*)	rma,rmb,ek0

    rat0=0.8*(1-exp(-10*Fn0))
	rat=rat0*(ek0/omeI)**2+1
	dvel=0.5
	av=omeI*omeI*temp1*temp1
	ds=av-dvel

	do
		dvel=dvel+0.5*ds
		b44=rma+8.0*rmb/(3.0*PI)*sqrt(dvel)
		B(2,5)=-omeI*b44*rat
		B(5,2)=-B(2,5)
		do i=1,3
			B(i,7)=F(1,2*i,9)
			B(i+3,7)=F(2,2*i,9)
		enddo
		A=B
		call Gauss(A,6,7)
		temp2=sqrt(A(2,7)**2+A(5,7)**2)
		av=omeI*omeI*temp2*temp2
		ds=av-dvel
		if(av==0.0.or.abs(ds)<0.005*dvel+esp)exit

	enddo
    m4x(1)=-b44*rat*A(5,7)*omeI+b440*rat*A(5,7)*omeI
	m4x(2)=b44*rat*A(2,7)*omeI-b440*rat*A(2,7)*omeI

elseif(ROllDAMPING==3)then !Empirical Method of Miller
    
	read(13,*)	Lbk,hbk,rb

	dvel=0.5
	av=temp1*temp1
	ds=av-dvel
	hx=C44/(rou*g0*vol)
	Cb=vol/(sl*sb*st)
	Abk=Lbk*hbk
	if(hx<0.0)hx=0.1
	Cv=4.85-3.0*dsqrt(hx)
	if(Cv<0.0)Cv=0.0
	k1=Cv*0.00085*sl/sb*dsqrt(sl/hx)*(Fn0/Cb+(Fn0/Cb)**2+2.0*(Fn0/Cb)**3)
	k2=19.25*(Abk*dsqrt(Lbk/rb+0.0024*sl*sb))*rb**3/(sl*sb**3*st*Cb)
	do

		dvel=dvel+0.5*ds
		b44=(k1+k2*dsqrt(dsqrt(dvel)))*C44*2.0/omeI
		B(2,5)=-omeI*b44
		B(5,2)=-B(2,5)
		do i=1,3
			B(i,7)=F(1,2*i,9)
			B(i+3,7)=F(2,2*i,9)
		enddo
		A=B
		call Gauss(A,6,7)
		temp2=sqrt(A(2,7)**2+A(5,7)**2)
		av=temp2*temp2
		ds=av-dvel
		if(av==0.0.or.abs(ds)<0.005*dvel+esp)exit

	enddo
	m4x(1)=-b44*A(5,7)*omeI+b440*A(5,7)*omeI
	m4x(2)=b44*A(2,7)*omeI-b440*A(2,7)*omeI
elseif(ROllDAMPING==4)then !Empirical Method for Barges
	dvel=0.5
	av=temp1*temp1
	ds=av-dvel
	k1=0.0013*(sb/st)**2
	k2=0.5
	do
		dvel=dvel+0.5*ds
		b44=(k1+k2*dsqrt(dvel))*C44*2.0/omeI
		B(2,5)=-omeI*b44
		B(5,2)=-B(2,5)
		do i=1,3
			B(i,7)=F(1,2*i,9)
			B(i+3,7)=F(2,2*i,9)
		enddo
		A=B
		call Gauss(A,6,7)
		temp2=sqrt(A(2,7)**2+A(5,7)**2)
		av=temp2*temp2
		ds=av-dvel
		if(av==0.0.or.abs(ds)<0.005*dvel+esp)exit
	enddo
	m4x(1)=-b44*A(5,7)*omeI+b440*A(5,7)*omeI
	m4x(2)=b44*A(2,7)*omeI-b440*A(2,7)*omeI
elseif(ROllDAMPING==5)	then
    read(13,*)	Bcoef
    
	b44=Bcoef*2.0*sqrt((M44+abs(F(1,4,4)))*C44)   !F(1,4,4)应为附加横摇惯性矩，原程序(a(4,4))有误--By SUN
! 	b44=Bcoef*2.0*sqrt((M44+abs(a(4,4)))*C44)
	B(2,5)=-omeI*b44+B(2,5)
	B(5,2)=-B(2,5)
	do i=1,3
		B(i,7)=F(1,2*i,9)
		B(i+3,7)=F(2,2*i,9)
	enddo
	A=B
	call Gauss(A,6,7)
	m4x(1)=-b44*A(5,7)*omeI
	m4x(2)=b44*A(2,7)*omeI

endif




end subroutine RollMotion
!--------------------------------------------------------------

END MODULE SHIPMOTION_MOD