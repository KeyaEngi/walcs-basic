!*************************************************************************************
!模块功能：计算格林函数G及其梯度、法向导数沿面元的积分
!输入参数：SRPHI为基本源第一下标(场点)在第二下标(源点)下的诱导速度及单层势，
!          DRPHI为上述诱导速度在场点面元的法向投影,wkjNum为波数数目,wkj为波数，Pnu为
!          w0**2/g,OrderNum为高斯积分节点数目
!		   id:是否计算积分项(1---是,0---否)
!输出参数：GcS(i,j,1:4),GsS(i,j,1:4)分别为格林函数梯度及其本身沿面元的积分的实部和虚部,
!          GcD和GsD分别为法向导数沿面元的积分的实部和虚部
!调用子程序：FG1,FG2,Flags,GR,GAUSSIN,Gcfifpoint,INFG,infiFxy,fiFxy,struveH0,struveH1
!模块编制：孙葳   2011.12.9
!*************************************************************************************
MODULE GREEN_MOD
USE INOUTACCESS_MOD
! use msimsl	!调用除J1外的函数，因为dbsJ1计算不准确,应用portlib中的dbesj1(x)
! use portlib	!调用dbesj1(x)  
USE BESSEL_MOD  !计算bessel函数J0,J1,Y0,Y1,K0,K1
USE DEI_MOD     !计算指数积分DEI(Y)

USE ENVIRONMENT_MOD,ONLY:G0,PI,H
USE MESH_MOD,ONLY:NPANEL,XAV,EA,E,NB,Factor,TarNum,Space
USE CAL_MOD,ONLY:WATERDEPTH
USE GAUSSIN_MOD

IMPLICIT NONE 
PRIVATE
PUBLIC::coef,FG1,FG2,INFG

CONTAINS
!--------------------------------------------------------------------
subroutine coef(wkjNum,wkj,Pnu,SRPHI,DRPHI,GcS,GsS,GcD,GsD,G0DC,G0DS)
implicit none
CHARACTER(LEN=500)::BASS
integer(4)::count    !为数组赋值的工具。edited by 王川    2012.09.16
integer(4),intent(in)::wkjNum
real(8),INTENT(IN),dimension(1:NPANEL,1:NPANEL,1:4,1:Factor)::SRPHI	!诱导速度及单层势
real(8),INTENT(IN),dimension(1:NPANEL,1:NPANEL,1:Factor)::DRPHI	  !诱导速度法向投影
real(8),intent(in),dimension(1:wkjNum)::wkj
real(8),intent(inout)::Pnu
real(8),intent(out),dimension(1:NPanel,1:NPanel)::GcD,GsD
real(8),intent(out),dimension(1:NPanel,1:NPanel,1:4)::GcS,GsS
real(8),intent(out),dimension(1:NPanel,1:NPanel,1:3,1:Factor)::G0DC,G0DS
INTEGER(4)::OrderNum
real(8),dimension(1:4)::Sc,Ss 
real(8),dimension(0:2)::Gc,Gs 
real(8),allocatable,dimension(:)::xk,wk
integer(4)::i,j,ij,midI,L     !用来控制对称性的循环执行
real(8)::x,y,R,pz,qz
real(8),allocatable::x00(:,:)
OrderNum=15
allocate(xk(1:OrderNum),wk(1:OrderNum))
allocate(x00(1:3,1:Npanel))
L=0
if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
	call GaussIn(OrderNum,xk,wk,1) !获得高斯积分的积分节点和权函数
endif

do count=1,TarNum,Space
select case(count)
       case(1)
                   x00=xav
        case(2)
              x00(1,:)=-xav(1,:)
              x00(2,:)=xav(2,:)   
       case(3)
             x00(1,:)=-xav(1,:)
             x00(2,:)=-xav(2,:)
        case(4)
            x00(1,:)=xav(1,:)
            x00(2,:)=-xav(2,:)
 end select
L=L+1
do i=1,NPanel
	if((trim(adjustl(WATERDEPTH))=="FINITE"))then  !有限水深   liu修改
		midI=1                   !liu修改
	else      !无限水深           liu修改
		midI=i                   !liu修改
	endif                        !liu修改
	do j=midI,NPanel

			Sc=0.0;Ss=0.0
			x=xav(1,i)-x00(1,j)
			y=xav(2,i)-x00(2,j)
			R=sqrt(x*x+y*y)
			pz=xav(3,i)
			qz=xav(3,j)
		!---------------------------------------------------------------------
			if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
				ij=1
				if(R/H>0.5)ij=2
				if(i==j.and.count==1)ij=1
				if(ij==1)then
					call FG1(OrderNum,xk,wk,wkj(1),Pnu,R,pz,qz,Gc,Gs) !积分形式
				elseif(ij==2)then
					call FG2(wkjNum,wkj,Pnu,R,pz,qz,Gc,Gs) !级数形式
				endif
		 !---------------------------------------------------------------------
			elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深Green函数计算
				call inFG(pnu,R,pz,qz,Gc,Gs)
			endif
		 !------------------------------------------------------------
			Sc(4)=Gc(0)   !Sc(4)=Gc(0)*ea(j)   liu修改
			Ss(4)=Gs(0)   !Ss(4)=Gs(0)*ea(j)   liu修改
			if(abs(R)>1e-6)then
				Sc(1)=Gc(1)*x/R  !Sc(1)=Gc(1)*ea(j)*x/R   liu修改  
				Sc(2)=Gc(1)*y/R  !Sc(2)=Gc(1)*ea(j)*y/R   liu修改
				Ss(1)=Gs(1)*x/R  !Ss(1)=Gs(1)*ea(j)*x/R   liu修改
				Ss(2)=Gs(1)*y/R  !Ss(2)=Gs(1)*ea(j)*y/R   liu修改 
			endif

			Sc(3)=Gc(2)  !Sc(3)=Gc(2)*ea(j)   liu修改
			Ss(3)=Gs(2)  !Ss(3)=Gs(2)*ea(j)   liu修改
			!---------------------------------------------------
			GcS(i,j,1:4)=Sc(1:4)*ea(j)  !GcS(i,j,1:4)=Sc(1:4)  liu修改
			GsS(i,j,1:4)=Ss(1:4)*ea(j)  !GsS(i,j,1:4)=Ss(1:4)  liu修改
			GcD(i,j)=dot_product(GcS(i,j,1:3),e(:,3,i)) !GcD(i,j)=dot_product(Sc(1:3),e(:,3,i))
			GsD(i,j)=dot_product(GsS(i,j,1:3),e(:,3,i)) !GsD(i,j)=dot_product(Ss(1:3),e(:,3,i))
			if(i/=j.and.(trim(adjustl(WATERDEPTH)))=="INFINITE")then
select case(COUNT)
       case(1)
				GcS(j,i,1:2)=-Sc(1:2)*ea(i);GcS(j,i,3:4)=Sc(3:4)*ea(i)         !liu修改
				GsS(j,i,1:2)=-Ss(1:2)*ea(i);GsS(j,i,3:4)=Ss(3:4)*ea(i)         !liu修改
				GcD(j,i)=dot_product(GcS(j,i,1:3),e(:,3,j))                    !liu修改
			    GsD(j,i)=dot_product(GsS(j,i,1:3),e(:,3,j))                    !liu修改
        case(2)
				GcS(j,i,1)=Sc(1)*ea(i);GcS(j,i,2)=-Sc(2)*ea(i)
				GsS(j,i,1)=Ss(1)*ea(i);GsS(j,i,2)=-Ss(2)*ea(i)
				GcS(j,i,3)=Sc(3)*ea(i);GcS(j,i,4)=Sc(4)*ea(i)
				GsS(j,i,3)=Ss(3)*ea(i);GsS(j,i,4)=Ss(4)*ea(i)
				GcD(j,i)=dot_product(GcS(j,i,1:3),e(:,3,j))                    !liu修改
			    GsD(j,i)=dot_product(GsS(j,i,1:3),e(:,3,j))                    !liu修改
       case(3)
				GcS(j,i,1)=Sc(1)*ea(i);GcS(j,i,2)=Sc(2)*ea(i)
				GsS(j,i,1)=Ss(1)*ea(i);GsS(j,i,2)=Ss(2)*ea(i)
				GcS(j,i,3)=Sc(3)*ea(i);GcS(j,i,4)=Sc(4)*ea(i)
				GsS(j,i,3)=Ss(3)*ea(i);GsS(j,i,4)=Ss(4)*ea(i)
				GcD(j,i)=dot_product(GcS(j,i,1:3),e(:,3,j))  
			    GsD(j,i)=dot_product(GsS(j,i,1:3),e(:,3,j))   
        case(4)
				GcS(j,i,1)=-Sc(1)*ea(i);GcS(j,i,2)=Sc(2)*ea(i)
				GsS(j,i,1)=-Ss(1)*ea(i);GsS(j,i,2)=Ss(2)*ea(i)
				GcS(j,i,3)=Sc(3)*ea(i);GcS(j,i,4)=Sc(4)*ea(i)
				GsS(j,i,3)=Ss(3)*ea(i);GsS(j,i,4)=Ss(4)*ea(i)
				GcD(j,i)=dot_product(GcS(j,i,1:3),e(:,3,j)) 
			    GsD(j,i)=dot_product(GsS(j,i,1:3),e(:,3,j))   			
 end select 

endif
			!-------------------------------------------
		if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
			if(ij==1)then
				GcS(i,j,1:4)=GcS(i,j,1:4)+SRPHI(i,j,1:4,L)
				GcD(i,j)=GcD(i,j)+DRPHI(i,j,L)
			elseif(i==j.and.ij==2.and.count==1)then
				GcD(i,j)=GcD(i,j)+2*PI
			endif
		!---------------------------------------------
		elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
			GcS(i,j,1:4)=GcS(i,j,1:4)+SRPHI(i,j,1:4,L)
			GcD(i,j)=GcD(i,j)+DRPHI(i,j,L)
			if(i/=j)then
				GcS(j,i,1:4)=GcS(j,i,1:4)+SRPHI(j,i,1:4,L)
				GcD(j,i)=GcD(j,i)+DRPHI(j,i,L)
			endif
		endif
		!------------------------------------------
		if(i==j.and.count==1.and.i>NB.and.(trim(adjustl(WATERDEPTH))=="INFINITE"))then
			GcD(i,j)=GcD(i,j)-6*PI
		elseif(i==j.and.count==1.and.i>NB.and.(trim(adjustl(WATERDEPTH))=="FINITE"))then
			GcD(i,j)=GcD(i,j)-6*PI
		endif
	enddo
enddo
			G0Dc(:,:,1,L)=GcS(:,:,4)				
			G0Ds(:,:,1,L)=GsS(:,:,4)
			G0Dc(:,:,3,L)=GcD(:,:)
			G0Ds(:,:,3,L)=GsD(:,:)
			G0Dc(:,:,2,L)=GcS(:,:,1)
			G0Ds(:,:,2,L)=GsS(:,:,1)
end do
!BASS=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.BASS'
!    OPEN (30021,FILE=BASS)	
!do i=1,1  !p
!if(XAV(2,i)>=0)then
!      do j=1,NB !q
!      if(XAV(2,J)>=0)then
!WRITE(30021,"(9F15.5)")XAV(1,I),XAV(2,I),XAV(3,I),XAV(1,J),XAV(2,J),XAV(3,J),G0Dc(i,j,3,1),G0Ds(i,j,3,1)
!end if
!END DO
!end if
!END DO
!
!do i=1,1 !p
!if(XAV(2,i)>=0)then
!      do j=1,NB !q
!      if(XAV(2,J)<0)then
!WRITE(30021,"(9F15.5)")XAV(1,I),XAV(2,I),XAV(3,I),XAV(1,J),-XAV(2,J),XAV(3,J),G0Dc(i,j,3,1),G0Ds(i,j,3,1)
!end if
!END DO
!end if
!END DO
deallocate(xk,wk)
end subroutine coef

!---------------------------------------------------------------------------------


!------------------------------------------------------------------
!按积分形式计算格林函数G(Gc+iGs)及其对R和z的偏导数
!------------------------------------------------------------------
subroutine FG1(OrderNum,xk,wk,k0,Pnu,R,pz,qz,Gc,Gs)
implicit none

integer(4),intent(in)::OrderNum	 !积分节点数目
real(8),intent(in),dimension(1:OrderNum)::xk,wk  !积分节点和权函数
real(8),intent(in)::k0,PNU,R,pz,qz
real(8),intent(out),dimension(0:2)::Gc,Gs !分别为(Gc+iGs)及其对R和z的偏导数
real(8)::C0,B1,B2,B3,B4,B5,midJ

B1=cosh(k0*qz)
B2=sinh(k0*qz)
B3=tanh(k0*H)
B4=cosh(k0*pz)
B5=sinh(k0*pz)
C0=H/((cosh(k0*H))**2)+Pnu/(k0**2)
C0=2*PI/C0*(B1+B2*B3)
midJ=dbsj0(k0*R)
Gs(0)=-C0*(B4+B5*B3)*midJ
Gs(1)=-C0*(B4+B5*B3)*dbsj1(k0*R)*(-k0)
Gs(2)=-C0*(B5+B4*B3)*midJ*k0
call Flags(OrderNum,xk,wk,k0,Pnu,R,pz,qz,Gc)   !应用高斯-拉盖尔积分计算实部
end subroutine FG1
!---------------------------------------------------------------------------------


!------------------------------------------------------------------
subroutine Flags(OrderNum,xk,wk,k0,PNU,R,pz,qz,Gc)                                             
implicit none

integer(4),intent(in)::OrderNum	   !积分节点数目
real(8),intent(in),dimension(1:OrderNum)::xk,wk	 !积分节点和权函数
real(8),intent(in)::k0,PNU,R,pz,qz
real(8),intent(out),dimension(0:2)::Gc	!分别为(Gc-1/r-1/r2)及其对R和z的偏导数
integer(4)::i
real(8)::x,Fpn,Fpn0,F0pn,F1pn,F2pn,xm0,Fm1,Fm2,GG,GG0,G00,G1,G2
real(8)::dbF0,dbF1,dbExf0,dbExf1,dbG0,dbG1,dbExg0,dbExg1,dbExg2,dbDeif
real(8)::Hw,Rpq,Zp,Zq,infiF(0:1)
!求解第一项
Gc(0)=+FR0((pz+qz),R)+FR0((qz+2*H-pz),R)+FR0((pz+2*H-qz),R)+FR0((pz+qz+4*H),R)
Gc(1)=-FR1(R,(pz+qz))-FR1(R,(qz+2*H-pz))-FR1(R,(pz+2*H-qz))-FR1(R,(pz+qz+4*H))
Gc(2)=-FR1((pz+qz),R)+FR1((qz+2*H-pz),R)-FR1((pz+2*H-qz),R)-FR1((pz+qz+4*H),R)
!求解第二项
Gc(1)=Gc(1)-(FR2(-pz-qz,R)+FR2(qz+2*H-pz,R)+FR2(pz+2*H-qz,R)+FR2(pz+qz+4*H,R))*2*PNU
Gc(2)=Gc(2)+(FR0(-pz-qz,R)+FR0(qz+2*H-pz,R)-FR0(pz+2*H-qz,R)-FR0(pz+qz+4*H,R))*2*PNU
!求解第三项 
xm0=(tanh(k0*H)+k0*H/(cosh(k0*H)**2))*(1+exp(-2*k0*H))   !分母导数
Fpn=2*PNU+(k0+PNU)*exp(-2*k0*H)
Fpn0=2*PNU*PNU+(k0*k0+3*k0*PNU+2*PNU*PNU)*exp(-2*k0*H)
dbF0=dbsj0(k0*R)
dbF1=dbsj1(k0*R)
dbExf0=exp(-k0*(qz+H-pz))+exp(-k0*(pz+H-qz))+exp(-k0*(qz+pz+3*H))
dbExf1=exp(-k0*(qz+H-pz))-exp(-k0*(pz+H-qz))-exp(-k0*(qz+pz+3*H))
do i=1,OrderNum
	x=xk(i)/H
	Fm1=(xk(i)*(1-exp(-2*xk(i)))/(1+exp(-2*xk(i)))-H*PNU)*(1+exp(-2*xk(i)))
	Fm2=(xk(i)-k0*H)*xm0
	if(abs(Fm1)<1e-8)Fm1=Fm1-0.5e-7
	if(abs(Fm2)<1e-8)Fm2=Fm2-0.5e-7
	GG=2*PNU+(x+PNU)*exp(-2*xk(i))
	GG0=2*PNU*PNU+(x**2+3*x*PNU+2*PNU*PNU)*exp(-2*xk(i))
	dbG0=dbsj0(x*R)
	dbG1=dbsj1(x*R)
	dbExg0=exp(-x*(qz+H-pz))+exp(-x*(pz+H-qz))+exp(-x*(pz+qz+3*H))
	dbExg1=exp(-x*(qz+H-pz))-exp(-x*(pz+H-qz))-exp(-x*(qz+pz+3*H))
	dbExg2=exp(x*(pz+qz+H))
	G00=GG*(dbExg2+dbExg0)*dbG0
	G1=-GG0*(dbExg2+dbExg0)*dbG1
	G2=GG0*(dbExg2+dbExg1)*dbG0
	F0pn=Fpn*(dbExg2*dbG0+dbExf0*dbF0)
	F1pn=-Fpn0*(dbExg2*dbG1+dbExf0*dbF1)
	F2pn=Fpn0*(dbExg2*dbG0+dbExf1*dbF0)
	Gc(0)=Gc(0)+wk(i)*(G00/Fm1-F0pn/Fm2)
 	Gc(1)=Gc(1)+wk(i)*(G1/Fm1-F1pn/Fm2)
	Gc(2)=Gc(2)+wk(i)*(G2/Fm1-F2pn/Fm2)
enddo
dbDeif=exp(-k0*H)*dei(k0*H)
Gc(0)=Gc(0)-dbDeif*Fpn*dbExf0*dbF0/xm0
Gc(1)=Gc(1)+dbDeif*Fpn0*dbExf0*dbF1/xm0
Gc(2)=Gc(2)-dbDeif*Fpn0*dbExf1*dbF0/xm0
Rpq=k0*R;Zp=-(pz+qz)*k0
call fiFxy(Rpq,Zp,infiF)
infiF(1)=-0.5*infiF(1)-FR2(Zp,Rpq)
Gc(0)=Gc(0)+0.5*infiF(0)*Fpn/xm0;Gc(1)=Gc(1)-infiF(1)*Fpn0/xm0;Gc(2)=Gc(2)+0.5*infiF(0)*Fpn0/xm0
!自定义函数
contains
real(8) function FR0(x,y)
real(8)::x,y
if(abs(x)<=1.0e-6.and.y<=1.0e-6)then
	FR0=0.0
else
	FR0=1/(SQRT(x**2+y**2))
endif
end function
real(8) function FR1(x,y)
real(8)::x,y
if(abs(x)<=1.0e-6.and.y<=1.0e-6)then
	FR1=0.0
else
	FR1=x/((SQRT(x**2+y**2))**3)
endif
end function
real(8) function FR2(x,y)
real(8)::x,y
if(y<+1.0e-6)then
	FR2=0.0
else
	FR2=1/y*(1-x/sqrt(x**2+y**2))
endif
end function
end subroutine Flags
!--------------------------------------------------                                           


!------------------------------------------------------------------------- 
!按级数形式计算格林函数G(Gc+iGs)及其对R和z的偏导数
!------------------------------------------------------------------------- 
subroutine FG2(wkjNum,wkj,Pnu,R,pz,qz,Gc,Gs)
implicit none

integer(4),intent(in)::wkjNum	!波数数目
real(8),intent(in),dimension(1:wkjNum)::wkj	 !波数
real(8),intent(in)::PNU,R,pz,qz
real(8),intent(out),dimension(0:2)::Gc,Gs !分别为(Gc+iGs)及其对R和z的偏导数
real(8)::C0,B1,B2,B3,B4,B5,w,CJ,k0,midJ
integer(4)::i,MJ
real(8),dimension(0:2)::G

k0=wkj(1)
B1=cosh(k0*qz)
B2=sinh(k0*qz)
B3=tanh(k0*H)
B4=cosh(k0*pz)
B5=sinh(k0*pz)
C0=H/((cosh(k0*H))**2)+Pnu/(k0**2)
C0=2*PI/C0*(B1+B2*B3)
midJ=dbsj0(k0*R)
Gs(0)=-C0*(B4+B5*B3)*midJ
Gs(1)=-C0*(B4+B5*B3)*dbsj1(k0*R)*(-k0)
Gs(2)=-C0*(B5+B4*B3)*midJ*k0
!-------------------------------------------------------------
midJ=dbsy0(k0*R)
Gc(0)=-C0*(B4+B5*B3)*midJ
Gc(1)=-C0*(B4+B5*B3)*dbsy1(k0*R)*(-k0)
Gc(2)=-C0*(B5+B4*B3)*midJ*k0
G=0.0
MJ=6*H/R
if(MJ>(wkjNum-1))MJ=wkjNum-1
do i=1,MJ
	w=wkj(i+1)
	midJ=dbsk0(w*R)
	CJ=4*(w*w+Pnu*Pnu)/(H*(w*w+Pnu*Pnu)-Pnu)*cos(w*(qz+H))
	G(0)=G(0)+CJ*cos(w*(pz+H))*midJ
	G(1)=G(1)+CJ*cos(w*(pz+H))*dbsk1(w*R)*(-w)
	G(2)=G(2)+CJ*sin(w*(pz+H))*midJ*(-w)
enddo
Gc=Gc+G
end subroutine FG2
!---------------------------------------------------------------------



!****************************************************************************
!本程序计算频域无限水深Green函数中G(p,q)-1/r-1/r1及其对R,z的偏导数部分
!输出参数:Gc(0:2),Gs(0:2)分别为函数本身及对R,z偏导数的实部和虚部
!****************************************************************************
subroutine inFG(pnu,R,pz,qz,Gc,Gs)
implicit none

real(8),intent(in)::pnu,R,pz,qz
real(8),intent(out),dimension(0:2)::Gc,Gs
real(8)::x,y

x=pnu*R
y=-pnu*(pz+qz)
call infiFxy(x,y,Gc)
Gc(0)=pnu*Gc(0)
Gc(1)=pnu*pnu*Gc(1)
Gc(2)=-pnu*pnu*Gc(2)
Gs(0)=-2.0*PI*pnu*dexp(pnu*(pz+qz))*DBSJ0(pnu*R)
Gs(1)=-2.0*PI*pnu*dexp(pnu*(pz+qz))*dbsj1(pnu*R)*(-pnu) !待修改
Gs(2)=-2.0*PI*pnu*dexp(pnu*(pz+qz))*DBSJ0(pnu*R)*pnu
end subroutine
!---------------------------------------------------------------------------------



!****************************************************************************
!本程序计算频域无限水深Green函数中F(x,y)部分
!输出参数:infiF(0:2)分别为函数本身及对x,y的导数
!****************************************************************************
subroutine infiFxy(x,y,infiF)
implicit none

real(8),intent(in)::x,y
real(8),dimension(0:2),intent(out)::infiF
integer(4),parameter::m=4
real(8),parameter::esp=1.0e-6  !for VF alter from 1.0e-12
real(8)::cta,polyvalue,R1,sta,cta2
real(8),dimension(1:m+2)::polycoef,ld
real(8),allocatable,dimension(:)::xx,yy
real(8),allocatable,dimension(:,:)::Fxy,dxFxy
real(8)::Fx0y1,Fx0y2,Fx0y0,slope,dxFx0y1,dxFx0y2,dxFx0y0,zero,Njc,infiFmid1,infiFmid2,infiFmid3,infiFmid(0:1),Cmn(0:6,14),Im(0:5),Im1(0:5),Im2(0:5)
integer(4)::k,i,j,countx,county,nx,ny,i1,i2
real(8)::coef0(0:5,0:5),coef1(0:5,0:5),m0x(0:5),m0y(0:5),ax(2),ay(2),x1,y1

infiF=0.0
zero=0.0    !实型数虚实结合时要用变量，以保持类型一致。(release)
R1=sqrt(x*x+y*y)

if(x>8.and.y>18)then
	cta=y/R1;sta=x/R1;cta2=cta**2
	infiF(0)=(((((9+cta2*(105*cta2-90))/R1+cta*(15*cta2-9))/R1+3*cta2-1.0)/R1+cta)/R1+1.0)/R1
	infiF(1)=((((((945*cta2-630)*cta2+45)/R1+(105*cta2-45)*cta)/R1+15*cta2-3)/R1+3*cta)/R1+1)*sta/(R1**2)
	infiF(0)=-PI*exp(-y)*(struveH0(x)+DBSY0(x))-2*infiF(0)
	infiF(1)=2*PI*exp(-y)*DBSY1(x)+2*infiF(1)
	infiF(2)=-2.0/R1-infiF(0)
	!求积分函数F(x,y)
elseif(abs(x)<1e-5.and.y/=0.0)then
	infiF(0)=-2*exp(-y)*dei(y)
	!求函数F(x,y)对x导数
	infiF(1)=0.0
	!求函数F(x,y)对y导数
	infiF(2)=-2.0/R1-infiF(0)
elseif(abs(x)<1e-5.and.abs(y)<1e-5)then
	infiF=0.0
elseif((x<3.7).and.(y<2))then
	infiF=0.0
	Cmn(0,1)= 0.2500000E+00;Cmn(0,2)= 0.5555556E-01;Cmn(0,3)= 0.1041667E-01;Cmn(0,4)= 0.1666667E-02;Cmn(0,5)= 0.2314815E-03;Cmn(0,6)= 0.2834467E-04;Cmn(0,7)= 0.3100198E-05;Cmn(0,8)= 0.3061924E-06;Cmn(0,9)= 0.2755732E-07;Cmn(0,10)= 0.2277464E-08;Cmn(0,11)= 0.1739730E-09;Cmn(0,12)= 0.1235311E-10;Cmn(0,13)=0.8193390E-12;Cmn(0,14)=0.5098109E-13
	Cmn(1,1)=-0.1562500E-01;Cmn(1,2)=-0.2222222E-02;Cmn(1,3)=-0.2893519E-03;Cmn(1,4)=-0.3401361E-04;Cmn(1,5)=-0.3616898E-05;Cmn(1,6)=-0.3499342E-06;Cmn(1,7)=-0.3100198E-07;Cmn(1,8)=-0.2530516E-08;Cmn(1,9)=-0.1913703E-09;Cmn(1,10)=-0.1347612E-10;Cmn(1,11)=-0.8876172E-12;Cmn(1,12)=-0.5490271E-13;Cmn(0,13)=0.0000000E+00;Cmn(0,14)=0.0000000E+00
	Cmn(2,1)= 0.4340278E-03;Cmn(2,2)= 0.4535147E-04;Cmn(2,3)= 0.4521123E-05;Cmn(2,4)= 0.4199211E-06;Cmn(2,5)= 0.3616898E-07;Cmn(2,6)= 0.2892018E-08;Cmn(2,7)= 0.2152916E-09;Cmn(2,8)= 0.1497347E-10;Cmn(2,9)= 0.9763789E-12;Cmn(2,10)= 0.5989387E-13;Cmn(2,11)= 0.0000000E+00;Cmn(2,12)= 0.0000000E+00;Cmn(0,13)=0.0000000E+00;Cmn(0,14)=0.0000000E+00
	Cmn(3,1)=-0.6781684E-05;Cmn(3,2)=-0.5598947E-06;Cmn(3,3)=-0.4521123E-07;Cmn(3,4)=-0.3470422E-08;Cmn(3,5)=-0.2511735E-09;Cmn(3,6)=-0.1711253E-10;Cmn(3,7)=-0.1098426E-11;Cmn(3,8)=-0.6654874E-13;Cmn(3,9)= 0.0000000E+00;Cmn(3,10)= 0.0000000E+00;Cmn(3,11)= 0.0000000E+00;Cmn(3,12)= 0.0000000E+00;Cmn(0,13)=0.0000000E+00;Cmn(0,14)=0.0000000E+00
	Cmn(4,1)= 0.6781684E-07;Cmn(4,2)= 0.4627229E-08;Cmn(4,3)= 0.3139669E-09;Cmn(4,4)= 0.2053504E-10;Cmn(4,5)= 0.1281497E-11;Cmn(4,6)= 0.7605571E-13;Cmn(4,7)= 0.0000000E+00;Cmn(4,8)= 0.0000000E+00;Cmn(4,9)= 0.0000000E+00;Cmn(4,10)= 0.0000000E+00;Cmn(4,11)= 0.0000000E+00;Cmn(4,12)= 0.0000000E+00;Cmn(0,13)=0.0000000E+00;Cmn(0,14)=0.0000000E+00
	Cmn(5,1)=-0.4709503E-09;Cmn(5,2)=-0.2738005E-10;Cmn(5,3)=-0.1601872E-11;Cmn(5,4)=-0.9126685E-13;Cmn(5,5)= 0.0000000E+00;Cmn(5,6)= 0.0000000E+00;Cmn(5,7)= 0.0000000E+00;Cmn(5,8)= 0.0000000E+00;Cmn(5,9)= 0.0000000E+00;Cmn(5,10)= 0.0000000E+00;Cmn(5,11)= 0.0000000E+00;Cmn(5,12)= 0.0000000E+00;Cmn(0,13)=0.0000000E+00;Cmn(0,14)=0.0000000E+00
	Cmn(6,1)= 0.2402808E-11;Cmn(6,2)= 0.1216891E-12;Cmn(6,3)= 0.0000000E+00;Cmn(6,4)= 0.0000000E+00;Cmn(6,5)= 0.0000000E+00;Cmn(6,6)= 0.0000000E+00;Cmn(6,7)= 0.0000000E+00;Cmn(6,8)= 0.0000000E+00;Cmn(6,9)= 0.0000000E+00;Cmn(6,10)= 0.0000000E+00;Cmn(6,11)= 0.0000000E+00;Cmn(6,12)= 0.0000000E+00;Cmn(0,13)=0.0000000E+00;Cmn(0,14)=0.0000000E+00
	infiFmid1=1.0;infiFmid2=1.0;infiFmid3=1.0;infiFmid=0.0
	do i=1,14
		infiFmid1=infiFmid1*y !infiFmid1表示y次幂
		infiFmid(0)=infiFmid(0)+infiFmid1*Cmn(0,i)
	enddo
	do i=1,6				
		infiFmid2=infiFmid2*x*x
		infiFmid1=1.0
		do j=1,14-2*i
			infiFmid1=infiFmid1*y
			infiFmid3=infiFmid1*infiFmid2*Cmn(i,j)
			infiFmid(0)=infiFmid(0)+infiFmid3
			infiFmid(1)=infiFmid(1)+infiFmid3*2*i
		enddo
	enddo
	infiFmid1=struveH0(x);infiFmid2=struveH1(x);infiFmid3=dbsj0(x)
	infiF(0)=exp(-y)*(-2*infiFmid3*log(y/x+R1/x)-PI*dbsy0(x)-PI*R1*infiFmid1/x-2*R1*infiFmid(0))
	infiF(1)=(-2*exp(-y))*(-dbsj1(x)*log(y/x+R1/x)-PI*dbsy1(x)/2-y*infiFmid3/(R1*x)+R1/x-PI*y*y*infiFmid1/(2*x*x*R1)-PI*infiFmid2*R1/(2*x)+x*infiFmid(0)/R1+R1*infiFmid(1)/x)
	infiF(2)=-2.0/R1-infiF(0)
elseif((2.0*x)<y)then
	infiF=0.0;infiFmid1=1.0;infiFmid2=-exp(-y)*dei(y);infiFmid=1.0
	do i=1,9
		do j=2*i-1,2*i
			if(j==1)THEN
				infiFmid1=1/y
			else
				infiFmid1=infiFmid1*(j-1)/y
			endif
			infiFmid2=infiFmid2+infiFmid1
		enddo
		infiFmid(0)=infiFmid(0)*(-1)*x*x/(4*i*i)
           if(i==1)THEN
				infiFmid(1)=1.0
			else
				infiFmid(1)=infiFmid(1)*(-1)*x*x/(4*(i-1)*i)
		endif
		infiF(0)=infiF(0)+infiFmid(0)*infiFmid2
		infiF(1)=infiF(1)+infiFmid(1)*infiFmid2
	enddo
	infiFmid2=infiFmid2+infiFmid1*18/y+infiFmid1*18*19/(y*y)
	infiFmid(1)=infiFmid(1)*(-1)*x*x/(4*9*10)
	infiF(1)=infiF(1)+infiFmid(1)*infiFmid2
	infiF(0)=infiF(0)*2-2*exp(-y)*dei(y)
	infiF(1)=-infiF(1)*x
	infiF(2)=-2.0/R1-infiF(0)
elseif((x>(2.0*y)).and.(x>=3.7))then
	infiF=0.0;infiFmid1=1.0;infiFmid2=1.0	
	Im(0)=1-exp(-y);Im(1)=y**2-2*y+2*Im(0);Im(2)=y**4-4*y**3+12*Im(1);Im(3)=y**6-6*y**5+30*Im(2);Im(4)=y**8-8*y**7+56*Im(3);Im(5)=y**10-10*y**9+90*Im(4)
	infiF(0)=Im(0)
	infiF(1)=Im(0)
	do i=1,5
		infiFmid1=-(2*i-1)/(2*i*x*x)*infiFmid1
		infiFmid2=-(2*i+1)/(2*i*x*x)*infiFmid2
		inFif(0)=inFif(0)+infiFmid1*Im(i)
		inFif(1)=inFif(1)+infiFmid2*Im(i)
	enddo
	inFif(0)=-PI*exp(-y)*(struveH0(x)+DBSY0(x))-2*inFif(0)/x
	inFif(1)=-2*exp(-y)+PI*exp(-y)*(struveH1(x)+DBSY1(x))+2*inFif(1)/(x*x)
	infiF(2)=-2.0/R1-infiF(0)
else
	call Gcfifpoint(coef0,coef1,x,y,ax,ay)
	x1=2.0*(x-ax(1))/(ax(2)-ax(1))-1.0
	y1=2.0*(y-ay(1))/(ay(2)-ay(1))-1.0
	infiF=0.0
	do i1=0,5
		m0x(i1)=Txy(i1,x1)
		m0y(i1)=Txy(i1,y1)
	enddo
	do i1=0,5
		do i2=0,5
			if(abs(coef0(i1,i2))>=1.0e-7)infiF(0)=infiF(0)+m0x(i1)*m0y(i2)*coef0(i1,i2)
			if(abs(coef1(i1,i2))>=1.0e-7)infiF(1)=infiF(1)+m0x(i1)*m0y(i2)*coef1(i1,i2)
		enddo
	enddo
	infiF(0)=-PI*exp(-y)*(struveH0(x)+DBSY0(x))-2*infiF(0)
	!求函数F(x,y)对x导数
	infiF(1)=-2*exp(-y)+PI*exp(-y)*(struveH1(x)+DBSY1(x))+2*infiF(1)
	!求函数F(x,y)对y导数
	infiF(2)=-2.0/R1-infiF(0)

endif
contains
real(8) Function Txy(n,x)
real(8)::x
integer(4)::n
if(n==0)then
	Txy=1
elseif(n==1)then
	Txy=x
elseif(n==2)then
	Txy=2*x**2-1
elseif(n==3)then
	Txy=4*x**3-3*x
elseif(n==4)then
	Txy=8*x**4-8*x**2+1
else
	Txy=16*x**5-20*x**3+5*x
endif
end function
end subroutine infiFxy
!---------------------------------------------------------



!****************************************************************************
!本程序计算频域有限水深Green函数中无限水深部分
!输出参数:infiF(0:1)分别为函数本身及对x的导数
!****************************************************************************
subroutine fiFxy(x,y,infiF)
implicit none

real(8),intent(in)::x,y
real(8),dimension(0:1),intent(out)::infiF
real(8),parameter::esp=1.0e-6  !for VF alter from 1.0e-12
real(8)::cta,R1,sta,cta2
real(8),allocatable,dimension(:)::xx,yy
real(8),allocatable,dimension(:,:)::Fxy,dxFxy
real(8)::Fx0y1,Fx0y2,Fx0y0,slope,dxFx0y1,dxFx0y2,dxFx0y0,zero,Njc,infiFmid1,infiFmid2,infiFmid3,infiFmid(0:1),Cmn(0:6,14),Im(0:5),Im1(0:5),Im2(0:5)
integer(4)::k,i,j,countx,county,nx,ny,i1,i2
real(8)::coef0(0:5,0:5),coef1(0:5,0:5),m0x(0:5),m0y(0:5),ax(2),ay(2),x1,y1

! PI=4.0*atan(1.0)
infiF=0.0
zero=0.0    !实型数虚实结合时要用变量，以保持类型一致。(release)
R1=sqrt(x*x+y*y)
if(x>8.and.y>18)then
	cta=y/R1;sta=x/R1;cta2=cta**2
	infiF(0)=(((((9+cta2*(105*cta2-90))/R1+cta*(15*cta2-9))/R1+3*cta2-1.0)/R1+cta)/R1+1.0)/R1
	infiF(1)=((((((945*cta2-630)*cta2+45)/R1+(105*cta2-45)*cta)/R1+15*cta2-3)/R1+3*cta)/R1+1)*sta/(R1**2)
	infiF(0)=-PI*exp(-y)*(struveH0(x)+DBSY0(x))-2*infiF(0)
	infiF(1)=2*PI*exp(-y)*DBSY1(x)+2*infiF(1)
elseif(abs(x)<1e-5.and.y/=0.0)then
	infiF(0)=-2*exp(-y)*dei(y)
	!求函数F(x,y)对x导数
	infiF(1)=0.0
	!求函数F(x,y)对y导数
elseif(abs(x)<1e-5.and.abs(y)<1e-5)then
	infiF=0.0
elseif((x<3.7).and.(y<2))then
	infiF=0.0
	Cmn(0,1)=+0.25000E+00;Cmn(0,2)=+0.55556E-01;Cmn(0,3)=+0.10417E-01;Cmn(0,4)=+0.16667E-02;Cmn(0,5)=+0.23148E-03;Cmn(0,6)=+0.28345E-04;Cmn(0,7)=+0.31002E-05;Cmn(0,8)=+0.30619E-06;Cmn(0,9)=+0.27557E-07;Cmn(0,10)=+0.22775E-08;Cmn(0,11)=+0.17397E-09;Cmn(0,12)=+0.12353E-10;Cmn(0,13)=0.81934E-12;Cmn(0,14)=0.50981E-13
	Cmn(1,1)=-0.15625E-01;Cmn(1,2)=-0.22222E-02;Cmn(1,3)=-0.28935E-03;Cmn(1,4)=-0.34014E-04;Cmn(1,5)=-0.36169E-05;Cmn(1,6)=-0.34993E-06;Cmn(1,7)=-0.31002E-07;Cmn(1,8)=-0.25305E-08;Cmn(1,9)=-0.19137E-09;Cmn(1,10)=-0.13476E-10;Cmn(1,11)=-0.88762E-12;Cmn(1,12)=-0.54903E-13;Cmn(0,13)=0.00000E+00;Cmn(0,14)=0.00000E+00
	Cmn(2,1)=+0.43403E-03;Cmn(2,2)=+0.45351E-04;Cmn(2,3)=+0.45211E-05;Cmn(2,4)=+0.41992E-06;Cmn(2,5)=+0.36169E-07;Cmn(2,6)=+0.28920E-08;Cmn(2,7)=+0.21529E-09;Cmn(2,8)=+0.14973E-10;Cmn(2,9)=+0.97638E-12;Cmn(2,10)=+0.59894E-13;Cmn(2,11)=+0.00000E+00;Cmn(2,12)=+0.00000E+00;Cmn(0,13)=0.00000E+00;Cmn(0,14)=0.00000E+00
	Cmn(3,1)=-0.67817E-05;Cmn(3,2)=-0.55989E-06;Cmn(3,3)=-0.45211E-07;Cmn(3,4)=-0.34704E-08;Cmn(3,5)=-0.25117E-09;Cmn(3,6)=-0.17113E-10;Cmn(3,7)=-0.10984E-11;Cmn(3,8)=-0.66549E-13;Cmn(3,9)=+0.00000E+00;Cmn(3,10)=+0.00000E+00;Cmn(3,11)=+0.00000E+00;Cmn(3,12)=+0.00000E+00;Cmn(0,13)=0.00000E+00;Cmn(0,14)=0.00000E+00
	Cmn(4,1)=+0.67817E-07;Cmn(4,2)=+0.46272E-08;Cmn(4,3)=+0.31397E-09;Cmn(4,4)=+0.20535E-10;Cmn(4,5)=+0.12815E-11;Cmn(4,6)=+0.76056E-13;Cmn(4,7)=+0.00000E+00;Cmn(4,8)=+0.00000E+00;Cmn(4,9)=+0.00000E+00;Cmn(4,10)=+0.00000E+00;Cmn(4,11)=+0.00000E+00;Cmn(4,12)=+0.00000E+00;Cmn(0,13)=0.00000E+00;Cmn(0,14)=0.00000E+00
	Cmn(5,1)=-0.47095E-09;Cmn(5,2)=-0.27380E-10;Cmn(5,3)=-0.16019E-11;Cmn(5,4)=-0.91267E-13;Cmn(5,5)=+0.00000E+00;Cmn(5,6)=+0.00000E+00;Cmn(5,7)=+0.00000E+00;Cmn(5,8)=+0.00000E+00;Cmn(5,9)=+0.00000E+00;Cmn(5,10)=+0.00000E+00;Cmn(5,11)=+0.00000E+00;Cmn(5,12)=+0.00000E+00;Cmn(0,13)=0.00000E+00;Cmn(0,14)=0.00000E+00
	Cmn(6,1)=+0.24028E-11;Cmn(6,2)=+0.12169E-12;Cmn(6,3)=+0.00000E+00;Cmn(6,4)=+0.00000E+00;Cmn(6,5)=+0.00000E+00;Cmn(6,6)=+0.00000E+00;Cmn(6,7)=+0.00000E+00;Cmn(6,8)=+0.00000E+00;Cmn(6,9)=+0.00000E+00;Cmn(6,10)=+0.00000E+00;Cmn(6,11)=+0.00000E+00;Cmn(6,12)=+0.00000E+00;Cmn(0,13)=0.00000E+00;Cmn(0,14)=0.00000E+00
	infiFmid1=1.0;infiFmid2=1.0;infiFmid3=1.0;infiFmid=0.0
	do i=1,14
		infiFmid1=infiFmid1*y !infiFmid1表示y次幂
		infiFmid(0)=infiFmid(0)+infiFmid1*Cmn(0,i)
	enddo
	do i=1,6				
		infiFmid2=infiFmid2*x*x
		infiFmid1=1.0
		do j=1,14-2*i
			infiFmid1=infiFmid1*y
			infiFmid3=infiFmid1*infiFmid2*Cmn(i,j)
			infiFmid(0)=infiFmid(0)+infiFmid3
			infiFmid(1)=infiFmid(1)+infiFmid3*2*i
		enddo
	enddo
	infiFmid1=struveH0(x);infiFmid2=struveH1(x);infiFmid3=dbsj0(x)
	infiF(0)=exp(-y)*(-2*infiFmid3*log(y/x+R1/x)-PI*dbsy0(x)-PI*R1*infiFmid1/x-2*R1*infiFmid(0))
	infiF(1)=(-2*exp(-y))*(-dbsj1(x)*log(y/x+R1/x)-PI*dbsy1(x)/2-y*infiFmid3/(R1*x)+R1/x-PI*y*y*infiFmid1/(2*x*x*R1)-PI*infiFmid2*R1/(2*x)+x*infiFmid(0)/R1+R1*infiFmid(1)/x)
elseif((2.0*x)<y)then
	infiF=0.0;infiFmid1=1.0;infiFmid2=-exp(-y)*dei(y);infiFmid=1.0
	do i=1,9
		do j=2*i-1,2*i
			if(j==1)THEN
				infiFmid1=1/y
			else
				infiFmid1=infiFmid1*(j-1)/y
			endif
			infiFmid2=infiFmid2+infiFmid1
		enddo
		infiFmid(0)=infiFmid(0)*(-1)*x*x/(4*i*i)
           if(i==1)THEN
				infiFmid(1)=1.0
			else
				infiFmid(1)=infiFmid(1)*(-1)*x*x/(4*(i-1)*i)
		endif
		infiF(0)=infiF(0)+infiFmid(0)*infiFmid2
		infiF(1)=infiF(1)+infiFmid(1)*infiFmid2
	enddo
	infiFmid2=infiFmid2+infiFmid1*18/y+infiFmid1*18*19/(y*y)
	infiFmid(1)=infiFmid(1)*(-1)*x*x/(4*9*10)
	infiF(1)=infiF(1)+infiFmid(1)*infiFmid2
	infiF(0)=infiF(0)*2-2*exp(-y)*dei(y)
	infiF(1)=-infiF(1)*x
elseif((x>(2.0*y)).and.(x>=3.7))then
	infiF=0.0;infiFmid1=1.0;infiFmid2=1.0	
	Im(0)=1-exp(-y);Im(1)=y**2-2*y+2*Im(0);Im(2)=y**4-4*y**3+12*Im(1);Im(3)=y**6-6*y**5+30*Im(2);Im(4)=y**8-8*y**7+56*Im(3);Im(5)=y**10-10*y**9+90*Im(4)
	infiF(0)=Im(0)
	infiF(1)=Im(0)
	do i=1,5
		infiFmid1=-(2*i-1)/(2*i*x*x)*infiFmid1
		infiFmid2=-(2*i+1)/(2*i*x*x)*infiFmid2
		inFif(0)=inFif(0)+infiFmid1*Im(i)
		inFif(1)=inFif(1)+infiFmid2*Im(i)
	enddo
	inFif(0)=-PI*exp(-y)*(struveH0(x)+DBSY0(x))-2*inFif(0)/x
	inFif(1)=-2*exp(-y)+PI*exp(-y)*(struveH1(x)+DBSY1(x))+2*inFif(1)/(x*x)
else
	call Gcfifpoint(coef0,coef1,x,y,ax,ay)
	x1=2.0*(x-ax(1))/(ax(2)-ax(1))-1.0
	y1=2.0*(y-ay(1))/(ay(2)-ay(1))-1.0
	infiF=0.0
	do i1=0,5
		m0x(i1)=Txy(i1,x1)
		m0y(i1)=Txy(i1,y1)
	enddo
	do i1=0,5
		do i2=0,5
			if(abs(coef0(i1,i2))>=1.0e-7)infiF(0)=infiF(0)+m0x(i1)*m0y(i2)*coef0(i1,i2)
			if(abs(coef1(i1,i2))>=1.0e-7)infiF(1)=infiF(1)+m0x(i1)*m0y(i2)*coef1(i1,i2)
		enddo
	enddo
	infiF(0)=-PI*exp(-y)*(struveH0(x)+DBSY0(x))-2*infiF(0)
	!求函数F(x,y)对x导数
	infiF(1)=-2*exp(-y)+PI*exp(-y)*(struveH1(x)+DBSY1(x))+2*infiF(1)
	!求函数F(x,y)对y导数
endif
contains
real(8) Function Txy(n,x)
real(8)::x
integer(4)::n
if(n==0)then
	Txy=1
elseif(n==1)then
	Txy=x
elseif(n==2)then
	Txy=2*x**2-1
elseif(n==3)then
	Txy=4*x**3-3*x
elseif(n==4)then
	Txy=8*x**4-8*x**2+1
else
	Txy=16*x**5-20*x**3+5*x
endif
end function
end subroutine fiFxy
!--------------------------------------------------------------------------



!*******************************************************************************
!用渐近表达式计算零阶和一阶struve函数:H0(x),H1(x)      
!*******************************************************************************
function struveH0(x) result(x_result)

implicit none
real(8)::x,x_result
real(8)::a0,a1,a2,a3,b1,b2,b3,tem1,tem2
real(8),dimension(1:6)::c
! PI=4.0*atan(1.0)
c=(/1.909859164,-1.909855001,0.687514637,-0.126164557,0.013828813,-0.000876918/)
a0=0.99999906;a1=4.77228920;a2=3.85542044;a3=0.32303607
b1=4.88331068;b2=4.28957333;b3=0.52120508
if(x<=3.0.and.x>=0)then
	x_result=(x/3)*(c(1)+(x/3)**2*(c(2)+(x/3)**2*(c(3)+(x/3)**2*c(4))))
else
	tem1=a0+(3/x)**2*(a1+(3/x)**2*(a2+(3/x)**2*a3))
	tem2=1+(3/x)**2*(b1+(3/x)**2*(b2+(3/x)**2*b3))
	x_result=DBSY0(x)+2*tem1/(PI*x*tem2) 
endif
end function
!------------------------------------------------------------------------
function struveH1(x) result(x_result)

implicit none
real(8)::x,x_result
real(8)::a0,a1,a2,a3,b1,b2,b3,tem1,tem2
real(8),dimension(1:6)::c
! PI=4.0*atan(1.0)
c=(/1.909859286,-1.145914713,0.294656958,-0.042070508,0.003785727,-0.000207183/)
a0=1.00000004;a1=3.92205313;a2=2.64893033;a3=0.27450895
b1=3.81095112;b2=2.26216956;b3=0.10885141
if(x<=3.0.and.x>=0)then
	x_result=(x/3)**2*(c(1)+(x/3)**2*(c(2)+(x/3)**2*(c(3)+(x/3)**2*(c(4)+(x/3)**2*c(5)))))
else
	tem1=a0+(3/x)**2*(a1+(3/x)**2*(a2+(3/x)**2*a3))
	tem2=1+(3/x)**2*(b1+(3/x)**2*(b2+(3/x)**2*b3))
	x_result=DBSY1(x)+2*tem1/(PI*tem2) 
endif
end function
!************************************************************************************

subroutine Gcfifpoint(a,b,x,y,ax,ay)
real(8)::a(0:5,0:5),b(0:5,0:5),x,y,ax(2),ay(2)
integer(4)::m,n
if((x>=1.0).and.(x<=2.0))then
	ax(1)=1.0;ax(2)=2.0
	if(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		a(0,0)= 0.410769;a(0,1)=-0.033887;a(0,2)=-0.002093;a(0,3)= 0.000438;a(0,4)=-0.000037;a(0,5)= 0.000002;a(1,0)=-0.070169;a(1,1)= 0.019589;a(1,2)=-0.000781;a(1,3)=-0.000084;a(1,4)= 0.000017;a(1,5)=-0.000002;a(2,0)= 0.004889;a(2,1)=-0.002619;a(2,2)= 0.000280;a(2,3)=-0.000010;a(2,4)=-0.000001;a(2,5)= 0.000000;a(3,0)=-0.000316;a(3,1)= 0.000242;a(3,2)=-0.000042;a(3,3)= 0.000004;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000028;a(4,1)=-0.000019;a(4,2)= 0.000004;a(4,3)=-0.000001;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000004;a(5,1)= 0.000002;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
		b(0,0)= 0.142271;b(0,1)=-0.040651;b(0,2)= 0.001814;b(0,3)= 0.000142;b(0,4)=-0.000032;b(0,5)= 0.000003;b(1,0)=-0.039565;b(1,1)= 0.021259;b(1,2)=-0.002301;b(1,3)= 0.000093;b(1,4)= 0.000009;b(1,5)=-0.000002;b(2,0)= 0.003866;b(2,1)=-0.002946;b(2,2)= 0.000504;b(2,3)=-0.000051;b(2,4)= 0.000003;b(2,5)= 0.000000;b(3,0)=-0.000450;b(3,1)= 0.000309;b(3,2)=-0.000064;b(3,3)= 0.000010;b(3,4)=-0.000001;b(3,5)= 0.000000;b(4,0)= 0.000083;b(4,1)=-0.000039;b(4,2)= 0.000006;b(4,3)=-0.000001;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)=-0.000016;b(5,1)= 0.000007;b(5,2)=-0.000001;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		a(0,0)= 0.337481;a(0,1)=-0.037130;a(0,2)= 0.000681;a(0,3)= 0.000089;a(0,4)=-0.000011;a(0,5)= 0.000001;a(1,0)=-0.038377;a(1,1)= 0.012129;a(1,2)=-0.000898;a(1,3)= 0.000028;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.001418;a(2,1)=-0.000978;a(2,2)= 0.000132;a(2,3)=-0.000010;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000045;a(3,1)= 0.000053;a(3,2)=-0.000011;a(3,3)= 0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000007;a(4,1)=-0.000004;a(4,2)= 0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000002;a(5,1)= 0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000	
		b(0,0)= 0.077041;b(0,1)=-0.024584;b(0,2)= 0.001861;b(0,3)=-0.000064;b(0,4)=-0.000002;b(0,5)= 0.000000;b(1,0)=-0.011459;b(1,1)= 0.007890;b(1,2)=-0.001071;b(1,3)= 0.000084;b(1,4)=-0.000004;b(1,5)= 0.000000;b(2,0)= 0.000576;b(2,1)=-0.000653;b(2,2)= 0.000131;b(2,3)=-0.000016;b(2,4)= 0.000001;b(2,5)= 0.000000;b(3,0)=-0.000118;b(3,1)= 0.000065;b(3,2)=-0.000011;b(3,3)= 0.000001;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000033;b(4,1)=-0.000014;b(4,2)= 0.000002;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)=-0.000006;b(5,1)= 0.000003;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
	elseif(y<=5.0)then
		ay(1)=4.0;ay(2)=5.0
		a(0,0)= 0.270513;a(0,1)=-0.029539;a(0,2)= 0.001056;a(0,3)=-0.000005;a(0,4)=-0.000002;a(0,5)= 0.000000;a(1,0)=-0.020143;a(1,1)= 0.006410;a(1,2)=-0.000534;a(1,3)= 0.000027;a(1,4)=-0.000001;a(1,5)= 0.000000;a(2,0)= 0.000211;a(2,1)=-0.000299;a(2,2)= 0.000048;a(2,3)=-0.000004;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000009;a(3,1)= 0.000008;a(3,2)=-0.000002;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
		b(0,0)= 0.040242;b(0,1)=-0.012868;b(0,2)= 0.001083;b(0,3)=-0.000057;b(0,4)= 0.000002;b(0,5)= 0.000000;b(1,0)=-0.001726;b(1,1)= 0.002411;b(1,2)=-0.000385;b(1,3)= 0.000035;b(1,4)=-0.000002;b(1,5)= 0.000000;b(2,0)=-0.000090;b(2,1)=-0.000096;b(2,2)= 0.000029;b(2,3)=-0.000004;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)=-0.000041;b(3,1)= 0.000019;b(3,2)=-0.000003;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000013;b(4,1)=-0.000006;b(4,2)= 0.000001;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)=-0.000002;b(5,1)= 0.000001;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
elseif(x<=3.0)then
	ax(1)=2.0;ax(2)=3.0
	if(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		a(0,0)= 0.300726;a(0,1)=-0.009009;a(0,2)=-0.002448;a(0,3)= 0.000269;a(0,4)=-0.000015;a(0,5)= 0.000000;a(1,0)=-0.041837;a(1,1)= 0.006722;a(1,2)= 0.000223;a(1,3)=-0.000073;a(1,4)= 0.000006;a(1,5)= 0.000000;a(2,0)= 0.002502;a(2,1)=-0.000876;a(2,2)= 0.000030;a(2,3)= 0.000006;a(2,4)=-0.000001;a(2,5)= 0.000000;a(3,0)=-0.000128;a(3,1)= 0.000079;a(3,2)=-0.000007;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000006;a(4,1)=-0.000005;a(4,2)= 0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
		b(0,0)= 0.084445;b(0,1)=-0.013920;b(0,2)=-0.000400;b(0,3)= 0.000146;b(0,4)=-0.000012;b(0,5)= 0.000000;b(1,0)=-0.020105;b(1,1)= 0.007095;b(1,2)=-0.000251;b(1,3)=-0.000047;b(1,4)= 0.000007;b(1,5)= 0.000000;b(2,0)= 0.001544;b(2,1)=-0.000951;b(2,2)= 0.000092;b(2,3)= 0.000001;b(2,4)=-0.000001;b(2,5)= 0.000000;b(3,0)=-0.000091;b(3,1)= 0.000087;b(3,2)=-0.000015;b(3,3)= 0.000001;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000005;b(4,1)=-0.000006;b(4,2)= 0.000002;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.057073;b(0,1)=-0.012731;b(0,2)= 0.000481;b(0,3)= 0.000022;b(0,4)=-0.000004;b(0,5)= 0.000000;b(1,0)=-0.008710;b(1,1)= 0.004209;b(1,2)=-0.000372;b(1,3)= 0.000011;b(1,4)= 0.000001;b(1,5)= 0.000000;b(2,0)= 0.000287;b(2,1)=-0.000339;b(2,2)= 0.000054;b(2,3)=-0.000004;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)= 0.000003;b(3,1)= 0.000015;b(3,2)=-0.000004;b(3,3)= 0.000001;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000000;b(4,1)= 0.000000;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.270905;a(0,1)=-0.019095;a(0,2)=-0.000372;a(0,3)= 0.000097;a(0,4)=-0.000007;a(0,5)= 0.000000;a(1,0)=-0.028465;a(1,1)= 0.006281;a(1,2)=-0.000227;a(1,3)=-0.000012;a(1,4)= 0.000002;a(1,5)= 0.000000;a(2,0)= 0.001089;a(2,1)=-0.000524;a(2,2)= 0.000046;a(2,3)=-0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000024;a(3,1)= 0.000028;a(3,2)=-0.000004;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.033326;b(0,1)=-0.010653;b(0,2)= 0.000835;b(0,3)=-0.000029;b(0,4)=-0.000001;b(0,5)= 0.000000;b(1,0)=-0.002327;b(1,1)= 0.002170;b(1,2)=-0.000334;b(1,3)= 0.000028;b(1,4)=-0.000001;b(1,5)= 0.000000;b(2,0)=-0.000100;b(2,1)=-0.000068;b(2,2)= 0.000026;b(2,3)=-0.000004;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)= 0.000012;b(3,1)=-0.000003;b(3,2)=-0.000001;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000000;b(4,1)= 0.000000;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.226418;a(0,1)=-0.024570;a(0,2)= 0.000593;a(0,3)= 0.000045;a(0,4)=-0.000007;a(0,5)= 0.000001;a(1,0)=-0.016688;a(1,1)= 0.005309;a(1,2)=-0.000411;a(1,3)= 0.000013;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000292;a(2,1)=-0.000272;a(2,2)= 0.000042;a(2,3)=-0.000003;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000008;a(3,1)= 0.000006;a(3,2)=-0.000002;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)=-0.000001;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.013894;b(0,1)=-0.008040;b(0,2)= 0.001372;b(0,3)=-0.000161;b(0,4)= 0.000012;b(0,5)= 0.000000;b(1,0)= 0.000363;b(1,1)= 0.000520;b(1,2)=-0.000245;b(1,3)= 0.000054;b(1,4)=-0.000008;b(1,5)= 0.000001;b(2,0)=-0.000102;b(2,1)= 0.000048;b(2,2)= 0.000000;b(2,3)=-0.000002;b(2,4)= 0.000001;b(2,5)= 0.000000;b(3,0)= 0.000004;b(3,1)=-0.000004;b(3,2)= 0.000001;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000000;b(4,1)= 0.000000;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.167139;a(0,1)=-0.032678;a(0,2)= 0.002555;a(0,3)=-0.000116;a(0,4)=-0.000006;a(0,5)= 0.000002;a(1,0)=-0.006972;a(1,1)= 0.004032;a(1,2)=-0.000686;a(1,3)= 0.000080;a(1,4)=-0.000006;a(1,5)= 0.000000;a(2,0)=-0.000045;a(2,1)=-0.000066;a(2,2)= 0.000031;a(2,3)=-0.000007;a(2,4)= 0.000001;a(2,5)= 0.000000;a(3,0)= 0.000009;a(3,1)=-0.000004;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
elseif(x<=3.75)then
	ax(1)=3.0;ax(2)=3.75
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.030762;b(0,1)= 0.026531;b(0,2)=-0.003994;b(0,3)= 0.000227;b(0,4)=-0.000009;b(0,5)= 0.000001;b(1,0)=-0.006634;b(1,1)=-0.005601;b(1,2)= 0.000992;b(1,3)=-0.000041;b(1,4)= 0.000000;b(1,5)= 0.000000;b(2,0)= 0.000533;b(2,1)= 0.000437;b(2,2)=-0.000094;b(2,3)= 0.000003;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)=-0.000038;b(3,1)=-0.000030;b(3,2)= 0.000008;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000002;b(4,1)= 0.000002;b(4,2)=-0.000001;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.104722;a(0,1)= 0.091638;a(0,2)=-0.012146;a(0,3)= 0.000883;a(0,4)=-0.000052;a(0,5)= 0.000003;a(1,0)=-0.011436;a(1,1)=-0.009867;a(1,2)= 0.001480;a(1,3)=-0.000085;a(1,4)= 0.000004;a(1,5)= 0.000000;a(2,0)= 0.000618;a(2,1)= 0.000522;a(2,2)=-0.000092;a(2,3)= 0.000004;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000033;a(3,1)=-0.000027;a(3,2)= 0.000006;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)= 0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.059214;b(0,1)= 0.003695;b(0,2)=-0.001856;b(0,3)= 0.000142;b(0,4)=-0.000004;b(0,5)= 0.000000;b(1,0)=-0.011553;b(1,1)= 0.000233;b(1,2)= 0.000456;b(1,3)=-0.000045;b(1,4)= 0.000001;b(1,5)= 0.000000;b(2,0)= 0.000816;b(2,1)=-0.000105;b(2,2)=-0.000037;b(2,3)= 0.000005;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)=-0.000050;b(3,1)= 0.000013;b(3,2)= 0.000002;b(3,3)=-0.000001;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000003;b(4,1)=-0.000001;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.216680;a(0,1)= 0.026134;a(0,2)=-0.005108;a(0,3)= 0.000372;a(0,4)=-0.000018;a(0,5)= 0.000001;a(1,0)=-0.022052;a(1,1)=-0.001405;a(1,2)= 0.000689;a(1,3)=-0.000052;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.001078;a(2,1)=-0.000021;a(2,2)=-0.000043;a(2,3)= 0.000004;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000051;a(3,1)= 0.000007;a(3,2)= 0.000002;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.056386;b(0,1)=-0.005423;b(0,2)=-0.000538;b(0,3)= 0.000077;b(0,4)=-0.000004;b(0,5)= 0.000000;b(1,0)=-0.008918;b(1,1)= 0.002057;b(1,2)= 0.000043;b(1,3)=-0.000022;b(1,4)= 0.000002;b(1,5)= 0.000000;b(2,0)= 0.000475;b(2,1)=-0.000200;b(2,2)= 0.000006;b(2,3)= 0.000002;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)=-0.000020;b(3,1)= 0.000014;b(3,2)=-0.000001;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000001;b(4,1)=-0.000001;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.239347;a(0,1)=-0.000888;a(0,2)=-0.001995;a(0,3)= 0.000171;a(0,4)=-0.000008;a(0,5)= 0.000000;a(1,0)=-0.021056;a(1,1)= 0.001996;a(1,2)= 0.000203;a(1,3)=-0.000029;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000834;a(2,1)=-0.000192;a(2,2)=-0.000004;a(2,3)= 0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000030;a(3,1)= 0.000012;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000001;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.043480;b(0,1)=-0.006990;b(0,2)= 0.000057;b(0,3)= 0.000026;b(0,4)=-0.000002;b(0,5)= 0.000000;b(1,0)=-0.005014;b(1,1)= 0.001739;b(1,2)=-0.000089;b(1,3)=-0.000003;b(1,4)= 0.000001;b(1,5)= 0.000000;b(2,0)= 0.000158;b(2,1)=-0.000113;b(2,2)= 0.000012;b(2,3)= 0.000000;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)=-0.000002;b(3,1)= 0.000005;b(3,2)=-0.000001;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000000;b(4,1)= 0.000000;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.226750;a(0,1)=-0.010541;a(0,2)=-0.000586;a(0,3)= 0.000074;a(0,4)=-0.000004;a(0,5)= 0.000000;a(1,0)=-0.016275;a(1,1)= 0.002600;a(1,2)=-0.000019;a(1,3)=-0.000010;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000470;a(2,1)=-0.000163;a(2,2)= 0.000008;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000010;a(3,1)= 0.000007;a(3,2)=-0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.028915;b(0,1)=-0.007281;b(0,2)= 0.000381;b(0,3)= 0.000003;b(0,4)=-0.000002;b(0,5)= 0.000000;b(1,0)=-0.001997;b(1,1)= 0.001227;b(1,2)=-0.000136;b(1,3)= 0.000006;b(1,4)= 0.000000;b(1,5)= 0.000000;b(2,0)= 0.000000;b(2,1)=-0.000045;b(2,2)= 0.000010;b(2,3)=-0.000001;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)= 0.000003;b(3,1)= 0.000000;b(3,2)= 0.000000;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000000;b(4,1)= 0.000000;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.198999;a(0,1)=-0.016687;a(0,2)= 0.000068;a(0,3)= 0.000055;a(0,4)=-0.000005;a(0,5)= 0.000000;a(1,0)=-0.010843;a(1,1)= 0.002722;a(1,2)=-0.000141;a(1,3)=-0.000001;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000187;a(2,1)=-0.000115;a(2,2)= 0.000013;a(2,3)=-0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000000;a(3,1)= 0.000003;a(3,2)=-0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.014035;b(0,1)=-0.006938;b(0,2)= 0.000960;b(0,3)=-0.000081;b(0,4)= 0.000002;b(0,5)= 0.000000;b(1,0)=-0.000160;b(1,1)= 0.000545;b(1,2)=-0.000164;b(1,3)= 0.000027;b(1,4)=-0.000003;b(1,5)= 0.000000;b(2,0)=-0.000035;b(2,1)= 0.000007;b(2,2)= 0.000004;b(2,3)=-0.000001;b(2,4)= 0.000000;b(2,5)= 0.000000;b(3,0)= 0.000002;b(3,1)=-0.000001;b(3,2)= 0.000000;b(3,3)= 0.000000;b(3,4)= 0.000000;b(3,5)= 0.000000;b(4,0)= 0.000000;b(4,1)= 0.000000;b(4,2)= 0.000000;b(4,3)= 0.000000;b(4,4)= 0.000000;b(4,5)= 0.000000;b(5,0)= 0.000000;b(5,1)= 0.000000;b(5,2)= 0.000000;b(5,3)= 0.000000;b(5,4)= 0.000000;b(5,5)= 0.000000
		a(0,0)= 0.154848;a(0,1)=-0.026061;a(0,2)= 0.001525;a(0,3)=-0.000010;a(0,4)=-0.000012;a(0,5)= 0.000002;a(1,0)=-0.005270;a(1,1)= 0.002603;a(1,2)=-0.000359;a(1,3)= 0.000030;a(1,4)=-0.000001;a(1,5)= 0.000000;a(2,0)= 0.000015;a(2,1)=-0.000051;a(2,2)= 0.000015;a(2,3)=-0.000003;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000002;a(3,1)= 0.000000;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00420839;b(0,1)=-0.00295043;b(0,2)= 0.00067739;b(0,3)=-0.00012157;b(0,4)= 0.00001790;b(0,5)=-0.00000206;b(1,0)= 0.00021999;b(1,1)=-0.00005694;b(1,2)=-0.00001773;b(1,3)= 0.00001029;b(1,4)=-0.00000295;b(1,5)= 0.00000057;b(2,0)=-0.00001192;b(2,1)= 0.00001068;b(2,2)=-0.00000257;b(2,3)= 0.00000034;b(2,4)= 0.00000001;b(2,5)=-0.00000002;b(3,0)= 0.00000014;b(3,1)=-0.00000023;b(3,2)= 0.00000010;b(3,3)=-0.00000003;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000000;b(4,1)= 0.00000000;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.102817;a(0,1)=-0.024509;a(0,2)= 0.002689;a(0,3)=-0.000257;a(0,4)= 0.000018;a(0,5)= 0.000000;a(1,0)=-0.001580;a(1,1)= 0.001108;a(1,2)=-0.000255;a(1,3)= 0.000046;a(1,4)=-0.000007;a(1,5)= 0.000001;a(2,0)=-0.000021;a(2,1)= 0.000005;a(2,2)= 0.000002;a(2,3)=-0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000001;a(3,1)=-0.000001;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00095060;b(0,1)=-0.00067508;b(0,2)= 0.00016367;b(0,3)=-0.00003312;b(0,4)= 0.00000614;b(0,5)=-0.00000100;b(1,0)= 0.00008676;b(1,1)=-0.00005425;b(1,2)= 0.00001064;b(1,3)=-0.00000149;b(1,4)= 0.00000011;b(1,5)= 0.00000001;b(2,0)=-0.00000138;b(2,1)= 0.00000145;b(2,2)=-0.00000049;b(2,3)= 0.00000013;b(2,4)=-0.00000003;b(2,5)= 0.00000001;b(3,0)=-0.00000001;b(3,1)= 0.00000001;b(3,2)= 0.00000000;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000000;b(4,1)= 0.00000000;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.063219;a(0,1)=-0.015446;a(0,2)= 0.001849;a(0,3)=-0.000216;a(0,4)= 0.000025;a(0,5)=-0.000003;a(1,0)=-0.000357;a(1,1)= 0.000253;a(1,2)=-0.000061;a(1,3)= 0.000012;a(1,4)=-0.000002;a(1,5)= 0.000000;a(2,0)=-0.000008;a(2,1)= 0.000005;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000000;a(3,1)= 0.000000;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000;
	else
		write(*,*)'error'
	endif
elseif(x<=5.0)then
	ax(1)=3.75;ax(2)=5.0
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.01875593;b(0,1)= 0.01631009;b(0,2)=-0.00228898;b(0,3)= 0.00014826;b(0,4)=-0.00000810;b(0,5)= 0.00000050;b(1,0)=-0.00526106;b(1,1)=-0.00451652;b(1,2)= 0.00070631;b(1,3)=-0.00003669;b(1,4)= 0.00000140;b(1,5)=-0.00000014;b(2,0)= 0.00055199;b(2,1)= 0.00046565;b(2,2)=-0.00008306;b(2,3)= 0.00000327;b(2,4)= 0.00000001;b(2,5)= 0.00000002;b(3,0)=-0.00005123;b(3,1)=-0.00004226;b(3,2)= 0.00000874;b(3,3)=-0.00000025;b(3,4)=-0.00000002;b(3,5)= 0.00000000;b(4,0)= 0.00000446;b(4,1)= 0.00000358;b(4,2)=-0.00000087;b(4,3)= 0.00000002;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000037;b(5,1)=-0.00000028;b(5,2)= 0.00000008;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.081447;a(0,1)= 0.071468;a(0,2)=-0.009230;a(0,3)= 0.000704;a(0,4)=-0.000043;a(0,5)= 0.000002;a(1,0)=-0.011550;a(1,1)=-0.010048;a(1,2)= 0.001405;a(1,3)=-0.000092;a(1,4)= 0.000005;a(1,5)= 0.000000;a(2,0)= 0.000814;a(2,1)= 0.000699;a(2,2)=-0.000109;a(2,3)= 0.000006;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000057;a(3,1)=-0.000048;a(3,2)= 0.000009;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000004;a(4,1)= 0.000003;a(4,2)=-0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.03752991;b(0,1)= 0.00349098;b(0,2)=-0.00104646;b(0,3)= 0.00007196;b(0,4)=-0.00000271;b(0,5)= 0.00000011;b(1,0)=-0.00986289;b(1,1)=-0.00038555;b(1,2)= 0.00034444;b(1,3)=-0.00002574;b(1,4)= 0.00000048;b(1,5)= 0.00000001;b(2,0)= 0.00095058;b(2,1)=-0.00003119;b(2,2)=-0.00004012;b(2,3)= 0.00000375;b(2,4)=-0.00000003;b(2,5)=-0.00000001;b(3,0)=-0.00007955;b(3,1)= 0.00000980;b(3,2)= 0.00000383;b(3,3)=-0.00000049;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000613;b(4,1)=-0.00000144;b(4,2)=-0.00000031;b(4,3)= 0.00000006;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000044;b(5,1)= 0.00000016;b(5,2)= 0.00000002;b(5,3)=-0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.170865;a(0,1)= 0.022465;a(0,2)=-0.003760;a(0,3)= 0.000276;a(0,4)=-0.000015;a(0,5)= 0.000001;a(1,0)=-0.023159;a(1,1)=-0.002192;a(1,2)= 0.000641;a(1,3)=-0.000044;a(1,4)= 0.000002;a(1,5)= 0.000000;a(2,0)= 0.001529;a(2,1)= 0.000062;a(2,2)=-0.000053;a(2,3)= 0.000004;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000098;a(3,1)= 0.000003;a(3,2)= 0.000004;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000006;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.03843191;b(0,1)=-0.00204777;b(0,2)=-0.00039507;b(0,3)= 0.00003917;b(0,4)=-0.00000162;b(0,5)= 0.00000003;b(1,0)=-0.00874887;b(1,1)= 0.00128868;b(1,2)= 0.00009192;b(1,3)=-0.00001565;b(1,4)= 0.00000073;b(1,5)= 0.00000000;b(2,0)= 0.00069494;b(2,1)=-0.00019381;b(2,2)=-0.00000357;b(2,3)= 0.00000207;b(2,4)=-0.00000014;b(2,5)= 0.00000000;b(3,0)=-0.00004522;b(3,1)= 0.00002086;b(3,2)=-0.00000054;b(3,3)=-0.00000020;b(3,4)= 0.00000002;b(3,5)= 0.00000000;b(4,0)= 0.00000251;b(4,1)=-0.00000181;b(4,2)= 0.00000013;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000012;b(5,1)= 0.00000013;b(5,2)=-0.00000002;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.193859;a(0,1)= 0.002367;a(0,2)=-0.001536;a(0,3)= 0.000118;a(0,4)=-0.000006;a(0,5)= 0.000000;a(1,0)=-0.023803;a(1,1)= 0.001219;a(1,2)= 0.000246;a(1,3)=-0.000024;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.001360;a(2,1)=-0.000198;a(2,2)=-0.000014;a(2,3)= 0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000072;a(3,1)= 0.000020;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000004;a(4,1)=-0.000002;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.03238130;b(0,1)=-0.00372614;b(0,2)=-0.00006191;b(0,3)= 0.00001781;b(0,4)=-0.00000104;b(0,5)= 0.00000003;b(1,0)=-0.00589615;b(1,1)= 0.00146268;b(1,2)=-0.00003019;b(1,3)=-0.00000531;b(1,4)= 0.00000050;b(1,5)=-0.00000002;b(2,0)= 0.00033271;b(2,1)=-0.00015768;b(2,2)= 0.00000945;b(2,3)= 0.00000028;b(2,4)=-0.00000007;b(2,5)= 0.00000000;b(3,0)=-0.00001189;b(3,1)= 0.00001185;b(3,2)=-0.00000132;b(3,3)= 0.00000003;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000006;b(4,1)=-0.00000064;b(4,2)= 0.00000012;b(4,3)=-0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000003;b(5,1)= 0.00000002;b(5,2)=-0.00000001;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.189845;a(0,1)=-0.005574;a(0,2)=-0.000562;a(0,3)= 0.000053;a(0,4)=-0.000003;a(0,5)= 0.000000;a(1,0)=-0.020134;a(1,1)= 0.002280;a(1,2)= 0.000042;a(1,3)=-0.000011;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000919;a(2,1)=-0.000227;a(2,2)= 0.000005;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000035;a(3,1)= 0.000016;a(3,2)=-0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000001;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.02381855;b(0,1)=-0.00466379;b(0,2)= 0.00013937;b(0,3)= 0.00001016;b(0,4)=-0.00000143;b(0,5)= 0.00000008;b(1,0)=-0.00303314;b(1,1)= 0.00133476;b(1,2)=-0.00010017;b(1,3)= 0.00000092;b(1,4)= 0.00000046;b(1,5)=-0.00000004;b(2,0)= 0.00007386;b(2,1)=-0.00009634;b(2,2)= 0.00001375;b(2,3)=-0.00000076;b(2,4)=-0.00000002;b(2,5)= 0.00000001;b(3,0)= 0.00000357;b(3,1)= 0.00000367;b(3,2)=-0.00000110;b(3,3)= 0.00000012;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)=-0.00000050;b(4,1)= 0.00000003;b(4,2)= 0.00000005;b(4,3)=-0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000003;b(5,1)=-0.00000001;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.172998;a(0,1)=-0.010974;a(0,2)=-0.000161;a(0,3)= 0.000047;a(0,4)=-0.000003;a(0,5)= 0.000000;a(1,0)=-0.014864;a(1,1)= 0.002885;a(1,2)=-0.000083;a(1,3)=-0.000007;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000474;a(2,1)=-0.000208;a(2,2)= 0.000015;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000008;a(3,1)= 0.000010;a(3,2)=-0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.01320314;b(0,1)=-0.00548042;b(0,2)= 0.00059164;b(0,3)=-0.00003004;b(0,4)=-0.00000175;b(0,5)= 0.00000054;b(1,0)=-0.00068451;b(1,1)= 0.00089324;b(1,2)=-0.00019538;b(1,3)= 0.00002295;b(1,4)=-0.00000121;b(1,5)=-0.00000008;b(2,0)=-0.00004131;b(2,1)=-0.00001627;b(2,2)= 0.00001243;b(2,3)=-0.00000276;b(2,4)= 0.00000034;b(2,5)=-0.00000002;b(3,0)= 0.00000448;b(3,1)=-0.00000212;b(3,2)=-0.00000015;b(3,3)= 0.00000015;b(3,4)=-0.00000003;b(3,5)= 0.00000000;b(4,0)=-0.00000020;b(4,1)= 0.00000020;b(4,2)=-0.00000004;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)=-0.00000001;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.141227;a(0,1)=-0.019948;a(0,2)= 0.000783;a(0,3)= 0.000039;a(0,4)=-0.000012;a(0,5)= 0.000001;a(1,0)=-0.008265;a(1,1)= 0.003420;a(1,2)=-0.000366;a(1,3)= 0.000018;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000108;a(2,1)=-0.000140;a(2,2)= 0.000031;a(2,3)=-0.000004;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000004;a(3,1)= 0.000002;a(3,2)=-0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00462172;b(0,1)=-0.00295826;b(0,2)= 0.00060000;b(0,3)=-0.00009162;b(0,4)= 0.00001068;b(0,5)=-0.00000083;b(1,0)= 0.00017839;b(1,1)= 0.00005880;b(1,2)=-0.00006054;b(1,3)= 0.00001933;b(1,4)=-0.00000410;b(1,5)= 0.00000062;b(2,0)=-0.00002594;b(2,1)= 0.00001906;b(2,2)=-0.00000296;b(2,3)=-0.00000011;b(2,4)= 0.00000019;b(2,5)=-0.00000005;b(3,0)= 0.00000080;b(3,1)=-0.00000107;b(3,2)= 0.00000037;b(3,3)=-0.00000008;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000000;b(4,1)= 0.00000001;b(4,2)=-0.00000001;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.098350;a(0,1)=-0.021534;a(0,2)= 0.002051;a(0,3)=-0.000152;a(0,4)= 0.000005;a(0,5)= 0.000001;a(1,0)=-0.002897;a(1,1)= 0.001855;a(1,2)=-0.000376;a(1,3)= 0.000057;a(1,4)=-0.000007;a(1,5)= 0.000001;a(2,0)=-0.000028;a(2,1)=-0.000009;a(2,2)= 0.000010;a(2,3)=-0.000003;a(2,4)= 0.000001;a(2,5)= 0.000000;a(3,0)= 0.000003;a(3,1)=-0.000002;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00115912;b(0,1)=-0.00079617;b(0,2)= 0.00018434;b(0,3)=-0.00003515;b(0,4)= 0.00000602;b(0,5)=-0.00000089;b(1,0)= 0.00011906;b(1,1)=-0.00006421;b(1,2)= 0.00000923;b(1,3)=-0.00000036;b(1,4)=-0.00000026;b(1,5)= 0.00000009;b(2,0)=-0.00000409;b(2,1)= 0.00000408;b(2,2)=-0.00000127;b(2,3)= 0.00000030;b(2,4)=-0.00000006;b(2,5)= 0.00000001;b(3,0)=-0.00000001;b(3,1)=-0.00000002;b(3,2)= 0.00000002;b(3,3)=-0.00000001;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000000;b(4,1)= 0.00000000;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.062147;a(0,1)=-0.014699;a(0,2)= 0.001672;a(0,3)=-0.000182;a(0,4)= 0.000019;a(0,5)=-0.000002;a(1,0)=-0.000726;a(1,1)= 0.000499;a(1,2)=-0.000116;a(1,3)= 0.000022;a(1,4)=-0.000004;a(1,5)= 0.000001;a(2,0)=-0.000019;a(2,1)= 0.000010;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000000;a(3,1)= 0.000000;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000;
	else
		write(*,*)'error'
	endif
elseif(x<=7.5)then
	ax(1)=5.0;ax(2)=7.5
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.00956080;b(0,1)= 0.00836666;b(0,2)=-0.00110877;b(0,3)= 0.00008028;b(0,4)=-0.00000484;b(0,5)= 0.00000026;b(1,0)=-0.00378862;b(1,1)=-0.00329406;b(1,2)= 0.00046304;b(1,3)=-0.00002974;b(1,4)= 0.00000168;b(1,5)=-0.00000010;b(2,0)= 0.00056501;b(2,1)= 0.00048697;b(2,2)=-0.00007378;b(2,3)= 0.00000405;b(2,4)=-0.00000019;b(2,5)= 0.00000001;b(3,0)=-0.00007481;b(3,1)=-0.00006377;b(3,2)= 0.00001055;b(3,3)=-0.00000048;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000939;b(4,1)= 0.00000789;b(4,2)=-0.00000144;b(4,3)= 0.00000005;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000110;b(5,1)=-0.00000091;b(5,2)= 0.00000018;b(5,3)=-0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.057771;a(0,1)= 0.050800;a(0,2)=-0.006428;a(0,3)= 0.000510;a(0,4)=-0.000032;a(0,5)= 0.000002;a(1,0)=-0.011598;a(1,1)=-0.010154;a(1,2)= 0.001340;a(1,3)=-0.000098;a(1,4)= 0.000006;a(1,5)= 0.000000;a(2,0)= 0.001161;a(2,1)= 0.001009;a(2,2)=-0.000141;a(2,3)= 0.000009;a(2,4)=-0.000001;a(2,5)= 0.000000;a(3,0)=-0.000116;a(3,1)=-0.000100;a(3,2)= 0.000015;a(3,3)=-0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000012;a(4,1)= 0.000010;a(4,2)=-0.000002;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)=-0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.01975483;b(0,1)= 0.00234799;b(0,2)=-0.00047808;b(0,3)= 0.00003260;b(0,4)=-0.00000167;b(0,5)= 0.00000009;b(1,0)=-0.00756019;b(1,1)=-0.00068045;b(1,2)= 0.00021738;b(1,3)=-0.00001400;b(1,4)= 0.00000053;b(1,5)=-0.00000003;b(2,0)= 0.00107661;b(2,1)= 0.00005515;b(2,2)=-0.00003689;b(2,3)= 0.00000244;b(2,4)=-0.00000005;b(2,5)= 0.00000000;b(3,0)=-0.00013465;b(3,1)=-0.00000032;b(3,2)= 0.00000543;b(3,3)=-0.00000040;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00001576;b(4,1)=-0.00000092;b(4,2)=-0.00000073;b(4,3)= 0.00000007;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000171;b(5,1)= 0.00000022;b(5,2)= 0.00000009;b(5,3)=-0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.122557;a(0,1)= 0.017212;a(0,2)=-0.002519;a(0,3)= 0.000191;a(0,4)=-0.000011;a(0,5)= 0.000001;a(1,0)=-0.024021;a(1,1)=-0.002901;a(1,2)= 0.000575;a(1,3)=-0.000039;a(1,4)= 0.000002;a(1,5)= 0.000000;a(2,0)= 0.002320;a(2,1)= 0.000213;a(2,2)=-0.000066;a(2,3)= 0.000004;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000221;a(3,1)=-0.000012;a(3,2)= 0.000008;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000021;a(4,1)= 0.000000;a(4,2)=-0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000002;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.02161140;b(0,1)=-0.00026540;b(0,2)=-0.00020498;b(0,3)= 0.00001538;b(0,4)=-0.00000067;b(0,5)= 0.00000003;b(1,0)=-0.00763127;b(1,1)= 0.00050242;b(1,2)= 0.00008857;b(1,3)=-0.00000807;b(1,4)= 0.00000028;b(1,5)= 0.00000000;b(2,0)= 0.00097502;b(2,1)=-0.00013630;b(2,2)=-0.00001237;b(2,3)= 0.00000163;b(2,4)=-0.00000006;b(2,5)= 0.00000000;b(3,0)=-0.00010615;b(3,1)= 0.00002526;b(3,2)= 0.00000119;b(3,3)=-0.00000028;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00001040;b(4,1)=-0.00000386;b(4,2)=-0.00000004;b(4,3)= 0.00000004;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000091;b(5,1)= 0.00000049;b(5,2)=-0.00000001;b(5,3)=-0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.142372;a(0,1)= 0.003835;a(0,2)=-0.001026;a(0,3)= 0.000075;a(0,4)=-0.000004;a(0,5)= 0.000000;a(1,0)=-0.026405;a(1,1)= 0.000247;a(1,2)= 0.000248;a(1,3)=-0.000018;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.002352;a(2,1)=-0.000149;a(2,2)=-0.000027;a(2,3)= 0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000201;a(3,1)= 0.000028;a(3,2)= 0.000003;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000017;a(4,1)=-0.000004;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.01991842;b(0,1)=-0.00131648;b(0,2)=-0.00007105;b(0,3)= 0.00000774;b(0,4)=-0.00000034;b(0,5)= 0.00000001;b(1,0)=-0.00617508;b(1,1)= 0.00089330;b(1,2)= 0.00001586;b(1,3)=-0.00000423;b(1,4)= 0.00000020;b(1,5)= 0.00000000;b(2,0)= 0.00065457;b(2,1)=-0.00017234;b(2,2)= 0.00000180;b(2,3)= 0.00000074;b(2,4)=-0.00000005;b(2,5)= 0.00000000;b(3,0)=-0.00005441;b(3,1)= 0.00002466;b(3,2)=-0.00000100;b(3,3)=-0.00000009;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000345;b(4,1)=-0.00000287;b(4,2)= 0.00000023;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000013;b(5,1)= 0.00000027;b(5,2)=-0.00000003;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.144050;a(0,1)=-0.001657;a(0,2)=-0.000422;a(0,3)= 0.000032;a(0,4)=-0.000002;a(0,5)= 0.000000;a(1,0)=-0.024489;a(1,1)= 0.001538;a(1,2)= 0.000090;a(1,3)=-0.000009;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.001913;a(2,1)=-0.000271;a(2,2)=-0.000005;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000136;a(3,1)= 0.000035;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000008;a(4,1)=-0.000004;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.01642763;b(0,1)=-0.00211690;b(0,2)= 0.00000178;b(0,3)= 0.00000743;b(0,4)=-0.00000055;b(0,5)= 0.00000002;b(1,0)=-0.00411869;b(1,1)= 0.00111886;b(1,2)=-0.00003573;b(1,3)=-0.00000296;b(1,4)= 0.00000036;b(1,5)=-0.00000002;b(2,0)= 0.00030580;b(2,1)=-0.00016728;b(2,2)= 0.00001186;b(2,3)= 0.00000013;b(2,4)=-0.00000007;b(2,5)= 0.00000000;b(3,0)=-0.00001119;b(3,1)= 0.00001737;b(3,2)=-0.00000222;b(3,3)= 0.00000007;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)=-0.00000074;b(4,1)=-0.00000124;b(4,2)= 0.00000029;b(4,3)=-0.00000002;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000019;b(5,1)= 0.00000004;b(5,2)=-0.00000003;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.136908;a(0,1)=-0.005364;a(0,2)=-0.000238;a(0,3)= 0.000030;a(0,4)=-0.000002;a(0,5)= 0.000000;a(1,0)=-0.020343;a(1,1)= 0.002542;a(1,2)= 0.000005;a(1,3)=-0.000009;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.001284;a(2,1)=-0.000344;a(2,2)= 0.000010;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000064;a(3,1)= 0.000035;a(3,2)=-0.000002;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)=-0.000003;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.01080872;b(0,1)=-0.00329878;b(0,2)= 0.00021662;b(0,3)= 0.00000233;b(0,4)=-0.00000212;b(0,5)= 0.00000025;b(1,0)=-0.00166657;b(1,1)= 0.00120201;b(1,2)=-0.00016210;b(1,3)= 0.00000834;b(1,4)= 0.00000068;b(1,5)=-0.00000017;b(2,0)= 0.00001750;b(2,1)=-0.00010267;b(2,2)= 0.00002712;b(2,3)=-0.00000309;b(2,4)= 0.00000009;b(2,5)= 0.00000002;b(3,0)= 0.00000977;b(3,1)= 0.00000256;b(3,2)=-0.00000257;b(3,3)= 0.00000054;b(3,4)=-0.00000005;b(3,5)= 0.00000000;b(4,0)=-0.00000132;b(4,1)= 0.00000056;b(4,2)= 0.00000010;b(4,3)=-0.00000006;b(4,4)= 0.00000001;b(4,5)= 0.00000000;b(5,0)= 0.00000010;b(5,1)=-0.00000009;b(5,2)= 0.00000001;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.119048;a(0,1)=-0.012210;a(0,2)= 0.000137;a(0,3)= 0.000052;a(0,4)=-0.000007;a(0,5)= 0.000001;a(1,0)=-0.013500;a(1,1)= 0.004059;a(1,2)=-0.000254;a(1,3)=-0.000005;a(1,4)= 0.000003;a(1,5)= 0.000000;a(2,0)= 0.000524;a(2,1)=-0.000375;a(2,2)= 0.000050;a(2,3)=-0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000004;a(3,1)= 0.000022;a(3,2)=-0.000006;a(3,3)= 0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)=-0.000002;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000;
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00474231;b(0,1)=-0.00253628;b(0,2)= 0.00040259;b(0,3)=-0.00004340;b(0,4)= 0.00000255;b(0,5)= 0.00000012;b(1,0)=-0.00008925;b(1,1)= 0.00036747;b(1,2)=-0.00013086;b(1,3)= 0.00002644;b(1,4)=-0.00000352;b(1,5)= 0.00000027;b(2,0)=-0.00005151;b(2,1)= 0.00001737;b(2,2)= 0.00000443;b(2,3)=-0.00000265;b(2,4)= 0.00000066;b(2,5)=-0.00000010;b(3,0)= 0.00000487;b(3,1)=-0.00000457;b(3,2)= 0.00000087;b(3,3)= 0.00000000;b(3,4)=-0.00000005;b(3,5)= 0.00000001;b(4,0)=-0.00000020;b(4,1)= 0.00000034;b(4,2)=-0.00000013;b(4,3)= 0.00000003;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)=-0.00000001;b(5,2)= 0.00000001;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.089450;a(0,1)=-0.016397;a(0,2)= 0.001143;a(0,3)=-0.000037;a(0,4)=-0.000005;a(0,5)= 0.000001;a(1,0)=-0.005960;a(1,1)= 0.003181;a(1,2)=-0.000500;a(1,3)= 0.000053;a(1,4)=-0.000003;a(1,5)= 0.000000;a(2,0)= 0.000029;a(2,1)=-0.000116;a(2,2)= 0.000041;a(2,3)=-0.000008;a(2,4)= 0.000001;a(2,5)= 0.000000;a(3,0)= 0.000011;a(3,1)=-0.000004;a(3,2)=-0.000001;a(3,3)= 0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)=-0.000001;a(4,1)= 0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00143228;b(0,1)=-0.00090964;b(0,2)= 0.00018940;b(0,3)=-0.00003152;b(0,4)= 0.00000450;b(0,5)=-0.00000054;b(1,0)= 0.00014291;b(1,1)=-0.00003982;b(1,2)=-0.00000638;b(1,3)= 0.00000430;b(1,4)=-0.00000126;b(1,5)= 0.00000025;b(2,0)=-0.00001501;b(2,1)= 0.00001300;b(2,2)=-0.00000323;b(2,3)= 0.00000053;b(2,4)=-0.00000004;b(2,5)= 0.00000000;b(3,0)= 0.00000030;b(3,1)=-0.00000051;b(3,2)= 0.00000023;b(3,3)=-0.00000007;b(3,4)= 0.00000002;b(3,5)= 0.00000000;b(4,0)= 0.00000002;b(4,1)=-0.00000001;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.059651;a(0,1)=-0.013061;a(0,2)= 0.001315;a(0,3)=-0.000119;a(0,4)= 0.000009;a(0,5)= 0.000000;a(1,0)=-0.001800;a(1,1)= 0.001145;a(1,2)=-0.000239;a(1,3)= 0.000040;a(1,4)=-0.000006;a(1,5)= 0.000001;a(2,0)=-0.000045;a(2,1)= 0.000012;a(2,2)= 0.000002;a(2,3)=-0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000003;a(3,1)=-0.000003;a(3,2)= 0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
elseif(x<=12.0)then
	ax(1)=7.5;ax(2)=12.0
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.00403518;b(0,1)= 0.00354542;b(0,2)=-0.00045216;b(0,3)= 0.00003530;b(0,4)=-0.00000219;b(0,5)= 0.00000011;b(1,0)=-0.00185504;b(1,1)=-0.00162550;b(1,2)= 0.00021276;b(1,3)=-0.00001577;b(1,4)= 0.00000097;b(1,5)=-0.00000005;b(2,0)= 0.00032226;b(2,1)= 0.00028135;b(2,2)=-0.00003811;b(2,3)= 0.00000263;b(2,4)=-0.00000016;b(2,5)= 0.00000001;b(3,0)=-0.00004983;b(3,1)=-0.00004330;b(3,2)= 0.00000611;b(3,3)=-0.00000039;b(3,4)= 0.00000002;b(3,5)= 0.00000000;b(4,0)= 0.00000735;b(4,1)= 0.00000635;b(4,2)=-0.00000094;b(4,3)= 0.00000005;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000100;b(5,1)=-0.00000086;b(5,2)= 0.00000013;b(5,3)=-0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.037363;a(0,1)= 0.032897;a(0,2)=-0.004110;a(0,3)= 0.000334;a(0,4)=-0.000021;a(0,5)= 0.000001;a(1,0)=-0.008717;a(1,1)=-0.007661;a(1,2)= 0.000974;a(1,3)=-0.000076;a(1,4)= 0.000005;a(1,5)= 0.000000;a(2,0)= 0.001015;a(2,1)= 0.000890;a(2,2)=-0.000116;a(2,3)= 0.000009;a(2,4)=-0.000001;a(2,5)= 0.000000;a(3,0)=-0.000118;a(3,1)=-0.000103;a(3,2)= 0.000014;a(3,3)=-0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000014;a(4,1)= 0.000012;a(4,2)=-0.000002;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000002;a(5,1)=-0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000;
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.00852069;b(0,1)= 0.00116356;b(0,2)=-0.00018135;b(0,3)= 0.00001327;b(0,4)=-0.00000079;b(0,5)= 0.00000004;b(1,0)=-0.00385817;b(1,1)=-0.00047841;b(1,2)= 0.00009065;b(1,3)=-0.00000613;b(1,4)= 0.00000034;b(1,5)=-0.00000002;b(2,0)= 0.00065676;b(2,1)= 0.00007031;b(2,2)=-0.00001731;b(2,3)= 0.00000109;b(2,4)=-0.00000005;b(2,5)= 0.00000000;b(3,0)=-0.00009901;b(3,1)=-0.00000849;b(3,2)= 0.00000295;b(3,3)=-0.00000018;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00001415;b(4,1)= 0.00000084;b(4,2)=-0.00000048;b(4,3)= 0.00000003;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000187;b(5,1)=-0.00000006;b(5,2)= 0.00000007;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.079836;a(0,1)= 0.011681;a(0,2)=-0.001558;a(0,3)= 0.000123;a(0,4)=-0.000008;a(0,5)= 0.000000;a(1,0)=-0.018433;a(1,1)=-0.002539;a(1,2)= 0.000389;a(1,3)=-0.000029;a(1,4)= 0.000002;a(1,5)= 0.000000;a(2,0)= 0.002115;a(2,1)= 0.000264;a(2,2)=-0.000049;a(2,3)= 0.000003;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000241;a(3,1)=-0.000026;a(3,2)= 0.000006;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000028;a(4,1)= 0.000002;a(4,2)=-0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000003;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.00978206;b(0,1)= 0.00018367;b(0,2)=-0.00007737;b(0,3)= 0.00000530;b(0,4)=-0.00000029;b(0,5)= 0.00000001;b(1,0)=-0.00427143;b(1,1)= 0.00002399;b(1,2)= 0.00004084;b(1,3)=-0.00000271;b(1,4)= 0.00000012;b(1,5)=-0.00000001;b(2,0)= 0.00069234;b(2,1)=-0.00002695;b(2,2)=-0.00000792;b(2,3)= 0.00000056;b(2,4)=-0.00000002;b(2,5)= 0.00000000;b(3,0)=-0.00009811;b(3,1)= 0.00000802;b(3,2)= 0.00000130;b(3,3)=-0.00000011;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00001297;b(4,1)=-0.00000178;b(4,2)=-0.00000019;b(4,3)= 0.00000002;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000157;b(5,1)= 0.00000032;b(5,2)= 0.00000002;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.094269;a(0,1)= 0.003534;a(0,2)=-0.000612;a(0,3)= 0.000046;a(0,4)=-0.000003;a(0,5)= 0.000000;a(1,0)=-0.021231;a(1,1)=-0.000444;a(1,2)= 0.000165;a(1,3)=-0.000011;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.002348;a(2,1)=-0.000009;a(2,2)=-0.000022;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000255;a(3,1)= 0.000009;a(3,2)= 0.000003;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000027;a(4,1)=-0.000002;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000003;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.00968832;b(0,1)=-0.00024155;b(0,2)=-0.00003400;b(0,3)= 0.00000236;b(0,4)=-0.00000011;b(0,5)= 0.00000001;b(1,0)=-0.00398078;b(1,1)= 0.00024706;b(1,2)= 0.00001718;b(1,3)=-0.00000141;b(1,4)= 0.00000005;b(1,5)= 0.00000000;b(2,0)= 0.00059313;b(2,1)=-0.00006797;b(2,2)=-0.00000274;b(2,3)= 0.00000033;b(2,4)=-0.00000001;b(2,5)= 0.00000000;b(3,0)=-0.00007520;b(3,1)= 0.00001403;b(3,2)= 0.00000028;b(3,3)=-0.00000007;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000855;b(4,1)=-0.00000248;b(4,2)= 0.00000000;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000086;b(5,1)= 0.00000037;b(5,2)=-0.00000001;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.097776;a(0,1)= 0.000269;a(0,2)=-0.000253;a(0,3)= 0.000018;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.021131;a(1,1)= 0.000467;a(1,2)= 0.000073;a(1,3)=-0.000005;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.002197;a(2,1)=-0.000131;a(2,2)=-0.000010;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000219;a(3,1)= 0.000025;a(3,2)= 0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000021;a(4,1)=-0.000004;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000002;a(5,1)= 0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.00885407;b(0,1)=-0.00058408;b(0,2)=-0.00002044;b(0,3)= 0.00000254;b(0,4)=-0.00000013;b(0,5)= 0.00000001;b(1,0)=-0.00328189;b(1,1)= 0.00044311;b(1,2)= 0.00000536;b(1,3)=-0.00000167;b(1,4)= 0.00000009;b(1,5)= 0.00000000;b(2,0)= 0.00042076;b(2,1)=-0.00010141;b(2,2)= 0.00000096;b(2,3)= 0.00000037;b(2,4)=-0.00000002;b(2,5)= 0.00000000;b(3,0)=-0.00004275;b(3,1)= 0.00001766;b(3,2)=-0.00000060;b(3,3)=-0.00000006;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000335;b(4,1)=-0.00000257;b(4,2)= 0.00000017;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000015;b(5,1)= 0.00000031;b(5,2)=-0.00000003;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.096347;a(0,1)=-0.001653;a(0,2)=-0.000174;a(0,3)= 0.000015;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.019448;a(1,1)= 0.001200;a(1,2)= 0.000047;a(1,3)=-0.000005;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.001822;a(2,1)=-0.000239;a(2,2)=-0.000003;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000157;a(3,1)= 0.000037;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000012;a(4,1)=-0.000005;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.00696482;b(0,1)=-0.00127337;b(0,2)= 0.00002023;b(0,3)= 0.00000573;b(0,4)=-0.00000069;b(0,5)= 0.00000005;b(1,0)=-0.00202634;b(1,1)= 0.00077077;b(1,2)=-0.00004108;b(1,3)=-0.00000258;b(1,4)= 0.00000056;b(1,5)=-0.00000004;b(2,0)= 0.00017199;b(2,1)=-0.00013455;b(2,2)= 0.00001425;b(2,3)=-0.00000002;b(2,4)=-0.00000014;b(2,5)= 0.00000001;b(3,0)=-0.00000619;b(3,1)= 0.00001636;b(3,2)=-0.00000307;b(3,3)= 0.00000016;b(3,4)= 0.00000002;b(3,5)= 0.00000000;b(4,0)=-0.00000090;b(4,1)=-0.00000133;b(4,2)= 0.00000048;b(4,3)=-0.00000005;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000023;b(5,1)= 0.00000004;b(5,2)=-0.00000006;b(5,3)= 0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.089387;a(0,1)=-0.005312;a(0,2)=-0.000129;a(0,3)= 0.000031;a(0,4)=-0.000003;a(0,5)= 0.000000;a(1,0)=-0.015477;a(1,1)= 0.002714;a(1,2)=-0.000029;a(1,3)=-0.000013;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.001136;a(2,1)=-0.000424;a(2,2)= 0.000021;a(2,3)= 0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000065;a(3,1)= 0.000050;a(3,2)=-0.000005;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)=-0.000005;a(4,2)= 0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00402754;b(0,1)=-0.00153693;b(0,2)= 0.00014813;b(0,3)=-0.00000520;b(0,4)=-0.00000085;b(0,5)= 0.00000018;b(1,0)=-0.00060094;b(1,1)= 0.00057790;b(1,2)=-0.00010934;b(1,3)= 0.00001054;b(1,4)=-0.00000012;b(1,5)=-0.00000013;b(2,0)=-0.00001327;b(2,1)=-0.00004320;b(2,2)= 0.00001807;b(2,3)=-0.00000327;b(2,4)= 0.00000028;b(2,5)= 0.00000001;b(3,0)= 0.00000789;b(3,1)=-0.00000144;b(3,2)=-0.00000143;b(3,3)= 0.00000054;b(3,4)=-0.00000009;b(3,5)= 0.00000001;b(4,0)=-0.00000102;b(4,1)= 0.00000081;b(4,2)=-0.00000004;b(4,3)=-0.00000005;b(4,4)= 0.00000002;b(4,5)= 0.00000000;b(5,0)= 0.00000008;b(5,1)=-0.00000011;b(5,2)= 0.00000003;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.074117;a(0,1)=-0.009584;a(0,2)= 0.000302;a(0,3)= 0.000023;a(0,4)=-0.000005;a(0,5)= 0.000001;a(1,0)=-0.009077;a(1,1)= 0.003409;a(1,2)=-0.000313;a(1,3)= 0.000008;a(1,4)= 0.000002;a(1,5)= 0.000000;a(2,0)= 0.000342;a(2,1)=-0.000326;a(2,2)= 0.000061;a(2,3)=-0.000006;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000005;a(3,1)= 0.000016;a(3,2)=-0.000007;a(3,3)= 0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)=-0.000002;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00161436;b(0,1)=-0.00085908;b(0,2)= 0.00014085;b(0,3)=-0.00001701;b(0,4)= 0.00000150;b(0,5)=-0.00000007;b(1,0)= 0.00002425;b(1,1)= 0.00009484;b(1,2)=-0.00004053;b(1,3)= 0.00000929;b(1,4)=-0.00000150;b(1,5)= 0.00000018;b(2,0)=-0.00002726;b(2,1)= 0.00001436;b(2,2)=-0.00000062;b(2,3)=-0.00000059;b(2,4)= 0.00000022;b(2,5)=-0.00000005;b(3,0)= 0.00000229;b(3,1)=-0.00000245;b(3,2)= 0.00000066;b(3,3)=-0.00000009;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)=-0.00000006;b(4,1)= 0.00000014;b(4,2)=-0.00000007;b(4,3)= 0.00000002;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000001;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.054170;a(0,1)=-0.009908;a(0,2)= 0.000739;a(0,3)=-0.000037;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.003663;a(1,1)= 0.001949;a(1,2)=-0.000318;a(1,3)= 0.000038;a(1,4)=-0.000003;a(1,5)= 0.000000;a(2,0)=-0.000012;a(2,1)=-0.000055;a(2,2)= 0.000023;a(2,3)=-0.000005;a(2,4)= 0.000001;a(2,5)= 0.000000;a(3,0)= 0.000010;a(3,1)=-0.000005;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)=-0.000001;a(4,1)= 0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
elseif(x<=22.0)then
	ax(1)=12.0;ax(2)=22.0
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.00140443;b(0,1)= 0.00123636;b(0,2)=-0.00015471;b(0,3)= 0.00001254;b(0,4)=-0.00000078;b(0,5)= 0.00000004;b(1,0)=-0.00082502;b(1,1)=-0.00072562;b(1,2)= 0.00009163;b(1,3)=-0.00000729;b(1,4)= 0.00000045;b(1,5)=-0.00000002;b(2,0)= 0.00018437;b(2,1)= 0.00016195;b(2,2)=-0.00002070;b(2,3)= 0.00000161;b(2,4)=-0.00000010;b(2,5)= 0.00000001;b(3,0)=-0.00003678;b(3,1)=-0.00003226;b(3,2)= 0.00000419;b(3,3)=-0.00000032;b(3,4)= 0.00000002;b(3,5)= 0.00000000;b(4,0)= 0.00000709;b(4,1)= 0.00000620;b(4,2)=-0.00000082;b(4,3)= 0.00000006;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000123;b(5,1)=-0.00000108;b(5,2)= 0.00000015;b(5,3)=-0.00000001;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.021835;a(0,1)= 0.019237;a(0,2)=-0.002388;a(0,3)= 0.000197;a(0,4)=-0.000012;a(0,5)= 0.000001;a(1,0)=-0.006561;a(1,1)=-0.005777;a(1,2)= 0.000722;a(1,3)=-0.000059;a(1,4)= 0.000004;a(1,5)= 0.000000;a(2,0)= 0.000985;a(2,1)= 0.000867;a(2,2)=-0.000109;a(2,3)= 0.000009;a(2,4)=-0.000001;a(2,5)= 0.000000;a(3,0)=-0.000148;a(3,1)=-0.000130;a(3,2)= 0.000017;a(3,3)=-0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000023;a(4,1)= 0.000020;a(4,2)=-0.000003;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000003;a(5,1)=-0.000003;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.00299799;b(0,1)= 0.00043618;b(0,2)=-0.00005900;b(0,3)= 0.00000463;b(0,4)=-0.00000029;b(0,5)= 0.00000001;b(1,0)=-0.00175186;b(1,1)=-0.00024722;b(1,2)= 0.00003590;b(1,3)=-0.00000270;b(1,4)= 0.00000017;b(1,5)=-0.00000001;b(2,0)= 0.00038870;b(2,1)= 0.00005254;b(2,2)=-0.00000839;b(2,3)= 0.00000060;b(2,4)=-0.00000004;b(2,5)= 0.00000000;b(3,0)=-0.00007685;b(3,1)=-0.00000981;b(3,2)= 0.00000176;b(3,3)=-0.00000012;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00001464;b(4,1)= 0.00000173;b(4,2)=-0.00000036;b(4,3)= 0.00000002;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000252;b(5,1)=-0.00000027;b(5,2)= 0.00000007;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.046822;a(0,1)= 0.006987;a(0,2)=-0.000889;a(0,3)= 0.000072;a(0,4)=-0.000005;a(0,5)= 0.000000;a(1,0)=-0.014018;a(1,1)=-0.002050;a(1,2)= 0.000274;a(1,3)=-0.000022;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.002094;a(2,1)= 0.000297;a(2,2)=-0.000043;a(2,3)= 0.000003;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000312;a(3,1)=-0.000042;a(3,2)= 0.000007;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000047;a(4,1)= 0.000006;a(4,2)=-0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000007;a(5,1)=-0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.00353091;b(0,1)= 0.00012601;b(0,2)=-0.00002355;b(0,3)= 0.00000173;b(0,4)=-0.00000010;b(0,5)= 0.00000001;b(1,0)=-0.00203684;b(1,1)=-0.00005496;b(1,2)= 0.00001506;b(1,3)=-0.00000103;b(1,4)= 0.00000006;b(1,5)= 0.00000000;b(2,0)= 0.00044408;b(2,1)= 0.00000671;b(2,2)=-0.00000370;b(2,3)= 0.00000024;b(2,4)=-0.00000001;b(2,5)= 0.00000000;b(3,0)=-0.00008588;b(3,1)= 0.00000000;b(3,2)= 0.00000081;b(3,3)=-0.00000005;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00001590;b(4,1)=-0.00000031;b(4,2)=-0.00000017;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000266;b(5,1)= 0.00000011;b(5,2)= 0.00000003;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.055755;a(0,1)= 0.002401;a(0,2)=-0.000337;a(0,3)= 0.000027;a(0,4)=-0.000002;a(0,5)= 0.000000;a(1,0)=-0.016544;a(1,1)=-0.000613;a(1,2)= 0.000108;a(1,3)=-0.000008;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.002439;a(2,1)= 0.000069;a(2,2)=-0.000018;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000357;a(3,1)=-0.000006;a(3,2)= 0.000003;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000053;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000008;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.00364441;b(0,1)=-0.00000144;b(0,2)=-0.00001015;b(0,3)= 0.00000067;b(0,4)=-0.00000004;b(0,5)= 0.00000000;b(1,0)=-0.00205636;b(1,1)= 0.00002871;b(1,2)= 0.00000691;b(1,3)=-0.00000042;b(1,4)= 0.00000002;b(1,5)= 0.00000000;b(2,0)= 0.00043496;b(2,1)=-0.00001422;b(2,2)=-0.00000176;b(2,3)= 0.00000011;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00008092;b(3,1)= 0.00000461;b(3,2)= 0.00000038;b(3,3)=-0.00000002;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00001425;b(4,1)=-0.00000127;b(4,2)=-0.00000008;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000226;b(5,1)= 0.00000028;b(5,2)= 0.00000001;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.058630;a(0,1)= 0.000642;a(0,2)=-0.000132;a(0,3)= 0.000010;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.017135;a(1,1)=-0.000028;a(1,2)= 0.000046;a(1,3)=-0.000003;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.002469;a(2,1)=-0.000030;a(2,2)=-0.000008;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000351;a(3,1)= 0.000011;a(3,2)= 0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000050;a(4,1)=-0.000003;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000007;a(5,1)= 0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.00355635;b(0,1)=-0.00008542;b(0,2)=-0.00000767;b(0,3)= 0.00000059;b(0,4)=-0.00000004;b(0,5)= 0.00000000;b(1,0)=-0.00193068;b(1,1)= 0.00009630;b(1,2)= 0.00000539;b(1,3)=-0.00000043;b(1,4)= 0.00000002;b(1,5)= 0.00000000;b(2,0)= 0.00038701;b(2,1)=-0.00003347;b(2,2)=-0.00000127;b(2,3)= 0.00000012;b(2,4)=-0.00000001;b(2,5)= 0.00000000;b(3,0)=-0.00006709;b(3,1)= 0.00000912;b(3,2)= 0.00000022;b(3,3)=-0.00000003;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00001073;b(4,1)=-0.00000221;b(4,2)=-0.00000002;b(4,3)= 0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000152;b(5,1)= 0.00000044;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.059071;a(0,1)=-0.000174;a(0,2)=-0.000088;a(0,3)= 0.000008;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.016814;a(1,1)= 0.000343;a(1,2)= 0.000035;a(1,3)=-0.000003;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.002330;a(2,1)=-0.000109;a(2,2)=-0.000006;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000314;a(3,1)= 0.000026;a(3,2)= 0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000042;a(4,1)=-0.000006;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000005;a(5,1)= 0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.00320405;b(0,1)=-0.00026919;b(0,2)=-0.00000790;b(0,3)= 0.00000141;b(0,4)=-0.00000011;b(0,5)= 0.00000001;b(1,0)=-0.00157972;b(1,1)= 0.00025405;b(1,2)= 0.00000299;b(1,3)=-0.00000130;b(1,4)= 0.00000009;b(1,5)=-0.00000001;b(2,0)= 0.00027559;b(2,1)=-0.00007653;b(2,2)= 0.00000070;b(2,3)= 0.00000041;b(2,4)=-0.00000003;b(2,5)= 0.00000000;b(3,0)=-0.00003924;b(3,1)= 0.00001803;b(3,2)=-0.00000063;b(3,3)=-0.00000010;b(3,4)= 0.00000001;b(3,5)= 0.00000000;b(4,0)= 0.00000461;b(4,1)=-0.00000367;b(4,2)= 0.00000025;b(4,3)= 0.00000002;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000040;b(5,1)= 0.00000062;b(5,2)=-0.00000007;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.057473;a(0,1)=-0.001453;a(0,2)=-0.000097;a(0,3)= 0.000012;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.015331;a(1,1)= 0.001155;a(1,2)= 0.000041;a(1,3)=-0.000006;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.001926;a(2,1)=-0.000295;a(2,2)=-0.000005;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000226;a(3,1)= 0.000061;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000024;a(4,1)=-0.000011;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000002;a(5,1)= 0.000002;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00240765;b(0,1)=-0.00050851;b(0,2)= 0.00001440;b(0,3)= 0.00000197;b(0,4)=-0.00000029;b(0,5)= 0.00000002;b(1,0)=-0.00090872;b(1,1)= 0.00039050;b(1,2)=-0.00002610;b(1,3)=-0.00000106;b(1,4)= 0.00000034;b(1,5)=-0.00000003;b(2,0)= 0.00010058;b(2,1)=-0.00008843;b(2,2)= 0.00001110;b(2,3)=-0.00000017;b(2,4)=-0.00000011;b(2,5)= 0.00000001;b(3,0)=-0.00000464;b(3,1)= 0.00001403;b(3,2)=-0.00000308;b(3,3)= 0.00000020;b(3,4)= 0.00000002;b(3,5)=-0.00000001;b(4,0)=-0.00000099;b(4,1)=-0.00000147;b(4,2)= 0.00000064;b(4,3)=-0.00000008;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000032;b(5,1)= 0.00000005;b(5,2)=-0.00000010;b(5,3)= 0.00000002;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.052382;a(0,1)=-0.003610;a(0,2)=-0.000041;a(0,3)= 0.000015;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.011787;a(1,1)= 0.002321;a(1,2)=-0.000044;a(1,3)=-0.000010;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.001130;a(2,1)=-0.000471;a(2,2)= 0.000029;a(2,3)= 0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000085;a(3,1)= 0.000072;a(3,2)=-0.000009;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000003;a(4,1)=-0.000009;a(4,2)= 0.000002;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000001;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00137080;b(0,1)=-0.00049533;b(0,2)= 0.00004718;b(0,3)=-0.00000212;b(0,4)=-0.00000014;b(0,5)= 0.00000004;b(1,0)=-0.00025992;b(1,1)= 0.00023749;b(1,2)=-0.00004332;b(1,3)= 0.00000426;b(1,4)=-0.00000012;b(1,5)=-0.00000003;b(2,0)=-0.00000969;b(2,1)=-0.00002182;b(2,2)= 0.00000905;b(2,3)=-0.00000164;b(2,4)= 0.00000015;b(2,5)= 0.00000000;b(3,0)= 0.00000631;b(3,1)=-0.00000146;b(3,2)=-0.00000085;b(3,3)= 0.00000034;b(3,4)=-0.00000006;b(3,5)= 0.00000001;b(4,0)=-0.00000104;b(4,1)= 0.00000083;b(4,2)=-0.00000007;b(4,3)=-0.00000003;b(4,4)= 0.00000001;b(4,5)= 0.00000000;b(5,0)= 0.00000010;b(5,1)=-0.00000014;b(5,2)= 0.00000004;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.043305;a(0,1)=-0.005280;a(0,2)= 0.000170;a(0,3)= 0.000008;a(0,4)=-0.000002;a(0,5)= 0.000000;a(1,0)=-0.006878;a(1,1)= 0.002422;a(1,2)=-0.000213;a(1,3)= 0.000006;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000333;a(2,1)=-0.000299;a(2,2)= 0.000053;a(2,3)=-0.000005;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)= 0.000007;a(3,1)= 0.000019;a(3,2)=-0.000008;a(3,3)= 0.000001;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)=-0.000004;a(4,1)= 0.000001;a(4,2)= 0.000001;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
elseif(x<=33.0)then
	ax(1)=22.0;ax(2)=33.0
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.00049874;b(0,1)= 0.00043939;b(0,2)=-0.00005456;b(0,3)= 0.00000449;b(0,4)=-0.00000028;b(0,5)= 0.00000001;b(1,0)=-0.00019940;b(1,1)=-0.00017561;b(1,2)= 0.00002188;b(1,3)=-0.00000179;b(1,4)= 0.00000011;b(1,5)=-0.00000001;b(2,0)= 0.00003009;b(2,1)= 0.00002649;b(2,2)=-0.00000332;b(2,3)= 0.00000027;b(2,4)=-0.00000002;b(2,5)= 0.00000000;b(3,0)=-0.00000404;b(3,1)=-0.00000356;b(3,2)= 0.00000045;b(3,3)=-0.00000004;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000052;b(4,1)= 0.00000045;b(4,2)=-0.00000006;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000006;b(5,1)=-0.00000005;b(5,2)= 0.00000001;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.013172;a(0,1)= 0.011607;a(0,2)=-0.001438;a(0,3)= 0.000119;a(0,4)=-0.000007;a(0,5)= 0.000000;a(1,0)=-0.002660;a(1,1)=-0.002344;a(1,2)= 0.000291;a(1,3)=-0.000024;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000269;a(2,1)= 0.000237;a(2,2)=-0.000029;a(2,3)= 0.000002;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000027;a(3,1)=-0.000024;a(3,2)= 0.000003;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000003;a(4,1)= 0.000002;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.00106935;b(0,1)= 0.00015946;b(0,2)=-0.00002032;b(0,3)= 0.00000165;b(0,4)=-0.00000010;b(0,5)= 0.00000001;b(1,0)=-0.00042670;b(1,1)=-0.00006295;b(1,2)= 0.00000824;b(1,3)=-0.00000066;b(1,4)= 0.00000004;b(1,5)= 0.00000000;b(2,0)= 0.00006423;b(2,1)= 0.00000934;b(2,2)=-0.00000127;b(2,3)= 0.00000010;b(2,4)=-0.00000001;b(2,5)= 0.00000000;b(3,0)=-0.00000860;b(3,1)=-0.00000123;b(3,2)= 0.00000017;b(3,3)=-0.00000001;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000109;b(4,1)= 0.00000015;b(4,2)=-0.00000002;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000013;b(5,1)=-0.00000002;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.028282;a(0,1)= 0.004252;a(0,2)=-0.000531;a(0,3)= 0.000044;a(0,4)=-0.000003;a(0,5)= 0.000000;a(1,0)=-0.005705;a(1,1)=-0.000851;a(1,2)= 0.000108;a(1,3)=-0.000009;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000575;a(2,1)= 0.000085;a(2,2)=-0.000011;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000058;a(3,1)=-0.000008;a(3,2)= 0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000006;a(4,1)= 0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.00127293;b(0,1)= 0.00005452;b(0,2)=-0.00000772;b(0,3)= 0.00000061;b(0,4)=-0.00000004;b(0,5)= 0.00000000;b(1,0)=-0.00050552;b(1,1)=-0.00002002;b(1,2)= 0.00000321;b(1,3)=-0.00000024;b(1,4)= 0.00000002;b(1,5)= 0.00000000;b(2,0)= 0.00007561;b(2,1)= 0.00000267;b(2,2)=-0.00000051;b(2,3)= 0.00000004;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00001005;b(3,1)=-0.00000030;b(3,2)= 0.00000007;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000127;b(4,1)= 0.00000003;b(4,2)=-0.00000001;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000015;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.033787;a(0,1)= 0.001529;a(0,2)=-0.000197;a(0,3)= 0.000016;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.006793;a(1,1)=-0.000293;a(1,2)= 0.000041;a(1,3)=-0.000003;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000681;a(2,1)= 0.000027;a(2,2)=-0.000004;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000068;a(3,1)=-0.000002;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000007;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.00133767;b(0,1)= 0.00001405;b(0,2)=-0.00000306;b(0,3)= 0.00000023;b(0,4)=-0.00000001;b(0,5)= 0.00000000;b(1,0)=-0.00052687;b(1,1)=-0.00000286;b(1,2)= 0.00000134;b(1,3)=-0.00000009;b(1,4)= 0.00000001;b(1,5)= 0.00000000;b(2,0)= 0.00007793;b(2,1)=-0.00000011;b(2,2)=-0.00000023;b(2,3)= 0.00000001;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00001021;b(3,1)= 0.00000010;b(3,2)= 0.00000003;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000126;b(4,1)=-0.00000003;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000015;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.035725;a(0,1)= 0.000511;a(0,2)=-0.000075;a(0,3)= 0.000006;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.007143;a(1,1)=-0.000078;a(1,2)= 0.000016;a(1,3)=-0.000001;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000710;a(2,1)= 0.000004;a(2,2)=-0.000002;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000070;a(3,1)= 0.000000;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000007;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.00134575;b(0,1)=-0.00000542;b(0,2)=-0.00000211;b(0,3)= 0.00000017;b(0,4)=-0.00000001;b(0,5)= 0.00000000;b(1,0)=-0.00052235;b(1,1)= 0.00000723;b(1,2)= 0.00000103;b(1,3)=-0.00000007;b(1,4)= 0.00000001;b(1,5)= 0.00000000;b(2,0)= 0.00007574;b(2,1)=-0.00000206;b(2,2)=-0.00000019;b(2,3)= 0.00000001;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000968;b(3,1)= 0.00000043;b(3,2)= 0.00000003;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000116;b(4,1)=-0.00000008;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000013;b(5,1)= 0.00000001;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.036335;a(0,1)= 0.000118;a(0,2)=-0.000046;a(0,3)= 0.000004;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.007193;a(1,1)= 0.000024;a(1,2)= 0.000011;a(1,3)=-0.000001;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000705;a(2,1)=-0.000009;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000068;a(3,1)= 0.000002;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000007;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)=-0.000001;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.00130230;b(0,1)=-0.00003917;b(0,2)=-0.00000279;b(0,3)= 0.00000026;b(0,4)=-0.00000003;b(0,5)= 0.00000000;b(1,0)=-0.00048665;b(1,1)= 0.00002921;b(1,2)= 0.00000155;b(1,3)=-0.00000014;b(1,4)= 0.00000001;b(1,5)= 0.00000000;b(2,0)= 0.00006698;b(2,1)=-0.00000683;b(2,2)=-0.00000028;b(2,3)= 0.00000003;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000800;b(3,1)= 0.00000126;b(3,2)= 0.00000004;b(3,3)=-0.00000001;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000088;b(4,1)=-0.00000021;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000009;b(5,1)= 0.00000003;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.036161;a(0,1)=-0.000300;a(0,2)=-0.000043;a(0,3)= 0.000006;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.006978;a(1,1)= 0.000197;a(1,2)= 0.000015;a(1,3)=-0.000001;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000658;a(2,1)=-0.000038;a(2,2)=-0.000002;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000061;a(3,1)= 0.000006;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000005;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00115748;b(0,1)=-0.00010637;b(0,2)=-0.00000283;b(0,3)= 0.00000044;b(0,4)=-0.00000002;b(0,5)= 0.00000000;b(1,0)=-0.00038670;b(1,1)= 0.00007034;b(1,2)= 0.00000065;b(1,3)=-0.00000035;b(1,4)= 0.00000002;b(1,5)= 0.00000000;b(2,0)= 0.00004534;b(2,1)=-0.00001446;b(2,2)= 0.00000020;b(2,3)= 0.00000008;b(2,4)=-0.00000001;b(2,5)= 0.00000000;b(3,0)=-0.00000431;b(3,1)= 0.00000231;b(3,2)=-0.00000010;b(3,3)=-0.00000001;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000033;b(4,1)=-0.00000031;b(4,2)= 0.00000003;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000002;b(5,1)= 0.00000004;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.034832;a(0,1)=-0.001041;a(0,2)=-0.000047;a(0,3)= 0.000004;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.006241;a(1,1)= 0.000545;a(1,2)= 0.000016;a(1,3)=-0.000002;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000526;a(2,1)=-0.000094;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000041;a(3,1)= 0.000013;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000003;a(4,1)=-0.000002;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00087213;b(0,1)=-0.00017323;b(0,2)= 0.00000404;b(0,3)= 0.00000064;b(0,4)=-0.00000007;b(0,5)= 0.00000000;b(1,0)=-0.00021869;b(1,1)= 0.00009177;b(1,2)=-0.00000589;b(1,3)=-0.00000023;b(1,4)= 0.00000007;b(1,5)= 0.00000000;b(2,0)= 0.00001551;b(2,1)=-0.00001390;b(2,2)= 0.00000173;b(2,3)=-0.00000003;b(2,4)=-0.00000002;b(2,5)= 0.00000000;b(3,0)=-0.00000034;b(3,1)= 0.00000144;b(3,2)=-0.00000032;b(3,3)= 0.00000002;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)=-0.00000009;b(4,1)=-0.00000009;b(4,2)= 0.00000004;b(4,3)=-0.00000001;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000002;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.031695;a(0,1)=-0.002085;a(0,2)=-0.000024;a(0,3)= 0.000007;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.004754;a(1,1)= 0.000915;a(1,2)=-0.000017;a(1,3)=-0.000004;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000300;a(2,1)=-0.000124;a(2,2)= 0.000008;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000014;a(3,1)= 0.000013;a(3,2)=-0.000002;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)=-0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
elseif(x<=45.0)then
	ax(1)=33.0;ax(2)=45.0
	if(y<=1.0)then
		ay(1)=0.0;ay(2)=1.0
		b(0,0)= 0.00024185;b(0,1)= 0.00021311;b(0,2)=-0.00002641;b(0,3)= 0.00000218;b(0,4)=-0.00000014;b(0,5)= 0.00000001;b(1,0)=-0.00007440;b(1,1)=-0.00006555;b(1,2)= 0.00000814;b(1,3)=-0.00000067;b(1,4)= 0.00000004;b(1,5)= 0.00000000;b(2,0)= 0.00000862;b(2,1)= 0.00000759;b(2,2)=-0.00000094;b(2,3)= 0.00000008;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000089;b(3,1)=-0.00000078;b(3,2)= 0.00000010;b(3,3)=-0.00000001;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000009;b(4,1)= 0.00000008;b(4,2)=-0.00000001;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000001;b(5,1)=-0.00000001;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.009211;a(0,1)= 0.008117;a(0,2)=-0.001005;a(0,3)= 0.000083;a(0,4)=-0.000005;a(0,5)= 0.000000;a(1,0)=-0.001425;a(1,1)=-0.001256;a(1,2)= 0.000156;a(1,3)=-0.000013;a(1,4)= 0.000001;a(1,5)= 0.000000;a(2,0)= 0.000110;a(2,1)= 0.000097;a(2,2)=-0.000012;a(2,3)= 0.000001;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000009;a(3,1)=-0.000008;a(3,2)= 0.000001;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000001;a(4,1)= 0.000001;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=2.0)then
		ay(1)=1.0;ay(2)=2.0
		b(0,0)= 0.00051914;b(0,1)= 0.00007790;b(0,2)=-0.00000977;b(0,3)= 0.00000080;b(0,4)=-0.00000005;b(0,5)= 0.00000000;b(1,0)=-0.00015954;b(1,1)=-0.00002382;b(1,2)= 0.00000303;b(1,3)=-0.00000025;b(1,4)= 0.00000002;b(1,5)= 0.00000000;b(2,0)= 0.00001845;b(2,1)= 0.00000274;b(2,2)=-0.00000035;b(2,3)= 0.00000003;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000190;b(3,1)=-0.00000028;b(3,2)= 0.00000004;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000018;b(4,1)= 0.00000003;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000002;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.019784;a(0,1)= 0.002980;a(0,2)=-0.000370;a(0,3)= 0.000031;a(0,4)=-0.000002;a(0,5)= 0.000000;a(1,0)=-0.003059;a(1,1)=-0.000459;a(1,2)= 0.000058;a(1,3)=-0.000005;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000236;a(2,1)= 0.000035;a(2,2)=-0.000004;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000018;a(3,1)=-0.000003;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000001;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=3.0)then
		ay(1)=2.0;ay(2)=3.0
		b(0,0)= 0.00061969;b(0,1)= 0.00002771;b(0,2)=-0.00000365;b(0,3)= 0.00000030;b(0,4)=-0.00000002;b(0,5)= 0.00000000;b(1,0)=-0.00019000;b(1,1)=-0.00000819;b(1,2)= 0.00000115;b(1,3)=-0.00000009;b(1,4)= 0.00000001;b(1,5)= 0.00000000;b(2,0)= 0.00002191;b(2,1)= 0.00000090;b(2,2)=-0.00000014;b(2,3)= 0.00000001;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000224;b(3,1)=-0.00000009;b(3,2)= 0.00000001;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000022;b(4,1)= 0.00000001;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000002;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.023656;a(0,1)= 0.001084;a(0,2)=-0.000137;a(0,3)= 0.000011;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.003652;a(1,1)=-0.000164;a(1,2)= 0.000022;a(1,3)=-0.000002;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000282;a(2,1)= 0.000012;a(2,2)=-0.000002;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000022;a(3,1)=-0.000001;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=4.0)then
		ay(1)=3.0;ay(2)=4.0
		b(0,0)= 0.00065434;b(0,1)= 0.00000880;b(0,2)=-0.00000140;b(0,3)= 0.00000011;b(0,4)=-0.00000001;b(0,5)= 0.00000000;b(1,0)=-0.00019981;b(1,1)=-0.00000219;b(1,2)= 0.00000045;b(1,3)=-0.00000003;b(1,4)= 0.00000000;b(1,5)= 0.00000000;b(2,0)= 0.00002291;b(2,1)= 0.00000017;b(2,2)=-0.00000006;b(2,3)= 0.00000000;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000233;b(3,1)=-0.00000001;b(3,2)= 0.00000001;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000022;b(4,1)= 0.00000000;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000002;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.025051;a(0,1)= 0.000381;a(0,2)=-0.000051;a(0,3)= 0.000004;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.003857;a(1,1)=-0.000052;a(1,2)= 0.000008;a(1,3)=-0.000001;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000296;a(2,1)= 0.000003;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000023;a(3,1)= 0.000000;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=5.333)then
		ay(1)=4.0;ay(2)=5.333
		b(0,0)= 0.00066385;b(0,1)= 0.00000104;b(0,2)=-0.00000089;b(0,3)= 0.00000008;b(0,4)=-0.00000001;b(0,5)= 0.00000000;b(1,0)=-0.00020126;b(1,1)= 0.00000066;b(1,2)= 0.00000031;b(1,3)=-0.00000003;b(1,4)= 0.00000000;b(1,5)= 0.00000000;b(2,0)= 0.00002285;b(2,1)=-0.00000022;b(2,2)=-0.00000004;b(2,3)= 0.00000000;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000230;b(3,1)= 0.00000004;b(3,2)= 0.00000000;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000022;b(4,1)=-0.00000001;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000002;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.025547;a(0,1)= 0.000128;a(0,2)=-0.000030;a(0,3)= 0.000003;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.003915;a(1,1)=-0.000007;a(1,2)= 0.000005;a(1,3)= 0.000000;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000298;a(2,1)=-0.000001;a(2,2)= 0.000000;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000023;a(3,1)= 0.000000;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.00000
	elseif(y<=8.0)then
		ay(1)=5.333;ay(2)=8.0
		b(0,0)= 0.00065629;b(0,1)=-0.00000890;b(0,2)=-0.00000098;b(0,3)= 0.00000010;b(0,4)=-0.00000002;b(0,5)= 0.00000000;b(1,0)=-0.00019521;b(1,1)= 0.00000556;b(1,2)= 0.00000043;b(1,3)=-0.00000004;b(1,4)= 0.00000000;b(1,5)= 0.00000000;b(2,0)= 0.00002161;b(2,1)=-0.00000106;b(2,2)=-0.00000007;b(2,3)= 0.00000000;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000210;b(3,1)= 0.00000016;b(3,2)= 0.00000001;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000019;b(4,1)=-0.00000002;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000002;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.025599;a(0,1)=-0.000077;a(0,2)=-0.000024;a(0,3)= 0.000004;a(0,4)=-0.000001;a(0,5)= 0.000000;a(1,0)=-0.003873;a(1,1)= 0.000050;a(1,2)= 0.000006;a(1,3)=-0.000001;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000290;a(2,1)=-0.000008;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000021;a(3,1)= 0.000001;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000002;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.00000
	elseif(y<=13.0)then
		ay(1)=8.0;ay(2)=13.0
		b(0,0)= 0.00061890;b(0,1)=-0.00002907;b(0,2)=-0.00000138;b(0,3)= 0.00000008;b(0,4)= 0.00000000;b(0,5)= 0.00000000;b(1,0)=-0.00017377;b(1,1)= 0.00001615;b(1,2)= 0.00000061;b(1,3)=-0.00000006;b(1,4)= 0.00000000;b(1,5)= 0.00000000;b(2,0)= 0.00001776;b(2,1)=-0.00000281;b(2,2)=-0.00000007;b(2,3)= 0.00000001;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000156;b(3,1)= 0.00000039;b(3,2)= 0.00000000;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000012;b(4,1)=-0.00000005;b(4,2)= 0.00000000;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)=-0.00000001;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.025141;a(0,1)=-0.000385;a(0,2)=-0.000023;a(0,3)= 0.000001;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.003660;a(1,1)= 0.000166;a(1,2)= 0.000008;a(1,3)= 0.000000;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000258;a(2,1)=-0.000024;a(2,2)=-0.000001;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000018;a(3,1)= 0.000003;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000001;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	elseif(y<=21.0)then
		ay(1)=13.0;ay(2)=21.0
		b(0,0)= 0.00052984;b(0,1)=-0.00005978;b(0,2)=-0.00000085;b(0,3)= 0.00000023;b(0,4)=-0.00000001;b(0,5)= 0.00000000;b(1,0)=-0.00012764;b(1,1)= 0.00002934;b(1,2)=-0.00000024;b(1,3)=-0.00000014;b(1,4)= 0.00000001;b(1,5)= 0.00000000;b(2,0)= 0.00001039;b(2,1)=-0.00000436;b(2,2)= 0.00000017;b(2,3)= 0.00000002;b(2,4)= 0.00000000;b(2,5)= 0.00000000;b(3,0)=-0.00000064;b(3,1)= 0.00000049;b(3,2)=-0.00000004;b(3,3)= 0.00000000;b(3,4)= 0.00000000;b(3,5)= 0.00000000;b(4,0)= 0.00000003;b(4,1)=-0.00000005;b(4,2)= 0.00000001;b(4,3)= 0.00000000;b(4,4)= 0.00000000;b(4,5)= 0.00000000;b(5,0)= 0.00000000;b(5,1)= 0.00000000;b(5,2)= 0.00000000;b(5,3)= 0.00000000;b(5,4)= 0.00000000;b(5,5)= 0.00000000
		a(0,0)= 0.023878;a(0,1)=-0.000889;a(0,2)=-0.000030;a(0,3)= 0.000002;a(0,4)= 0.000000;a(0,5)= 0.000000;a(1,0)=-0.003148;a(1,1)= 0.000346;a(1,2)= 0.000006;a(1,3)=-0.000001;a(1,4)= 0.000000;a(1,5)= 0.000000;a(2,0)= 0.000191;a(2,1)=-0.000043;a(2,2)= 0.000000;a(2,3)= 0.000000;a(2,4)= 0.000000;a(2,5)= 0.000000;a(3,0)=-0.000010;a(3,1)= 0.000004;a(3,2)= 0.000000;a(3,3)= 0.000000;a(3,4)= 0.000000;a(3,5)= 0.000000;a(4,0)= 0.000000;a(4,1)= 0.000000;a(4,2)= 0.000000;a(4,3)= 0.000000;a(4,4)= 0.000000;a(4,5)= 0.000000;a(5,0)= 0.000000;a(5,1)= 0.000000;a(5,2)= 0.000000;a(5,3)= 0.000000;a(5,4)= 0.000000;a(5,5)= 0.000000
	else
		write(*,*)'error'
	endif
else
	write(*,*)'error'
endif
end subroutine
!---------------------------------------------------


!----------------------------------------------------------------------
END MODULE GREEN_MOD