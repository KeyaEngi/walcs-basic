!****************************************************************************************************
!模块功能：获得高斯积分的积分节点OrderNum多项式的阶数和权函数
!****************************************************************************************************
MODULE GAUSSIN_MOD

USE ENVIRONMENT_MOD,ONLY:PI

IMPLICIT NONE
PRIVATE
PUBLIC::GaussIn
CONTAINS 
!--------------------------------------------------------------------
!获得高斯积分的积分节点和权函数
subroutine GaussIn(OrderNum,xk,wk,flag)
implicit none 
integer(4),intent(in)::OrderNum,flag
real(8),intent(out),dimension(1:OrderNum)::xk,wk
real(8),allocatable,dimension(:)::polycoef !多项式的系数
allocate(polycoef(1:OrderNum+1))
call polynome(OrderNum,polycoef,flag)
call polyroot(OrderNum,polycoef,xk)
call polywk(OrderNum,polycoef,xk,wk,flag)
deallocate(polycoef)
end subroutine GaussIn


!*********************************************************************
!本程序计算正交多项式的系数
!输入参数:OrderNum为多项式的阶数,flag为控制参数:
!         1--------laguerre多项式
!         2--------legendre多项式
!         3--------hermite多项式
!输出参数:polycoef(1:OrderNum+1)为多项式的系数(升幂排列)
!应用范围:若要计算其它正交多项式的系数,需扩充 select case模块(给出三项递推关系)
!程序编制：张海彬     时间：1999年12月11日
!****************************************************************************
subroutine polynome(OrderNum,polycoef,flag)
implicit none
integer,intent(in)::OrderNum	!多项式的阶数
integer,intent(in)::flag        !控制参数
real(8),intent(out),dimension(1:OrderNum+1)::polycoef !多项式的系数
integer::i
real(8),allocatable,dimension(:)::polycoef0,polycoef1,ltem
allocate(polycoef0(1:OrderNum+1),polycoef1(1:OrderNum+1),ltem(1:OrderNum+1))
polycoef0=0.0;polycoef1=0.0;polycoef=0.0;ltem=0.0	  !清零
select case(flag)
case(1)	 !计算laguerre多项式系数
	polycoef0(1)=1.0;polycoef1(1)=1.0		!赋前两项系数
	if(OrderNum/=0)polycoef1(2)=-1.0
	select case(OrderNum)
	case(0)
		polycoef=polycoef0
	case(1)
		polycoef=polycoef1
	case default	!递推过程
		do i=1,OrderNum-1
			ltem(1)=0
			ltem(2:OrderNum+1)=polycoef1(1:OrderNum) !x项系数
			polycoef=(1+2*i)*polycoef1-ltem-i**2*polycoef0
			polycoef0=polycoef1
			polycoef1=polycoef
		enddo
	end select
case(2)	 !计算legendre多项式系数
	polycoef0(1)=1.0;polycoef1(1)=0.0		!赋前两项系数
	if(OrderNum/=0)polycoef1(2)=1.0
	select case(OrderNum)
	case(0)
		polycoef=polycoef0
	case(1)
		polycoef=polycoef1
	case default	!递推过程
		do i=1,OrderNum-1
			ltem(1)=0
			ltem(2:OrderNum+1)=polycoef1(1:OrderNum) !x项系数
			polycoef=((1+2*i)*ltem-i*polycoef0)/(i+1)
			polycoef0=polycoef1
			polycoef1=polycoef
		enddo
	end select
case(3)	 !计算hermite多项式系数
	polycoef0(1)=1.0;polycoef1(1)=0.0		!赋前两项系数
	if(OrderNum/=0)polycoef1(2)=2.0
	select case(OrderNum)
	case(0)
		polycoef=polycoef0
	case(1)
		polycoef=polycoef1
	case default	!递推过程
		do i=1,OrderNum-1
			ltem(1)=0
			ltem(2:OrderNum+1)=polycoef1(1:OrderNum) !x项系数
			polycoef=2*ltem-2*i*polycoef0
			polycoef0=polycoef1
			polycoef1=polycoef
		enddo
	end select
end select
deallocate(polycoef0,polycoef1,ltem)
end subroutine polynome


!***************************************************************************
!本程序用劈因子(林士谔方法)迭代法计算正交多项式的零点
!输入参数：OrderNum为多项式的阶数，polycoef(1:OrderNum+1)为多项式的系数(升幂排列)
!输出参数：root(1:OrderNum)为多项式的零点
!应用范围：本程序可以计算所有无复数零点多项式的零点(正交多项式满足此要求)，若要计算多项式
!          的复数零点，需要对程序作出相应的修改。
!程序编制：张海彬     时间：1999年12月11日
!***************************************************************************            
subroutine polyroot(OrderNum,polycoef,root)
implicit none
integer,intent(in)::OrderNum  !多项式的阶数
real(8),intent(in),dimension(1:OrderNum+1)::polycoef !多项式的系数(升幂排列)
real(8),intent(out),dimension(1:OrderNum)::root  !多项式的零点
real(8),allocatable,dimension(:)::a	  !劈因子前多项式的系数
real(8),allocatable,dimension(:)::b	  !劈因子后多项式的系数
integer::aOrderNum,bOrderNum		  !劈因子前后多项式的阶数
integer::i,k,i0
real(8)::p0,q0,p1,q1,deta,esp1	!劈因子系数、Δ 和误差控制
allocate(a(1:OrderNum+1),b(1:OrderNum))
root=0.0							!赋初值
esp1=1.0e-10
aOrderNum=OrderNum
i0=0	
a=polycoef							
do
!判断零根
	if(abs(a(1))<=esp1)then
		i0=i0+1
		a(1:aOrderNum)=a(2:aOrderNum+1)
		aOrderNum=aOrderNum-1
		root(i0)=0.0
	else
		exit
	endif
enddo
bOrderNum=aOrderNum-2
a=a/a(aOrderNum+1)			!使多项式最高次项系数为1
p0=3;q0=1						!劈因子迭代初值
!p0=0;q0=-1						!劈因子迭代初值
i=0
do
! 劈因子求零点
	i=i+1						!劈因子计数器
!-----------------------------------------------------------
    !循环出口
	if(aOrderNum<=0)then
		exit
	elseif(aOrderNum==1)then
		root(2*i-1+i0)=-a(1)/a(2)
		exit
	elseif(aOrderNum==2)then
		deta=a(2)**2-4*a(1)*a(3)
		if(deta>=0.0)then
			root(2*i-1+i0)=(-a(2)-sqrt(deta))/(2*a(3))
			root(2*i+i0)=(-a(2)+sqrt(deta))/(2*a(3))		 
		endif
		exit
	endif
!-----------------------------------------------------
	b(aOrderNum)=0;b(aOrderNum-1)=1	 !赋劈因子后多项式最高次项的系数
	do
	!迭代劈因子系数
  		do k=bOrderNum,1,-1
		!计算劈因子后多项式其余项的系数
			b(k)=a(k+2)-(p0*b(k+1)+q0*b(k+2))  
		enddo
		p1=(a(2)-q0*b(2))/b(1)				   !修正的劈因子系数
		q1=a(1)/b(1)
		if(abs(p1-p0)<esp1.and.abs(q1-q0)<esp1)then	!劈因子系数精度控制
			deta=p1**2-4.0*q1
			if(deta>=0.0)then
				root(2*i-1)=(-p1-sqrt(deta))/2.0		!求劈因子零点
				root(2*i)=(-p1+sqrt(deta))/2.0
			endif
			a(1:bOrderNum+1)=b(1:bOrderNum+1)		!进行下一次劈因子
			aOrderNum=aOrderNum-2
			bOrderNum=aOrderNum-2
			exit					!迭代出口
		else
			p0=p1;q0=q1
		endif
	enddo
enddo
deallocate(a,b)
end subroutine polyroot


!***********************************************************************************
!本程序计算高斯型积分中正交多项式的权系数
!输入参数：OrderNum为正交多项式的阶数，polycoef(1:OrderNum+1)正交多项式的系数(升幂排列)
!          xk(1:OrderNum)为正交多项式的零点
!          flag为控制参数:
!                         1--------laguerre多项式
!                         2--------legendre多项式
!                         3--------hermite多项式
!输出参数：wk(1:OrderNum)高斯型积分相应积分节点的权系数
!应用范围：若要计算其它正交多项式的权系数,需扩充 select case模块(给出权系数表达式)
!程序编制：张海彬     时间：1999年12月11日
!***********************************************************************************
subroutine polywk(OrderNum,polycoef,xk,wk,flag)
implicit none
integer,intent(in)::OrderNum
real(8),intent(in),dimension(1:OrderNum+1)::polycoef
real(8),intent(in),dimension(1:OrderNum)::xk
integer,intent(in)::flag
real(8),intent(out),dimension(1:OrderNum)::wk
real(8),allocatable,dimension(:)::ld !正交多项式(OrderNum>0)导数或前一阶多项式系数
real(8)::ldvalue				  !导数或前一阶多项式值
integer::i,j
!阶乘计算模块接口
! interface
! 	recursive function fac(n) result(fac_result)
! 	integer,intent(in)::n
! 	real(8)::fac_result
! 	end function fac
! end interface
allocate(ld(1:OrderNum))

select case(flag)
case(1)	 !求Gauss-Laguerre积分权系数
	!求正交多项式(OrderNum>0)导数多项式系数
	do i=1,OrderNum
		ld(i)=polycoef(i+1)*i
	enddo
	do i=1,OrderNum
		ldvalue=0.0
		do j=1,OrderNum
			ldvalue=ldvalue+ld(j)*(xk(i)**(j-1)) !计算导数值
		enddo
		wk(i)=fac(OrderNum)**2/(xk(i)*(ldvalue**2))
	enddo
case(2)  !求Gauss-Legendre积分权系数
	!计算OrderNum-1阶多项式系数
    call polynome(OrderNum-1,ld,flag)
	do i=1,OrderNum
		ldvalue=0.0
		do j=1,OrderNum
			ldvalue=ldvalue+ld(j)*(xk(i)**(j-1)) !计算OrderNum-1阶多项式值
		enddo
		wk(i)=(2.0*(1.0-xk(i)**2))/((OrderNum*ldvalue)**2)
	enddo
case(3)  !求Gauss-Hermite积分权系数
	!求正交多项式(OrderNum>0)导数多项式系数
	do i=1,OrderNum
		ld(i)=polycoef(i+1)*i
	enddo
	do i=1,OrderNum
		ldvalue=0.0
		do j=1,OrderNum
			ldvalue=ldvalue+ld(j)*(xk(i)**(j-1)) !计算导数值
		enddo
		wk(i)=2**(OrderNum+1)*fac(OrderNum)*sqrt(PI)/(ldvalue**2)
	enddo
end select
deallocate(ld)
end subroutine polywk
!-----------------------------------------------------------------------------------

recursive function fac(n) result(fac_result)
!本子程序用来计算阶乘
integer,intent(in)::n
real(8)::fac_result	  !由于存在整数越界的可能，计算结果保存为实型
if(n==0)then
	fac_result=1
else
	fac_result=n*fac(n-1)
endif
end function fac



END MODULE GAUSSIN_MOD
!**************************************************************************************