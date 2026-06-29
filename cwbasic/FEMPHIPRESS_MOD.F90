MODULE FEMPHIPRESS_MOD

USE MESH_MOD,ONLY:XAV,NPANEL,EA,Factor,TarNum,Space
USE ENVIRONMENT_MOD,ONLY:H

USE GREEN_MOD
USE GAUSSIN_MOD
USE PVPOINT_MOD,ONLY:NPPOINT,PRESSPOINT
USE CAL_MOD,ONLY:WATERDEPTH


IMPLICIT NONE
PRIVATE
PUBLIC::GetPressPhi


CONTAINS
!**************************************************************************
!程序功能：压力计算点速度势梯度及速度势
!输出参数：Phi_press(1:7,1:NPPoint,1:2,1:4)为各压力计算点速度势梯度及速度势
!调用子程序：GR_P,coef_P              
!程序编制：孙葳  
!***************************************************************************
subroutine GetPressPhi(wkjNum,wkj,Pnu,S_press,SSTRENGTH,Phi_press)
implicit none
integer(4),intent(in)::wkjNum
real(8),intent(inout)::Pnu
real(8),intent(IN),dimension(1:NPPoint,1:NPanel,1:4,1:Factor)::S_PRESS
real(8),intent(in),dimension(1:wkjNum)::wkj
real(8),intent(in),dimension(1:2,1:7,1:Factor,1:NPANEL)::SSTRENGTH
real(8),intent(out),dimension(1:2,1:2,1:7,1:NPPoint)::Phi_press
integer(4)::k,L,i,j
real(8),allocatable,dimension(:,:,:,:)::GcS_press,GsS_press

allocate(GcS_press(1:NPPoint,1:NPanel,1:2,1:Factor),GsS_press(1:NPPoint,1:NPanel,1:2,1:Factor))

Phi_press=0.
call coef_P(S_press,wkjNum,wkj,Pnu,GcS_press,GsS_press)  !计算格林函数及其偏导数
   
!计算速度势

do L=1,7      
   do k=1,2
       do i=1,Factor  !两个对称面的修改  edited by 王川 2012.07.30
        Phi_press(1,k,L,1:NPpoint)=Phi_press(1,k,L,1:NPpoint)-matmul(GsS_press(1:NPpoint,1:Npanel,k,i),SSTRENGTH(2,L,i,1:Npanel))+matmul(GcS_press(1:NPpoint,1:Npanel,k,i),SSTRENGTH(1,L,i,1:Npanel))
        Phi_press(2,k,L,1:NPpoint)=Phi_press(2,k,L,1:NPpoint)+matmul(GcS_press(1:NPpoint,1:Npanel,k,i),SSTRENGTH(2,L,i,1:Npanel))+matmul(GsS_press(1:NPpoint,1:Npanel,k,i),SSTRENGTH(1,L,i,1:Npanel))
      end do
   enddo
enddo          
deallocate(GcS_press,GsS_press)
end subroutine GetPressPhi
!-----------------------------------------------------------------------------------------------




!**************************************************************************************************
!程序功能：计算格林函数G及其梯度、法向导数沿面元的积分(对压力计算点)
!         s_press为压力计算点在各面元下的诱导速度及单层势
!         NPoint为压力计算点数目，PressPoint为压力计算点坐标  
!         OrderNum为高斯积分节点数目
!输出参数：GcS_press(i,j,1:4),GsS_press(i,j,1:4)分别为格林函数梯度及其本身沿面元的积分的实部和虚部
!调用子程序：FG1,FG2,Flags
!程序编制：张海彬    时间：2002年4月16日        
!**************************************************************************************************
subroutine coef_P(S_press,wkjNum,wkj,Pnu,GcS_press,GsS_press)
implicit none
integer(4),intent(in)::wkjNum
real(8),intent(in),dimension(1:NPPoint,1:NPanel,1:4,1:Factor)::S_press
real(8),intent(in),dimension(1:wkjNum)::wkj
real(8),intent(inout)::Pnu
real(8),intent(out),dimension(1:NPPoint,1:NPanel,1:2,1:Factor)::GcS_press,GsS_press
real(8),dimension(1:4)::Sc,Ss !,Sc1,Ss1
real(8),dimension(0:2)::Gc,Gs !,Gc1,Gs1
real(8),allocatable,dimension(:)::xk,wk
real(8),allocatable,dimension(:,:)::x00
integer(4)::i,j,ij,count,L
real(8)::x,y,R,pz,qz,tem,eps
INTEGER(4)::OrderNum

OrderNum=15
L=0
! PI=4.0*atan(1.0)
eps=1.0e-6
allocate(xk(1:OrderNum),wk(1:OrderNum),x00(1:3,1:Npanel))
if((trim(adjustl(WATERDEPTH))=="FINITE")) then !有限水深
	call GaussIn(OrderNum,xk,wk,1) !获得高斯积分的积分节点和权函数
endif
!采用对称性计算,需要进行四部分矩阵的修复计算。edited by 王川
!A矩阵生成

do count=1,TarNum,Space
select   case(count)
               case(1)
                      x00(1,:)=xav(1,:)
                      x00(2,:)=xav(2,:)
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
do i=1,NPPoint
 
	do j=1,NPanel
			Sc=0.0;Ss=0.0
			x=PressPoint(1,i)-x00(1,j)
			y=PressPoint(2,i)-x00(2,j)
			R=sqrt(x*x+y*y)
			pz=PressPoint(3,i)
			qz=xav(3,j)
			tem=sqrt(R**2+(pz-qz)**2)
		!---------------------------------------------------------------------
			if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
				ij=1
				if(R/H>0.5)ij=2
				if(tem<eps)ij=1
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

			Sc(4)=Gc(0)
			Ss(4)=Gs(0)

			if(abs(R)>1e-6)then

				Sc(1)=Gc(1)*x/R
				Ss(1)=Gs(1)*x/R

			endif
			
		GcS_press(i,j,1,L)=Sc(4)*ea(j)+S_press(i,j,4,L)
		GsS_press(i,j,1,L)=Ss(4)*ea(j)
		GcS_press(i,j,2,L)=Sc(1)*ea(j)+S_press(i,j,1,L)
		GsS_press(i,j,2,L)=Ss(1)*ea(j)

		!------------------------------------------
	enddo
enddo
end do

deallocate(xk,wk,x00)
end subroutine coef_P



END MODULE FEMPHIPRESS_MOD