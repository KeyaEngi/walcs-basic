MODULE SECTLOAD_MOD

USE MESH_MOD,ONLY:NB,HRM,SECTHRM ,sectEA ,wetea,Factor,TarNum,Space
USE ENVIRONMENT_MOD,ONLY:ROU,G0,ZG,MLN1,XML,XG,SML,NBSECT,SECTMASSMATRIX,MASSMATRIX,SECTPSN,SECTCOG,PI,SECTTYPE,ZSC
USE PRINT_MOD,ONLY:ANSMETHOD
USE WAVEFORCE_MOD,ONLY:SECTLOAD 

IMPLICIT NONE
PRIVATE
PUBLIC::SHIPSECTLOAD, ShipSectLoad_lt

CONTAINS
!-----------------------------------------------------------


!***************************************************************************************************
!程序功能：求解船舶各站剖面载荷
!输出参数：sectLoad(1:2,1:6,1:NBSECT)为船舶各站剖面载荷的实部和虚部(N,SFy,SFz,TM,BMy,BMz)
!***************************************************************************************************
subroutine ShipSectLoad	(Beta,e0,ea0,pall,xav0,OMEI0,omei,omeI1,WAVEPHASE,AMPI,motion,M4x,IP)
implicit	none

real(8)	,	intent(in)	::	OMEI0,omei,Beta,WAVEPHASE,AMPI,omeI1
real(8)	,	intent(in)	,	dimension(1:3,1:NB)			::xav0
real(8)	,	intent(in)	,	dimension(1:3,1:3,1:NB)		::e0
real(8)	,	intent(in)	,	dimension(1:NB)				::ea0
real(8)	,	intent(in)	,	dimension(1:2,1:Factor,1:NB)			::pall
real(8)	,	intent(in)	,	dimension(1:2,1:6)			::motion
real(8)	,	allocatable	,	dimension(:,:,:)				::n
real(8)	,	allocatable	,	dimension(:,:)				::X00
real(8)	,	dimension(1:2,1:6)	::	temp,tempS
real(8)	,	dimension(1:3)		::	temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8 !高斯公式值
real(8)	,	dimension(1:8)		::	adjust_m
real(8)	,	dimension(1:2)		::	x10
integer(4)  ::  IP
integer(4)	::	i,j,k,L,count
real(8)::Aall,SectA(1:NBSECT),MS(1:2,1:NBSECT),M4x(1:2)
xMl(:)=xMl(:)-XG
ALLOCATE(SECTLOAD(1:2,1:6,1:NBSECT),N(1:6,1:Factor,1:NB),x00(1:3,1:NB))

do i=1,NB
	do j=1,3
		n(J,1,I)=e0(j,3,i)
	enddo
	n(4,1,I)=xav0(2,i)*e0(3,3,i)-(xav0(3,i)-ZG)*e0(2,3,i)
	n(5,1,I)=(xav0(3,i)-ZG)*e0(1,3,i)-xav0(1,i)*e0(3,3,i)
	n(6,1,I)=xav0(1,i)*e0(2,3,i)-xav0(2,i)*e0(1,3,i)
enddo
!初始化
L=0
!重量调整-----------------------------------------------
temp1=0.0;	temp2=0.0;	temp3=0.0;	temp4=0.0
temp5=0.0;	temp6=0.0;	temp7=0.0;	temp8=0.0
do count=1,TarNum,Space
    L=L+1
    select   case(count)
       case(1)
           x00(1,:)=xav0(1,:)
           x00(2,:)=xav0(2,:)
           n(1,L,1:NB)=n(1,1,1:NB)
           n(2,L,1:NB)=n(2,1,1:NB)
           n(3,L,1:NB)=n(3,1,1:NB)
           n(4,L,1:NB)=n(4,1,1:NB)
           n(5,L,1:NB)=n(5,1,1:NB)
           n(6,L,1:NB)=n(6,1,1:NB)             
        case(2)
           x00(1,:)=-xav0(1,:)
           x00(2,:)=xav0(2,:)
           n(1,L,1:NB)=-n(1,1,1:NB)
           n(2,L,1:NB)=n(2,1,1:NB)
           n(3,L,1:NB)=n(3,1,1:NB)
           n(4,L,1:NB)=n(4,1,1:NB)
           n(5,L,1:NB)=-n(5,1,1:NB)
           n(6,L,1:NB)=-n(6,1,1:NB)
       case(3)
           x00(1,:)=-xav0(1,:)
           x00(2,:)=-xav0(2,:)
           n(1,L,1:NB)=-n(1,1,1:NB)
           n(2,L,1:NB)=-n(2,1,1:NB)
           n(3,L,1:NB)=n(3,1,1:NB)
           n(4,L,1:NB)=-n(4,1,1:NB)
           n(5,L,1:NB)=-n(5,1,1:NB)
           n(6,L,1:NB)=n(6,1,1:NB)
        case(4)
           x00(1,:)=xav0(1,:)
           x00(2,:)=-xav0(2,:)
           n(1,L,1:NB)=n(1,1,1:NB)
           n(2,L,1:NB)=-n(2,1,1:NB)
           n(3,L,1:NB)=n(3,1,1:NB)
           n(4,L,1:NB)=-n(4,1,1:NB)
           n(5,L,1:NB)=n(5,1,1:NB)
           n(6,L,1:NB)=-n(6,1,1:NB)
    end select
    do i=	1,	NB
	    temp1(:)	=temp1(:)+n(1:3,L,i)*ea0(i)
	    temp2(:)	=temp2(:)+x00(1,i)*n(1:3,L,i)*ea0(i)
	    temp3(:)	=temp3(:)+x00(2,i)*n(1:3,L,i)*ea0(i)
	    temp4(:)	=temp4(:)+xav0(3,i)*n(1:3,L,i)*ea0(i)
	    temp5(:)	=temp5(:)+n(4:6,L,i)*ea0(i)
	    temp6(:)	=temp6(:)+x00(1,i)*n(4:6,L,i)*ea0(i)
	    temp7(:)	=temp7(:)+x00(2,i)*n(4:6,L,i)*ea0(i)
	    temp8(:)	=temp8(:)+xav0(3,i)*n(4:6,L,i)*ea0(i)
    end	do
end do

temp=0.0
sectload = 0.0
!-------------------------------------------------
Aall=0.0; SECTA=0.0
do j=1,NB
	Aall=Aall+ea0(j)
end do
Aall=Factor*Aall
!----------------------------------------------------------
!-----程序调试：法向量和面积输出---------------------------
open(55,file='AreaN.dat')

do  i=1,    NB

    write(55,'(7(2x,f10.5))')(n(j,1,i),j=1,6),ea0(i)
    
enddo
!------------------------------------------------------------
!--------------------------------------------------
do	i=	1,	nbsect

    L=0
	temp(1,1:6)	=0.0
	temp(2,1:6)	=0.0
	
	do count=1,TarNum,Space
	    L=L+1
	    do	j=	1,	NB
	        select case(count)
                case(1)
                    x10(1)=xav0(1,j)
                case(2)
                    x10(1)=-xav0(1,j)      
                case(3)
                    x10(1)=-xav0(1,j)
                case(4)
                    x10(1)=xav0(1,j)
            end select
        
		    if	(x10(1)<sectpsn(1,i))	then
		
			    temp(1,1:6)	=temp(1,1:6)-pall(1,L,j)*n(1:6,L,j)*ea0(j)
			    temp(2,1:6)	=temp(2,1:6)-pall(2,L,j)*n(1:6,L,j)*ea0(j)
			    SectA(i)    =SectA(i)+ea0(j)
			
		    end	if
	    end do
    end	do
	SectLoad(1,1:6,i)	=temp(1,1:6)-omei*omei*Matmul(SECTMASSMATRIX(1:6,1:6,i),motion(1,1:6))
	SectLoad(2,1:6,i)	=temp(2,1:6)-omei*omei*Matmul(SECTMASSMATRIX(1:6,1:6,i),motion(2,1:6))
	!去出重力分量
	adjust_m(1)=	-rou*g0*temp2(1)*SECTMASSMATRIX(1,1,I)/MASSMATRIX(1,1) !沿x
	adjust_m(2)=	-rou*g0*temp3(2)*SECTMASSMATRIX(1,1,I)/MASSMATRIX(1,1)         !沿y
	!沿z
	adjust_m(3)=	rou*g0*(-motion(1,3)*(temp1(3)-HRM(3,3)/rou/g0)+	&
					motion(1,5)*(temp2(3)+HRM(3,5)/rou/g0))
	adjust_m(4)=	rou*g0*(-motion(2,3)*(temp1(3)-HRM(3,3)/rou/g0)+	&
					motion(2,5)*(temp2(3)+HRM(3,5)/rou/g0))
	adjust_m(3:4)=	adjust_m(3:4)*SECTMASSMATRIX(1,1,I)/MASSMATRIX(1,1)  !sectm(ns+2,1,1)
	!绕x
	adjust_m(5)=	-rou*g0*(temp7(1)-HRM(4,4)/rou/g0)*SECTMASSMATRIX(1,1,I)/MASSMATRIX(1,1) !sectm(ns+2,1,1) 
	!绕y
	adjust_m(6)=	rou*g0*(-motion(1,3)*(temp5(2)-HRM(3,5)/rou/g0)+	&
					motion(1,5)*(temp6(2)+HRM(5,5)/rou/g0))
	adjust_m(7)=	rou*g0*(-motion(2,3)*(temp5(2)-HRM(3,5)/rou/g0)+&
					motion(2,5)*(temp6(2)+HRM(5,5)/rou/g0))
	adjust_m(6:7)=	adjust_m(6:7)*SECTMASSMATRIX(1,1,I)/MASSMATRIX(1,1)  !sectm(ns+2,1,1)
	!绕z
	adjust_m(8)=	-rou*g0*temp7(3)*SECTMASSMATRIX(1,1,I)/MASSMATRIX(1,1) !sectm(ns+2,1,1) 

	SectLoad(1:2,1,i)=	SectLoad(1:2,1,i)-adjust_m(1)*motion(1:2,5)  
	SectLoad(1:2,2,i)=	SectLoad(1:2,2,i)+adjust_m(2)*motion(1:2,4)
	SectLoad(1:2,3,i)=	SectLoad(1:2,3,i)+adjust_m(3:4)
	SectLoad(1:2,4,i)=	SectLoad(1:2,4,i)+adjust_m(5)*motion(1:2,4) !是否有问题？符号应该是正还是负？
	SectLoad(1:2,5,i)=	SectLoad(1:2,5,i)+adjust_m(6:7)
	SectLoad(1:2,6,i)=	SectLoad(1:2,6,i)+adjust_m(8)*motion(1:2,4)+&
						adjust_m(2)*motion(1:2,4)*SECTCOG(1,I)     

	!考虑系泊力
	do	j=1,	mlN1
		if	(SECTPSN(1,I)>XML(j))	then
			SectLoad(1,1:6,i)	=SectLoad(1,1:6,i)+Matmul(SML(1:6,1:6,j),motion(1,1:6))
			SectLoad(2,1:6,i)	=SectLoad(2,1:6,i)+Matmul(SML(1:6,1:6,j),motion(2,1:6))
		end	if
	end	do
	!消除剪力引起的弯矩
	SectLoad(1:2,4,i)=	SectLoad(1:2,4,i)-SectLoad(1:2,2,i)*(sectpsn(3,i)-ZSC(i)) !修改剪力修正前符号 By LiHui 2014.4.18
	SectLoad(1:2,5,i)=	SectLoad(1:2,5,i)+SectLoad(1:2,3,i)*sectpsn(1,i) !xs(i)
	SectLoad(1:2,6,i)=	SectLoad(1:2,6,i)-SectLoad(1:2,2,i)*sectpsn(1,i) !xs(i)
	
	!横摇阻尼力矩修正
	Ms(1:2,I)=M4x(1:2)*secta(i)/Aall	
	sectload(1:2,4,i)=sectload(1:2,4,i)+Ms(1:2,I) !正负号
	
end	do 

if(abs(sin(Beta))<0.001)then
	SectLoad(:,2,:)=0.0;SectLoad(:,4,:)=0.0;SectLoad(:,6,:)=0.0
endif

deallocate(n)

if(IP>0) then
    call dsls_putouts	(BETA,omeI0,omeI,omeI1,WAVEPHASE,AMPI)
endif
!----------------------------------------------------------------------------
DEALLOCATE(SECTLOAD)

end subroutine ShipSectLoad


!***************************************************************************************************
!程序功能：求解浮体纵向（横向）剖面载荷
!输出参数：sectLoad(1:2,1:6,1:NBSECT)为船舶各站剖面载荷的实部和虚部(N,SFy,SFz,TM,BMy,BMz)
!***************************************************************************************************
subroutine ShipSectLoad_lt	(Beta,OMEI0,omei,omeI1,WAVEPHASE,AMPI,motion,bv44,sectf,IP)
implicit	none
real(8)	,	intent(in)	::	OMEI0,omei,Beta,WAVEPHASE,AMPI,omeI1
real(8)	,	intent(in)	,	dimension(1:2,1:6)			::motion
real(8),	intent(in)::bv44
real(8),	intent(in)::sectf(1:2,1:6,1:nbsect,1:9) 	
real(8)::A(1:6,1:6),B(1:6,1:6) ! 分段附加质量和阻尼
real(8):: CWMA(1:6,1:6),WB(1:6,1:6)!运动方程系数
real(8):: EF(1:2,1:6,1:nbsect)      ! 波浪激励力
real(8):: FG(1:2,1:6)               !未去除剪力贡献
real(8):: FC(1:2,1:6,1:7,1:nbsect)  !求解剖面载荷方程系数
integer(4):: IP
integer(4)	::	i,j,k,II,JJ,KK
FC=0.0
ALLOCATE(SECTLOAD(1:2,1:6,1:NBSECT))

do k=1,nbsect
    DO I=1,6
        DO J=1,6
            A(I,J)=sectf(1,I,K,J)
            B(I,J)=sectf(2,I,K,J)
            if(i==4.and.j==4) B(4,4)=B(4,4)+BV44/WETEA*SECTEA(K)
            CWMA(I,J)=SECTHRM(I,J,K)-OMEI*OMEI*(SECTMASSMATRIX(I,J,K)+A(I,J))
            WB(I,J)=OMEI*B(I,J)  

        ENDDO
    ENDDO

    FG(1,1:6)=Matmul(CWMA(1:6,1:6),motion(1,1:6))-Matmul(WB(1:6,1:6),motion(2,1:6))-sectf(1,1:6,K,9)
    FG(2,1:6)=Matmul(CWMA(1:6,1:6),motion(2,1:6))+Matmul(WB(1:6,1:6),motion(1,1:6))-sectf(2,1:6,K,9)
    sectload(:,1,K)=FG(:,1);sectload(:,2,K)=FG(:,2);sectload(:,3,K)=FG(:,3)
    sectload(:,4,k)=FG(:,4)+FG(:,2)*(SECTPSN(3,K)-ZSC(k))-FG(:,3)*SECTPSN(2,K)  !剪力的调整
    sectload(:,5,k)=FG(:,5)+FG(:,3)*SECTPSN(1,K)-FG(:,1)*(SECTPSN(3,K)-ZSC(k))
    !write(*,*)SECTPSN(3,K)-ZSC(k)
    sectload(:,6,k)=FG(:,6)+FG(:,1)*SECTPSN(2,K)-FG(:,2)*SECTPSN(1,K)
!    IF(SECTTYPE(K)==1) SECTLOAD(:,4,K)=SECTLOAD(:,4,K)-SECTLOAD(:,2,K)*(SECTPSN(3,K)-ZSC(K)) !横剖面扭矩调整到扭心
!    IF(SECTTYPE(K)==2) SECTLOAD(:,5,K)=SECTLOAD(:,5,K)+SECTLOAD(:,1,K)*(SECTPSN(3,K)-ZSC(K)) !纵剖面扭矩调整到扭心          
enddo
if(IP>0)   then
    call dsls_putouts	(BETA,omeI0,omeI,omeI1,WAVEPHASE,AMPI)
endif

DEALLOCATE(SECTLOAD)

end subroutine ShipSectLoad_LT

!*************************************************************************************
!输出船舶各站剖面载荷结果
!*************************************************************************************
subroutine dsls_putouts(BETA,omeI0,omeI,omeI1,WAVEPHASE,AMPI)
implicit none
! integer(4),intent(in)::NP
! real(8),intent(in)::sl,sb,rou,g0
real(8),intent(in)::omeI0,omeI,omeI1,wavephase,AMPI,BETA
! real(8),intent(in),dimension(1:2,1:6,1:Nbsect)::SectLoad
! real(8),intent(in),dimension(0:NP+1,1:6,1:2)::MidSL
real(8),dimension(1:6,1:2)::temp
integer(4)::j,k,kk
! real(8)::PI 
! PI=4.0*atan(1.0)

	if(TRIM(ADJUSTL(ansMethod))=="YES")	then
		write(132,'(a,f5.1,a,f6.3)')"# HEAD= ",Beta*180.0/PI,"    Ome0= ",omei0
		write(132,'(a)')"#  ome0  ome"
		write(132,'(2f7.3)')omei0,omeI1
		write(132,'(a,6x,a,2x,3x,a,1x,2(3x,a,2x,3x,a,1x),5x,a,2x,3x,a,1x,2(3x,a,2x,3x,a,1x))')"#",&
						    "  F.X.","PHA"," F.Y. ","PHA"," F.Z. ","PHA","M.X.","PHA","  M.Y.","PHA","  M.Z.","PHA"
		do j=1,Nbsect
			do k=1,6
				call comp_to_ampha(SectLoad(1,k,j),SectLoad(2,k,j),temp(k,1),temp(k,2))
			enddo
			temp(1:6,1)=temp(1:6,1)*cos((temp(1:6,2)+WavePhase)/180.0*PI)
			write(132,'(i4,6(e13.4,f7.2))')j,((temp(k,kk),kk=1,2),k=1,6)
		enddo
! 	enddo
	!传递函数输出
	elseif(TRIM(ADJUSTL(ansMethod))=="NO")	then 
! 		do i=1,NP
			write(132,'(a,f5.1,a,f8.3)')"# HEAD= ",Beta*180.0/PI," Ome0= ",omei0
			write(132,'(a)')"#  ome0  ome"
			write(132,'(2f7.3)')omei0,omeI1
			write(132,'(a,6x,a,2x,3x,a,1x,2(3x,a,2x,3x,a,1x),5x,a,2x,3x,a,1x,2(3x,a,2x,3x,a,1x))')"#",&
						    "  F.X.","PHA"," F.Y. ","PHA"," F.Z. ","PHA","M.X.","PHA","  M.Y.","PHA","  M.Z.","PHA"
			do j=1,Nbsect
				do k=1,6
					call comp_to_ampha(SectLoad(1,k,j),SectLoad(2,k,j),temp(k,1),temp(k,2))
				enddo
				temp(1:6,1)=temp(1:6,1)/(ampI)
				write(132,'(i4,6(e11.4,f7.2))')j,((temp(k,kk),kk=1,2),k=1,6)
			enddo
! 		enddo
	endif

end subroutine dsls_putouts
!-----------------------------------------------------------------------------------
!------------------------------------------------
END MODULE SECTLOAD_MOD