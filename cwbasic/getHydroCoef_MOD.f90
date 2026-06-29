MODULE getHydroCoef_MOD

USE INOUTACCESS_MOD
USE ENVIRONMENT_MOD,ONLY:ZG,ROU,G0,H,PI,NBSECT ,SECTPSN,SECTTYPE
USE MESH_MOD,ONLY:NB,VOL,SL,Factor,TarNum,Space
USE CAL_MOD,ONLY:WATERDEPTH


IMPLICIT NONE
PRIVATE
PUBLIC::getHydroCoef,getWaveExistForce,getHydroCoef_sect,getWaveExistForce_sect


CONTAINS
!**************************************************************************
!程序功能：计算三维水动力系数及单位波幅的绕射力         
!输出参数：F(1:2,1:6,1:7,1:NBOME)为三维水动力系数及单位波幅的绕射力         
!**************************************************************************
subroutine getHydroCoef(U,ea0,e0,xav0,omeI,phi,F) !UKN(IU),ea(1:NB),e(:,:,1:NB),xav(:,1:NB),ome(IP),phi(:,:,:,:,1:NB),F(:,:,1:7) 
implicit none

real(8),intent(in)::U
real(8),intent(in),dimension(1:NB)::ea0
real(8),intent(in),dimension(1:3,1:3,1:NB)::e0 !面元局部坐标系下法相坐标（1~3xyz,？？，1~湿表面）
real(8),intent(in),dimension(1:3,1:NB)::xav0 !中心点坐标
! real(8),intent(in),dimension(1:NBOME)::omeI
real(8),intent(in)::omeI
real(8),intent(in),dimension(1:2,1:2,1:7,1:Factor,1:NB)::phi !(Phi(实部、虚部，？？，辐射1~6绕射7，对称性因子，湿表面序号）
real(8),intent(out),dimension(1:2,1:6,1:7)::F !六自由度，7个速度势
real(8),allocatable,dimension(:,:,:)::n 
integer(4)::i,j,k,L,count
allocate(n(1:6,1:Factor,1:NB)) 

do i=1,NB
	do j=1,3
		n(J,1,I)=e0(j,3,i)
	enddo
	n(4,1,I)=xav0(2,i)*e0(3,3,i)-(xav0(3,i)-ZG)*e0(2,3,i) !将原点移动到重心处
	n(5,1,I)=(xav0(3,i)-ZG)*e0(1,3,i)-xav0(1,i)*e0(3,3,i)
	n(6,1,I)=xav0(1,i)*e0(2,3,i)-xav0(2,i)*e0(1,3,i)
enddo
L=0
do count=1,TarNum,Space !不对称，TarNum=1,Space=1
      L=L+1
select   case(count)

               case(1)
                       n(1,L,1:NB)=n(1,1,1:NB)
                       n(2,L,1:NB)=n(2,1,1:NB)
                       n(3,L,1:NB)=n(3,1,1:NB)
                       n(4,L,1:NB)=n(4,1,1:NB)
                       n(5,L,1:NB)=n(5,1,1:NB)
                       n(6,L,1:NB)=n(6,1,1:NB)             
                case(2)
                       n(1,L,1:NB)=-n(1,1,1:NB)
                       n(2,L,1:NB)=n(2,1,1:NB)
                       n(3,L,1:NB)=n(3,1,1:NB)
                       n(4,L,1:NB)=n(4,1,1:NB)
                       n(5,L,1:NB)=-n(5,1,1:NB)
                       n(6,L,1:NB)=-n(6,1,1:NB)
               case(3)
                        n(1,L,1:NB)=-n(1,1,1:NB)
                        n(2,L,1:NB)=-n(2,1,1:NB)
                        n(3,L,1:NB)=n(3,1,1:NB)
                        n(4,L,1:NB)=-n(4,1,1:NB)
                        n(5,L,1:NB)=-n(5,1,1:NB)
                        n(6,L,1:NB)=n(6,1,1:NB)
                case(4)
                        n(1,L,1:NB)=n(1,1,1:NB)
                        n(2,L,1:NB)=-n(2,1,1:NB)
                        n(3,L,1:NB)=n(3,1,1:NB)
                        n(4,L,1:NB)=-n(4,1,1:NB)
                        n(5,L,1:NB)=n(5,1,1:NB)
                        n(6,L,1:NB)=-n(6,1,1:NB)
 end select
end do

!!!!!!!!!!!!!!!!!!!!!!!!!          1~4            辐射    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
F=0.0 !F(1:2,1:6,1:7) 实部虚部，6个自由度运动模态，，辐射绕射  !(Phi(实部、虚部，？1：不求导、2：负的对x求导？？，辐射1~6绕射7，对称性因子，湿表面序号）
! do IP=1,NBOME
	do I=1,6  !6个运动方向
        do j=1,4	!1~4辐射
		      do L=1,Factor !1~1对称因子
				   do k=1,NB
					F(1,i,j)=F(1,i,j)+ea0(k)*phi(1,1,j,L,k)*n(I,L,K)&
								-ea0(k)*phi(2,2,j,L,k)*n(I,L,K)*U/omeI ! 附加质量系数
					F(2,i,j)=F(2,i,j)-ea0(k)*phi(2,1,j,L,k)*n(I,L,K)&
								-ea0(k)*phi(1,2,j,L,k)*n(I,L,K)*U/omeI
				 enddo
			end do
				F(2,i,j)=F(2,i,j)*omeI !阻尼系数
		enddo
 !!!!!!!!!!!!!!!!!!!!!!!!!!       5 辐射          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
      do L=1,Factor
			do k=1,NB
				F(1,i,5)=F(1,i,5)+ea0(k)*phi(1,1,5,L,k)*n(I,L,K)&
							-ea0(k)*phi(2,2,5,L,k)*n(I,L,K)*U/omeI&
							+ea0(k)*phi(2,1,3,L,k)*n(I,L,K)*U/omeI&
							+ea0(k)*phi(1,2,3,L,k)*n(I,L,K)*(U/omeI)**2
				F(2,i,5)=F(2,i,5)-ea0(k)*phi(2,1,5,L,k)*n(I,L,K)&
							-ea0(k)*phi(1,2,5,L,k)*n(I,L,K)*U/omeI&
							+ea0(k)*phi(1,1,3,L,k)*n(I,L,K)*U/omeI&
							-ea0(k)*phi(2,2,3,L,k)*n(I,L,K)*(U/omeI)**2
			enddo
	 end do
			F(2,i,5)=F(2,i,5)*omeI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     6 辐射      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do L=1,Factor
			do k=1,NB
				F(1,i,6)=F(1,i,6)+ea0(k)*phi(1,1,6,L,k)*n(I,L,K)&
							-ea0(k)*phi(2,2,6,L,k)*n(I,L,K)*U/omeI&
							-ea0(k)*phi(2,1,2,L,k)*n(I,L,K)*U/omeI&
							-ea0(k)*phi(1,2,2,L,k)*n(I,L,K)*(U/omeI)**2
				F(2,i,6)=F(2,i,6)-ea0(k)*phi(2,1,6,L,k)*n(I,L,K)&
							-ea0(k)*phi(1,2,6,L,k)*n(I,L,K)*U/omeI&
							-ea0(k)*phi(1,1,2,L,k)*n(I,L,K)*U/omeI&
							+ea0(k)*phi(2,2,2,L,k)*n(I,L,K)*(U/omeI)**2
			enddo
	  end do
			F(2,i,6)=F(2,i,6)*omeI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      7辐射      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do L=1,Factor
			do k=1,NB
				F(1,i,7)=F(1,i,7)+ea0(k)*phi(2,1,7,L,k)*n(I,L,K)*omeI&
							+ea0(k)*phi(1,2,7,L,k)*n(I,L,K)*U !单位波幅绕射力 实部
				F(2,i,7)=F(2,i,7)-ea0(k)*phi(1,1,7,L,k)*n(I,L,K)*omeI&
							+ea0(k)*phi(2,2,7,L,k)*n(I,L,K)*U !单位波幅绕射力 虚部
			enddo
	end do
enddo

F=F*rou

deallocate(n)

end subroutine getHydroCoef
!---------------------------------------------------------------------


!**************************************************************************
!程序功能：计算波浪主干扰力和波浪干扰力         
!输出参数：F(1:NP,1:6,8:9,1:2)为波浪主干扰力(8)和波浪干扰力(9)
!改动：  对频率循环放在外面，以解决VF6.6堆栈溢出问题。 2002.3.29
!***************************************************************************
subroutine getWaveExistForce(OMEI,OMEI0,AMPI,k0,beta,ea0,e0,xav0,F,IP)
! call getWaveExistForce (k0,HEAD(IB),ea(1:NB),e(:,:,1:NB),xav(:,1:NB),F(:,:,7:9,:)) !计算波浪干扰力
implicit none

real(8),intent(in)::beta,K0,AMPI,OMEI,OMEI0
real(8),intent(in),dimension(1:NB)::ea0
real(8),intent(in),dimension(1:3,1:3,1:NB)::e0
real(8),intent(in),dimension(1:3,1:NB)::xav0
! real(8),intent(in),dimension(1:2,1:4,1:NB)::xq0
real(8),intent(inout),dimension(1:2,1:6,1:9)::F
! real(8),intent(out),dimension(1:2,1:6)::Fks
real(8),allocatable,dimension(:,:)::x00(:,:)
real(8),allocatable,dimension(:,:,:)::n
real(8)::temp
integer(4)::IP
integer(4)::i,j,k,L,count
allocate(n(1:NB,1:Factor,1:6),X00(1:3,1:NB))
do i=1,NB
	do j=1,3
		n(i,1,j)=e0(j,3,i)
	enddo
	n(i,1,4)=xav0(2,i)*e0(3,3,i)-(xav0(3,i)-ZG)*e0(2,3,i)
	n(i,1,5)=(xav0(3,i)-ZG)*e0(1,3,i)-xav0(1,i)*e0(3,3,i)
	n(i,1,6)=xav0(1,i)*e0(2,3,i)-xav0(2,i)*e0(1,3,i)
enddo

!计算波浪主干扰力

F(1:2,1:6,8:9)=0.0
L=0
do count=1,TarNum,Space
      L=L+1
select   case(count)

               case(1)
                       x00=xav0
                        n(1:NB,L,1)=n(1:NB,1,1)
                        n(1:NB,L,2)=n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=n(1:NB,1,4)
                        n(1:NB,L,5)=n(1:NB,1,5)
                        n(1:NB,L,6)=n(1:NB,1,6)        
                case(2)

                                x00(1,:)=-xav0(1,:)
                                x00(2,:)=xav0(2,:)
                                x00(3,:)=xav0(3,:)
                        n(1:NB,L,1)=-n(1:NB,1,1)
                        n(1:NB,L,2)=n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=n(1:NB,1,4)
                        n(1:NB,L,5)=-n(1:NB,1,5)
                        n(1:NB,L,6)=-n(1:NB,1,6)
               case(3)
                                x00(1,:)=-xav0(1,:)
                                x00(2,:)=-xav0(2,:)
                                x00(3,:)=xav0(3,:)
                        n(1:NB,L,1)=-n(1:NB,1,1)
                        n(1:NB,L,2)=-n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=-n(1:NB,1,4)
                        n(1:NB,L,5)=-n(1:NB,1,5)
                        n(1:NB,L,6)=n(1:NB,1,6)
                case(4)
                                x00(1,:)=xav0(1,:)
                                x00(2,:)=-xav0(2,:)
                                x00(3,:)=xav0(3,:)
                        n(1:NB,L,1)=n(1:NB,1,1)
                        n(1:NB,L,2)=-n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=-n(1:NB,1,4)
                        n(1:NB,L,5)=n(1:NB,1,5)
                        n(1:NB,L,6)=-n(1:NB,1,6)
 end select




	if((trim(adjustl(WATERDEPTH))=="FINITE"))then  !有限水深
		temp=rou*g0*ampI/cosh(k0*H)
	elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		temp=rou*g0*ampI
	endif

	do j=1,6
		do k=1,NB
			if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
			!1部分
				F(1,j,8)=F(1,j,8)+temp*cosh(k0*(x00(3,k)+H))*cos(k0*(x00(1,k)*cos(beta)-&
										x00(2,k)*sin(beta)))*n(k,L,j)*ea0(k)	
			!1部分					
				F(2,j,8)=F(2,j,8)+temp*cosh(k0*(x00(3,k)+H))*sin(k0*(x00(1,k)*cos(beta)-&
										(x00(2,k))*sin(beta)))*n(k,L,j)*ea0(k)	

										
			elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
			!1部分	
				F(1,j,8)=F(1,j,8)+temp*exp(k0*x00(3,k))*cos(k0*(x00(1,k)*cos(beta)-&
										x00(2,k)*sin(beta)))*n(k,L,j)*ea0(k)	

			!1部分				
				F(2,j,8)=F(2,j,8)+temp*exp(k0*x00(3,k))*sin(k0*(x00(1,k)*cos(beta)-&
										(x00(2,k))*sin(beta)))*n(k,L,j)*ea0(k)	

			endif
		enddo
	enddo
end do
!计算波浪干扰力
F(:,:,7)=ampI*F(:,:,7)
F(:,:,9)=F(:,:,7)+F(:,:,8)

IF(IP>0)THEN
    call f_putouts(F,BETA,OMEI,OMEI0)
END IF


deallocate(n,x00)
end subroutine
!----------------------------------------------------------------------

!**************************************************************************
!程序功能：计算三维水动力系数及单位波幅的绕射力 (含局部分段对应项)        
!输出参数：F(1:2,1:6,1:7,1:NBOME)为三维水动力系数及单位波幅的绕射力         
!**************************************************************************
subroutine getHydroCoef_sect(U,ea0,e0,xav0,omeI,phi,F,sectf)
implicit none

real(8),intent(in)::U
real(8),intent(in),dimension(1:NB)::ea0
real(8),intent(in),dimension(1:3,1:3,1:NB)::e0
real(8),intent(in),dimension(1:3,1:NB)::xav0
real(8),intent(in)::omeI
real(8),intent(in),dimension(1:2,1:2,1:7,1:Factor,1:NB)::phi
real(8),intent(out),dimension(1:2,1:6,1:7)::F
real(8),intent(out),dimension(1:2,1:6,1:nbsect,1:7)::sectF
real(8),allocatable,dimension(:,:,:)::n
real(8),allocatable,dimension(:,:)::x00
real(8)::tempc(1:7),temps(1:7)
integer(4)::i,j,k,l,ii,jj,LK,count
allocate(n(1:6,1:Factor,1:NB),X00(1:3,1:NB))
F=0.0
sectf=0.0
do i=1,NB
	do j=1,3
		n(J,1,I)=e0(j,3,i)
	enddo
	n(4,1,I)=xav0(2,i)*e0(3,3,i)-(xav0(3,i)-ZG)*e0(2,3,i)
	n(5,1,I)=(xav0(3,i)-ZG)*e0(1,3,i)-xav0(1,i)*e0(3,3,i)
	n(6,1,I)=xav0(1,i)*e0(2,3,i)-xav0(2,i)*e0(1,3,i)
enddo
L=0
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

!计算波浪主干扰力

do I=1,6
    do k=1,NB
	    do J=1,4
	        tempc(j)=ea0(k)*phi(1,1,j,L,k)*n(I,L,K)-ea0(k)*phi(2,2,j,L,k)*n(I,L,K)*U/omeI
			temps(j)=(-ea0(k)*phi(2,1,j,L,k)*n(I,L,K)-ea0(k)*phi(1,2,j,L,k)*n(I,L,K)*U/omeI)*omeI
		    F(1,i,j)=F(1,i,j)+tempc(j)
		    F(2,i,j)=F(2,i,j)+temps(j)
		    
	    enddo

	    
	    tempc(5)=ea0(k)*phi(1,1,5,L,k)*n(I,L,K)&
				    -ea0(k)*phi(2,2,5,L,k)*n(I,L,K)*U/omeI&
				    +ea0(k)*phi(2,1,3,L,k)*n(I,L,K)*U/omeI&
				    +ea0(k)*phi(1,2,3,L,k)*n(I,L,K)*(U/omeI)**2
		temps(5)=(-ea0(k)*phi(2,1,5,L,k)*n(I,L,K)&
				    -ea0(k)*phi(1,2,5,L,k)*n(I,L,K)*U/omeI&
				    +ea0(k)*phi(1,1,3,L,k)*n(I,L,K)*U/omeI&
				    -ea0(k)*phi(2,2,3,L,k)*n(I,L,K)*(U/omeI)**2)*omeI


	    F(1,i,5)=F(1,i,5)+tempc(5)
	    F(2,i,5)=F(2,i,5)+temps(5)
	    

        tempc(6)=ea0(k)*phi(1,1,6,L,k)*n(I,L,K)&
				    -ea0(k)*phi(2,2,6,L,k)*n(I,L,K)*U/omeI&
				    -ea0(k)*phi(2,1,2,L,k)*n(I,L,K)*U/omeI&
				    -ea0(k)*phi(1,2,2,L,k)*n(I,L,K)*(U/omeI)**2
        temps(6)=(-ea0(k)*phi(2,1,6,L,k)*n(I,L,K)&
				    -ea0(k)*phi(1,2,6,L,k)*n(I,L,K)*U/omeI&
				    -ea0(k)*phi(1,1,2,L,k)*n(I,L,K)*U/omeI&
				    +ea0(k)*phi(2,2,2,L,k)*n(I,L,K)*(U/omeI)**2)*omeI
	    F(1,i,6)=F(1,i,6)+tempc(6)
	    F(2,i,6)=F(2,i,6)+temps(6)
	    


        tempc(7)=ea0(k)*phi(2,1,7,L,k)*n(I,L,K)*omeI&
				    +ea0(k)*phi(1,2,7,L,k)*n(I,L,K)*U
        temps(7)=-ea0(k)*phi(1,1,7,L,k)*n(I,L,K)*omeI&
				    +ea0(k)*phi(2,2,7,L,k)*n(I,L,K)*U


	    F(1,i,7)=F(1,i,7)+tempc(7)
	    F(2,i,7)=F(2,i,7)+temps(7)

    DO ii=1,NBSECT
            IF(SECTTYPE(II)==1)THEN
                LK=1
            ELSEIF(SECTTYPE(II)==2)THEN
                LK=2
            ENDIF
            IF(x00(LK,K)<=SECTPSN(LK,ii))THEN
                sectf(1,i,ii,1:7)=sectf(1,i,ii,1:7)+tempc(1:7)
                sectf(2,i,ii,1:7)=sectf(2,i,ii,1:7)+temps(1:7)
            ENDIF
    ENDDO 
     enddo
  enddo
end do
F=F*rou
sectf=sectf*rou

deallocate(n,x00)

end subroutine getHydroCoef_sect
!---------------------------------------------------------------------



!**************************************************************************
!程序功能：计算波浪主干扰力和波浪干扰力 (含局部)        
!输出参数：F(1:NP,1:6,8:9,1:2)为波浪主干扰力(8[入射波力])和波浪干扰力(9[入射波力与绕射波力之和])
!***************************************************************************
subroutine getWaveExistForce_sect(OMEI,OMEI0,AMPI,k0,beta,ea0,e0,xav0,F,sectf,IP)
implicit none

real(8),intent(in)::beta,K0,AMPI,OMEI,OMEI0
real(8),intent(in),dimension(1:NB)::ea0
real(8),intent(in),dimension(1:3,1:3,1:NB)::e0
real(8),intent(in),dimension(1:3,1:NB)::xav0
real(8),intent(inout),dimension(1:2,1:6,1:9)::F
real(8),intent(inout),dimension(1:2,1:6,1:nbsect,1:9)::sectF
real(8),allocatable,dimension(:,:,:)::n
real(8),allocatable,dimension(:,:)::x00
real(8)::temp
real(8)::tempc(8:9),temps(8:9)
integer(4)::IP
integer(4)::i,j,k,l,ii,LK,count
allocate(n(1:NB,1:4,1:6),x00(1:3,1:NB))

do i=1,NB
	do j=1,3
		n(i,1,j)=e0(j,3,i)
	enddo
	n(i,1,4)=xav0(2,i)*e0(3,3,i)-(xav0(3,i)-ZG)*e0(2,3,i)
	n(i,1,5)=(xav0(3,i)-ZG)*e0(1,3,i)-xav0(1,i)*e0(3,3,i)
	n(i,1,6)=xav0(1,i)*e0(2,3,i)-xav0(2,i)*e0(1,3,i)
enddo
L=0
!计算波浪主干扰力
F(1:2,1:6,8)=0.0
sectf(1:2,1:6,:,8)=0.0


if((trim(adjustl(WATERDEPTH))=="FINITE"))then  !有限水深
		temp=rou*g0*ampI/cosh(k0*H)
elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		temp=rou*g0*ampI
endif
do count=1,TarNum,Space
      L=L+1
select   case(count)

               case(1)
                       x00=xav0
                        n(1:NB,L,1)=n(1:NB,1,1)
                        n(1:NB,L,2)=n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=n(1:NB,1,4)
                        n(1:NB,L,5)=n(1:NB,1,5)
                        n(1:NB,L,6)=n(1:NB,1,6)        
                case(2)

                                x00(1,:)=-xav0(1,:)
                                x00(2,:)=xav0(2,:)
                                x00(3,:)=xav0(3,:)
                        n(1:NB,L,1)=-n(1:NB,1,1)
                        n(1:NB,L,2)=n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=n(1:NB,1,4)
                        n(1:NB,L,5)=-n(1:NB,1,5)
                        n(1:NB,L,6)=-n(1:NB,1,6)
               case(3)
                                x00(1,:)=-xav0(1,:)
                                x00(2,:)=-xav0(2,:)
                                x00(3,:)=xav0(3,:)
                        n(1:NB,L,1)=-n(1:NB,1,1)
                        n(1:NB,L,2)=-n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=-n(1:NB,1,4)
                        n(1:NB,L,5)=-n(1:NB,1,5)
                        n(1:NB,L,6)=n(1:NB,1,6)
                case(4)
                                x00(1,:)=xav0(1,:)
                                x00(2,:)=-xav0(2,:)
                                x00(3,:)=xav0(3,:)
                        n(1:NB,L,1)=n(1:NB,1,1)
                        n(1:NB,L,2)=-n(1:NB,1,2)
                        n(1:NB,L,3)=n(1:NB,1,3)
                        n(1:NB,L,4)=-n(1:NB,1,4)
                        n(1:NB,L,5)=n(1:NB,1,5)
                        n(1:NB,L,6)=-n(1:NB,1,6)
 end select
	do j=1,6
		do k=1,NB

			if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
			    tempc(8)=temp*cosh(k0*(xav0(3,k)+H))*cos(k0*(x00(1,K)*cos(beta)-&
										x00(2,K)*sin(beta)))*n(k,L,j)*ea0(k)	
			    temps(8)=temp*cosh(k0*(xav0(3,k)+H))*sin(k0*(x00(1,K)*cos(beta)-&
										x00(2,K)*sin(beta)))*n(k,L,j)*ea0(k)	
				F(1,j,8)=F(1,j,8)+tempc(8)
				F(2,j,8)=F(2,j,8)+temps(8)
			elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
			
				tempc(8)=temp*exp(k0*xav0(3,k))*cos(k0*(x00(1,K)*cos(beta)-&
										x00(2,K)*sin(beta)))*n(k,L,j)*ea0(k)	
				temps(8)=temp*exp(k0*xav0(3,k))*sin(k0*(x00(1,K)*cos(beta)-&
										x00(2,K)*sin(beta)))*n(k,L,j)*ea0(k)	
				F(1,j,8)=F(1,j,8)+tempc(8)
				F(2,j,8)=F(2,j,8)+temps(8)
			endif
            DO ii=1,NBSECT
                IF(SECTTYPE(II)==1)THEN
                    LK=1
                ELSEIF(SECTTYPE(II)==2)THEN
                    LK=2
                ENDIF
                IF(x00(LK,K)<=SECTPSN(LK,ii))THEN
                    sectf(1,J,ii,8)=sectf(1,J,ii,8)+tempc(8)
                    sectf(2,J,ii,8)=sectf(2,J,ii,8)+temps(8)
                ENDIF
            ENDDO 			
		enddo
	enddo
end do
!计算波浪干扰力
F(:,:,7)=ampI*F(:,:,7)
sectf(:,:,:,7)=ampI*sectf(:,:,:,7)
F(:,:,9)=F(:,:,7)+F(:,:,8)
sectf(:,:,:,9)=sectf(:,:,:,7)+sectf(:,:,:,8)  
if(IP>0)   then           
    call f_putouts(F,BETA,OMEI,OMEI0)
endif
deallocate(n,x00)
end subroutine getWaveExistForce_sect
!----------------------------------------------------------------------

!-----------------------------------------------------------
!输出水动力系数结果
!-----------------------------------------------------------
subroutine f_putouts(F,BETA,OMEI,OMEI0)
implicit none
! integer::flag
real(8),intent(in)::beta,OMEI,OMEI0
real(8),intent(in),dimension(1:2,1:6,1:9)::F
real(8),dimension(1:6,1:6)::NDF	   !无因次水动力系数
real(8),dimension(1:2,1:6)::temp
real(8)::NDome
! real(8)::PI
integer::i,j,k
! CHARACTER(LEN=80)::OUTPUT3HC

! allocate(NDome(1:NBOME))
! PI=4.0*atan(1.0)
! 		do k=1,NBOME-----------------------------------
		    write(5,'(a,F5.1,a,f8.3)')"HEAD= ",BETA*180.0/PI,"  Ome0= ",omeI0
			write(5,*)
			write(5,*)"附加质量A, (波浪绕射力Fd, 波浪主干扰力Fw, 波浪干扰力F)实部:"
			write(5,'(9(4x,a))')"a(1:6,1)","a(1:6,2)","a(1:6,3)","a(1:6,4)","a(1:6,5)","a(1:6,6)","   Fd   ","   Fw   ","    F   "
			do i=1,6
				write(5,'(9e12.4)')(F(1,i,j),j=1,9)
			enddo

			write(5,*)"阻尼系数B, (波浪绕射力Fd, 波浪主干扰力Fw, 波浪干扰力F)虚部:"
			write(5,'(9(4x,a))')"b(1:6,1)","b(1:6,2)","b(1:6,3)","b(1:6,4)","b(1:6,5)","b(1:6,6)","   Fd   ","   Fw   ","    F   "
			do i=1,6
				write(5,'(9e12.4)')(F(2,i,j),j=1,9)
			enddo
! 		enddo-------------------------------------------
!---------------------------------------------------------------------------------------
		write(5,*)
		write(5,*)"绘图数据部分："
! 		select case(flag)
! 		case(1)
! 			write(5,*)"附加质量(ome*sqrt(l/g),a22,a24,a33,a55,a44,a66):" !方箱
! 		case(2)
! 			write(5,*)"附加质量(ome*sqrt(l/g),a22,a33):"  !半球
! 		case(3,4)
! 			write(5,*)"附加质量(ome*sqrt(l/g),a22,a24,a26,a33,a35,a42,a44,a46,a53,a55,a62,a64,a66):" !实船
! 		case(5)
! 			write(5,*)"附加质量(ome*sqrt(l/g),a22,a33,a44):" !圆盘
! 		case default
			write(5,'(a)')"附加质量(ome0,a11,a13,a15,a22,a24,a26,a31,a33,a35,a42,a44,a46,a51,a53,a55,a62,a64,a66):" !实船
			write(5,'(3x,a,2x,18(5x,a,4x))')"ome0","a11","a13","a15","a22","a24","a26","a31","a33","a35","a42","a44","a46","a51","a53","a55","a62","a64","a66" !实船
! 		end select
! 		do i=1,NBOME----------------------------------------------
! 			select case(flag)
! 			case(1)
! 				call NonDimension(ome(i),F(:,:,1:6,I),NDome(i),NDF(:,:,:))	
! 				write(5,'(f8.3,6e12.4)')NDome(i),NDF(1,2,2),NDF(1,2,4),NDF(1,3,3),NDF(1,5,5),NDF(1,4,4),NDF(1,6,6)
! 			case(2)
! 				call NonDimension(ome(i),F(:,:,1:6,I),NDome(i),NDF(:,:,:))	
! 				write(5,'(f5.3,2e10.3)')NDome(i),NDF(1,2,2),NDF(1,3,3) !半球
! 			case(3,4)
! 				call NonDimension(ome(i),F(:,:,1:6,I),NDome(i),NDF(:,:,:))
! 				write(5,'(f8.3,13e12.4)')NDome(i),NDF(1,2,2),NDF(1,2,4),NDF(1,2,6),NDF(1,3,3),&
! 						  NDF(1,3,5),NDF(1,4,2),NDF(1,4,4),NDF(1,4,6),NDF(1,5,3),NDF(1,5,5),&
! 						  NDF(1,6,2),NDF(1,6,4),NDF(1,6,6)
! 			case(5)
! 				call NonDimension(ome(i),F(:,:,1:6,I),NDome(i),NDF(:,:,:))
! 				write(5,'(f5.3,3e10.3)')NDome(i),NDF(1,2,2),NDF(1,3,3),NDF(1,4,4) !圆盘	
! 			case default
				call NonDimension(omeI,F(1,:,1:6),NDome,NDF(:,:))
				write(5,'(f8.3,18e12.4)')omeI0,NDF(1,1),NDF(1,3),NDF(1,5),NDF(2,2),NDF(2,4),NDF(2,6),NDF(3,1),NDF(3,3),&
						  NDF(3,5),NDF(4,2),NDF(4,4),NDF(4,6),NDF(5,1),NDF(5,3),NDF(5,5),NDF(6,2),NDF(6,4),NDF(6,6)
! 			end select
! 		enddo----------------------------------------------------

! 		select case(flag)
! 		case(1)
! 			write(5,*)"阻尼系数(ome*sqrt(l/g),b22,b24,b33,b55,b44,b66):" !方箱
! 		case(2)
! 			write(5,*)"阻尼系数(ome*sqrt(l/g),b22,b33):"  !半球
! 		case(3,4)
! 			write(5,*)"阻尼系数(ome*sqrt(l/g),b22,b24,b26,b33,b35,b42,b44,b46,b53,b55,b62,b64,b66):" !实船
! 		case(5)
! 			write(5,*)"阻尼系数(ome*sqrt(l/g),b22,b33,b44):" !圆盘
! 		case default
			write(5,'(a)')"阻尼系数(om0,b11,b13,b15,b22,b24,b26,b31,b33,b35,b42,b44,b46,b51,b53,b55,b62,b64,b66):" !实船
			write(5,'(3x,a,2x,18(5x,a,4x))')"ome0","b11","b13","b15","b22","b24","b26","b31","b33","b35","b42","b44","b46","b51","b53","b55","b62","b64","b66" !实船
! 		end select
! 		do i=1,NBOME
! 			select case(flag)
! 			case(1)
! 				write(5,'(f8.3,6e12.4)')NDome(i),NDF(2,2,2),NDF(2,2,4),NDF(2,3,3),NDF(2,5,5),NDF(2,4,4),NDF(2,6,6)
! 			case(2)
! 				write(5,'(f5.3,2e10.3)')NDome(i),NDF(2,2,2),NDF(2,3,3) !半球
! 			case(3,4)
! 				write(5,'(f8.3,13e12.4)')NDome(i),NDF(2,2,2),NDF(2,2,4),NDF(2,2,6),NDF(2,3,3),&
! 						  NDF(2,3,5),NDF(2,4,2),NDF(2,4,4),NDF(2,4,6),NDF(2,5,3),NDF(2,5,5),&
! 						  NDF(2,6,2),NDF(2,6,4),NDF(2,6,6)
! 			case(5)
! 				write(5,'(f5.3,3e10.3)')NDome(i),NDF(2,2,2),NDF(2,3,3),NDF(2,4,4) !圆盘
! 			case default
                call NonDimension(omeI,F(2,:,1:6),NDome,NDF(:,:))
				write(5,'(f8.3,18e12.4)')omeI0,NDF(1,1),NDF(1,3),NDF(1,5),NDF(2,2),NDF(2,4),NDF(2,6),NDF(3,1),NDF(3,3),&
						  NDF(3,5),NDF(4,2),NDF(4,4),NDF(4,6),NDF(5,1),NDF(5,3),NDF(5,5),NDF(6,2),NDF(6,4),NDF(6,6)
				write(5,*)
! 			end select
! 		enddo

write(5,'(a)')"波浪干扰力：" !实船
write(5,'(3x,a,2x,6(5x,a,4x,2x,a,2x))')"ome0","Fa1","PHA","Fa2","PHA","Fa3","PHA","Fa4","PHA","Fa5","PHA","Fa6","PHA"
! do i=1,NBOME
	do j=1,6
		call comp_to_ampha(F(1,j,9),F(2,j,9),temp(1,J),temp(2,J))
	enddo
	write(5,'(1x,f7.3,6(e12.4,f7.2))')omeI0,((temp(K,J),k=1,2),j=1,6)
	write(5,*)
! enddo
! deallocate(NDome)
end subroutine f_putouts
!-------------------------------------------------------------------



!**************************************************************************
!程序功能：对三维水动力系数无因次化处理
!输入参数：ome为波浪遭遇频率,Vol为排水体积，sl为船长，F为三维水动力系数         
!输出参数：NDome为无因次化波浪遭遇频率，NDF为无因次化三维水动力系数
!***************************************************************************
subroutine NonDimension(omeI,F,NDome,NDF)
implicit none
real(8),intent(in)::omeI
real(8),intent(in),dimension(1:6,1:6)::F
real(8),intent(out)::NDome
real(8),intent(out),dimension(1:6,1:6)::NDF
integer(4)::i,j,mi,mj
NDome=omeI*sqrt(sl/(g0))
do i=1,6
	do j=1,6
		mi=0;mj=0
		if(i>3)mi=1
		if(j>3)mj=1
		NDF(i,j)=F(i,j)  !/(vol*sl**(mi+mj))
! 		NDF(2,i,j)=F(2,i,j)  !/(vol*sl**(mi+mj))*sqrt(sl/(g0))
	enddo
enddo
end subroutine NonDimension
!----------------------------------------------------------------------



!--------------------------------------------------------------
END MODULE getHydroCoef_MOD
