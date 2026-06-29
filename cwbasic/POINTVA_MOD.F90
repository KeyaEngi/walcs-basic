MODULE POINTVA_MOD

USE INOUTACCESS_MOD
USE PVPOINT_MOD,ONLY:NPVA,VAPOINT
USE ENVIRONMENT_MOD,ONLY:ZG,PI,XG,g0
USE MESH_MOD,ONLY:DELTAX,ST,SL,TY,Ta
USE PRINT_MOD,ONLY:ANSMETHOD


IMPLICIT NONE
PRIVATE
PUBLIC::GETVAPOINT,shippointva

CONTAINS
!---------------------------------------------------------------
subroutine GETVAPOINT()
IMPLICIT NONE

CHARACTER(LEN=300)::PVAINPUT
real(8),dimension(1:3)::tempcoor
REAL::ratio
INTEGER::TEMP,I,J

PVAINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.PVA'
open(28,file=PVAINPUT,status='old')
call IUTMP(28)
read(28,*)(tempcoor(i),i=1,3),ratio		!随船平动坐标系原点在压力计算点所在坐标系下坐标
tempcoor(:)=tempcoor(:)/ratio
read(28,*)NPva							!计算点数(分别为总数，X,Y,Z向加速度计算点数)
allocate(vaPoint(1:3,1:NPva))

!do i=1,NPva
!	read(28,*)TEMP,(vaPoint(j,i),j=1,3)	
!	vaPoint(1,i)=vaPoint(1,i)/ratio+tempcoor(1)-XG
!	vaPoint(2,i)=vaPoint(2,i)/ratio+tempcoor(2)
!	vaPoint(3,i)=vaPoint(3,i)/ratio+tempcoor(3)-ST
!enddo
do i=1,NPva
	read(28,*)TEMP,(vaPoint(j,i),j=1,3)	
	vaPoint(1,i)=vaPoint(1,i)/ratio+tempcoor(1)-XG
	vaPoint(2,i)=vaPoint(2,i)/ratio+tempcoor(2)
	vaPoint(3,i)=vaPoint(3,i)/ratio+tempcoor(3)-ST
enddo
vaPoint(3,:)=vaPoint(3,:)-ZG
!write(*,*)vaPoint(1,1),vaPoint(2,1),vaPoint(3,1)
!----------纵倾坐标变换----------
do i=1,NPva
	vaPoint(1:3,i)=Matmul(Ty,vaPoint(1:3,i))
enddo
close(28)

END subroutine GETVAPOINT
!-------------------------------------------------------------------------


!************************************************************************************
!程序功能：求解计算点处船舶表面速度和加速度
!输入参数：NPva计算点数目;vaPoint(1:3,1:NPVA)为计算点坐标
!输出参数：pva(1:Npva,1:6,1:2)分别为Vx,Vy,Vz,Ax,Ay,Az的实部和虚部
!程序编制：张海彬    时间：2002年5月23日      
!程序修改：孙葳       时间 ：2012年7月16日
!************************************************************************************
subroutine shippointva(BETA,OMEI0,OMEI,OMEI1,k0,WAVEPHASE,motion,IP)
implicit none

real(8),INTENT(IN)::OMEI0,OMEI,OMEI1,k0,WAVEPHASE,BETA
real(8),intent(in),dimension(1:2,1:6)::motion
real(8),ALLOCATABLE,dimension(:,:,:)::pva
integer(4)::IP
integer::i 
!----------------------------------------------------------------------
allocate(pva(1:2,1:6,1:NPva))

do i=1,NPva
	!速度
	pva(1,1,I)=-omeI*(motion(2,1)+vaPoint(3,i)*motion(2,5)-vaPoint(2,i)*motion(2,6))
	pva(2,1,I)=omeI*(motion(1,1)+vaPoint(3,i)*motion(1,5)-vaPoint(2,i)*motion(1,6))
	pva(1,2,I)=-omeI*(motion(2,2)-vaPoint(3,i)*motion(2,4)+vaPoint(1,i)*motion(2,6))
	pva(2,2,I)=omeI*(motion(1,2)-vaPoint(3,i)*motion(1,4)+vaPoint(1,i)*motion(1,6))
	pva(1,3,I)=-omeI*(motion(2,3)-vaPoint(1,i)*motion(2,5)+vaPoint(2,i)*motion(2,4))
	pva(2,3,I)=omeI*(motion(1,3)-vaPoint(1,i)*motion(1,5)+vaPoint(2,i)*motion(1,4))
	!加速度
!	pva(1,4,I)=-omeI*omeI*(motion(1,1)+vaPoint(3,i)*motion(1,5)-vaPoint(2,i)*motion(1,6))
!	pva(2,4,I)=-omeI*omeI*(motion(2,1)+vaPoint(3,i)*motion(2,5)-vaPoint(2,i)*motion(2,6))
!	pva(1,5,I)=-omeI*omeI*(motion(1,2)-vaPoint(3,i)*motion(1,4)+vaPoint(1,i)*motion(1,6))
!	pva(2,5,I)=-omeI*omeI*(motion(2,2)-vaPoint(3,i)*motion(2,4)+vaPoint(1,i)*motion(2,6))
!	pva(1,6,I)=-omeI*omeI*(motion(1,3)-vaPoint(1,i)*motion(1,5)+vaPoint(2,i)*motion(1,4))
!	pva(2,6,I)=-omeI*omeI*(motion(2,3)-vaPoint(1,i)*motion(2,5)+vaPoint(2,i)*motion(2,4))
	
    pva(1,4,I)=-omeI*omeI*(motion(1,1)+vaPoint(3,i)*motion(1,5)-vaPoint(2,i)*motion(1,6))&
	                +omeI*omeI*(-vaPoint(2,i)*(motion(1,4)*motion(1,5)-motion(2,4)*motion(2,5)))&
	                +omeI*omeI*(vaPoint(1,i)*(motion(1,5)*motion(1,5)-motion(2,5)*motion(2,5)))&
	                +omeI*omeI*(vaPoint(1,i)*(motion(1,6)*motion(1,6)-motion(2,6)*motion(2,6)))&
	                +omeI*omeI*(-vaPoint(3,i)*(motion(1,4)*motion(1,6)-motion(2,4)*motion(2,6))) !&
!	                -g0*motion(1,5)
	            	            
	pva(2,4,I)=-omeI*omeI*(motion(2,1)+vaPoint(3,i)*motion(2,5)-vaPoint(2,i)*motion(2,6))&
	                +omeI*omeI*(-vaPoint(2,i)*(motion(1,4)*motion(2,5)+motion(2,4)*motion(1,5)))&
	                +omeI*omeI*(vaPoint(1,i)*(motion(1,5)*motion(2,5)+motion(2,5)*motion(1,5)))&
	                +omeI*omeI*(vaPoint(1,i)*(motion(1,6)*motion(2,6)+motion(1,6)*motion(2,6)))&
	                +omeI*omeI*(-vaPoint(3,i)*(motion(1,4)*motion(2,6)+motion(2,4)*motion(1,6))) !&
	!                -g0*motion(2,5) 
	                
	pva(1,5,I)=-omeI*omeI*(motion(1,2)-vaPoint(3,i)*motion(1,4)+vaPoint(1,i)*motion(1,6))&
	                +omeI*omeI*(-vaPoint(3,i)*(motion(1,5)*motion(1,6)-motion(2,5)*motion(2,6)))&
	                +omeI*omeI*(vaPoint(2,i)*(motion(1,6)*motion(1,6)-motion(2,6)*motion(2,6)))&
	                +omeI*omeI*(vaPoint(2,i)*(motion(1,4)*motion(1,4)-motion(2,4)*motion(2,4)))&
	                +omeI*omeI*(-vaPoint(1,i)*(motion(1,4)*motion(1,5)-motion(2,4)*motion(2,5)))  !&
	 !           	+g0*motion(1,4)    
	            	        	
     pva(2,5,I)=-omeI*omeI*(motion(2,2)-vaPoint(3,i)*motion(2,4)+vaPoint(1,i)*motion(2,6))&            	
	                +omeI*omeI*(-vaPoint(3,i)*(motion(1,5)*motion(2,6)+motion(2,5)*motion(1,6)))&
	                +omeI*omeI*(vaPoint(2,i)*(motion(1,6)*motion(2,6)+motion(1,6)*motion(2,6)))&
	                +omeI*omeI*(vaPoint(2,i)*(motion(1,4)*motion(2,4)+motion(2,4)*motion(1,4)))&
	                +omeI*omeI*(-vaPoint(1,i)*(motion(1,4)*motion(2,5)+motion(2,4)*motion(1,5))) ! &
	!                +g0*motion(2,4)
	                
     pva(1,6,I)=-omeI*omeI*(motion(1,3)-vaPoint(1,i)*motion(1,5)+vaPoint(2,i)*motion(1,4))&
	                +omeI*omeI*(-vaPoint(1,i)*(motion(1,4)*motion(1,6)-motion(2,4)*motion(2,6)))&
	                +omeI*omeI*(vaPoint(3,i)*(motion(1,4)*motion(1,4)-motion(2,4)*motion(2,4)))&
	                +omeI*omeI*(vaPoint(3,i)*(motion(1,5)*motion(1,5)-motion(2,5)*motion(2,5)))&
	                +omeI*omeI*(-vaPoint(2,i)*(motion(1,5)*motion(1,6)-motion(2,5)*motion(2,6))) 
	            
	             
     pva(2,6,I)=-omeI*omeI*(motion(2,3)-vaPoint(1,i)*motion(2,5)+vaPoint(2,i)*motion(2,4))&
	                +omeI*omeI*(-vaPoint(1,i)*(motion(1,4)*motion(2,6)+motion(2,4)*motion(1,6)))&
	                +omeI*omeI*(vaPoint(3,i)*(motion(1,4)*motion(2,4)+motion(1,4)*motion(2,4)))&
	                +omeI*omeI*(vaPoint(3,i)*(motion(1,5)*motion(2,5)+motion(2,5)*motion(1,5)))&
	                +omeI*omeI*(-vaPoint(2,i)*(motion(1,5)*motion(2,6)+motion(2,5)*motion(1,6)))
enddo
if(IP>0)   then
    call va_putouts(omeI0,omeI,OMEI1,k0,pva,WavePhase,BETA)
endif
DEallocate(pva)

end subroutine shippointva
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------
!输出计算点速度加速度结果
!-----------------------------------------------------------------
subroutine va_putouts(omeI0,omeI,OMEI1,k0,pva,WavePhase,BETA)
implicit none

real(8),intent(in)::omeI0,k0,WavePhase,BETA,omeI,OMEI1
real(8),intent(in),dimension(1:2,1:6,1:NPva)::pva
real(8),dimension(1:2,1:6)::temp
integer(4)::i,j,k,jj
! real(8)::PI
! PI=4.0*atan(1.0)

if(TRIM(ADJUSTL(ansMethod))=="NO") then !谱分析值

		write(134,'(a,F5.1,a,f7.3,a,f7.3,a,e9.2)')"#HEAD=  ",BETA*180/PI,"  wl/sl= ",1.0/(k0*sl/(2*PI))," Ome0= ",omeI0," k0= ",k0
		write(134,'(a)')"#  ome0   ome"
		write(134,'(2f7.3)')omei0,OMEI1
		write(134,'(a,5x,a,3x,2(5x,a,3x),6(4x,a,4x,a,1x))')"#   Pt","x","y","z","Vx ","PHA","Vy ","PHA","Vz ","PHA","Ax ","PHA","Ay ","PHA","Az ","PHA"
		do i=1,NPva
			do k=1,6
				call comp_to_ampha(pva(1,K,I),pva(2,K,I),temp(1,K),temp(2,K))
			enddo
			write(134,'(i5,3f9.3,6(f8.4,f7.2))')I,(vaPoint(jj,i),jj=1,3),((temp(K,JJ),k=1,2),jj=1,6)
!			WRITE(134,*)
		enddo
else if(TRIM(ADJUSTL(ansMethod))=="YES") then   !确定值 !总垂向加速度

		write(19,'(a,F5.1,a,f7.3,a,f7.3,a,e9.2)')"#HEAD= ",BETA*180/PI,"  wl/sl= ",1.0/(k0*sl/(2*PI))," Ome0= ",omeI0," k0= ",k0
		write(19,'(a,5x,a,3x,2(5x,a,3x),5x,a,5x,a,8x,a,6x,a,7x,a,5x,a,7x,a,6x,a,8x,a,6x,a,7x,a,5x,a,1x)')"#   Pt","x","y","z","Vx ","PHA","Vy ","PHA","Vz ","PHA","Ax ","PHA","Ay ","PHA","Az ","PHA"

		do i=1,NPva
			do k=1,6
				call comp_to_ampha(pva(1,K,I),pva(2,K,I),temp(1,K),temp(2,K))
				temp(1,K)=temp(1,K)*cos((temp(2,K)+WavePhase)/180.0*PI)   !此处原程序有误temp(2,4)
			end	do
			write(19,'(i5,3f9.3,6(e12.4,f7.2))')I,(vaPoint(jj,i),jj=1,3),((temp(K,JJ),k=1,2),jj=1,6)
!			WRITE(19,*)
		end	do

endif
end subroutine va_putouts
!---------------------------------------------------------------------------------------------------







!-----------------------------------------------------------------
END MODULE POINTVA_MOD