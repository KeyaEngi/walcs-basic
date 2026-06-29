MODULE POINTPRESS_MOD

USE INOUTACCESS_MOD
USE PVPOINT_MOD,ONLY:NPPOINT,PRESSPOINT
USE MESH_MOD,ONLY:DELTAX,ST,SL,TY,Factor,TarNum,Space,Ta
USE ENVIRONMENT_MOD,ONLY:H,ROU,G0,PI
USE CAL_MOD,ONLY:WATERDEPTH
USE PRINT_MOD,ONLY:ANSMETHOD



IMPLICIT NONE
PRIVATE
PUBLIC::READPPOINT,shippointpress

CONTAINS
!----------------------------------------------------------
SUBROUTINE READPPOINT(RATIOPTS)
IMPLICIT NONE

CHARACTER(LEN=300)::PTSINPUT
INTEGER::I,J,TEMP
REAL::RATIOPTS
real(8),dimension(1:3)::tempcoor

PTSINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.pts'
open(29,file=PTSINPUT,status='old')
call IUTMP(29)	
read(29,*)(tempcoor(i),i=1,3),ratioPTS !压力计算点所在坐标系在初始坐标系下的坐标
tempcoor(:)=tempcoor(:)/ratioPTS       !单位转化为米
read(29,*)NPpoint

allocate	(PressPoint(1:3,1:NPPoint))
! allocate	(ElmNum(1:NPPoint))
!do i=1,NPPoint
!		read(29,*)TEMP,(PressPoint(j,i),j=1,3)	
!		PressPoint(1,i)=PressPoint(1,i)/ratioPTS+tempcoor(1)-deltax           !压力计算点的坐标转化到随船平动坐标系下
!		PressPoint(2,i)=PressPoint(2,i)/ratioPTS+tempcoor(2)
!        PressPoint(3,i)=PressPoint(3,i)/ratioPTS+tempcoor(3)-ST        
!enddo
do i=1,NPPoint
		read(29,*)TEMP,(PressPoint(j,i),j=1,3)	
		PressPoint(1,i)=PressPoint(1,i)/ratioPTS+tempcoor(1)-deltax           !压力计算点的坐标转化到随船平动坐标系下
		PressPoint(2,i)=PressPoint(2,i)/ratioPTS+tempcoor(2)
        PressPoint(3,i)=PressPoint(3,i)/ratioPTS+tempcoor(3)-ST       
enddo
!----------纵倾坐标变换----------
do i=1,NPPoint
	PressPoint(1:3,i)=Matmul(Ty,PressPoint(1:3,i))
enddo

close(29)

END SUBROUTINE READPPOINT
!------------------------------------------------------------------


!************************************************************************************
!程序功能：求解压力计算点处船舶表面压力
!输入参数：PressPoint(1:3,1:NPPoint)为压力计算点坐标,phi_press(1:2,1:4,1:7,1:NPPoint)       
!************************************************************************************
subroutine shippointpress(omeI0,omeI,omeI1,U,Beta,ampI,k0,phi_press,motion,RATIOPTS,WAVEPHASE,IP)
implicit none

real(8),intent(in)::omeI0,omeI,omeI1,U,Beta,k0,ampI,WAVEPHASE
REAL,INTENT(IN)::RATIOPTS
real(8),intent(in),dimension(1:2,1:2,1:7,1:NPPoint)::phi_press
real(8),intent(in),dimension(1:2,1:6)::motion
real(8),ALLOCATABLE,dimension(:,:)::pps,ppw,ppd,ppr,ppall
integer(4)  ::  IP
integer(4)::i,j
real(8)::temp1,temp2,temp3,temp4
real(8),dimension(1:2,1:2)::tphi !tphi(1,1:2)是入射波势对x导数的实虚部，tphi(2,1:2)是入射波势的实虚部

ALLOCATE(pps(1:2,1:NPPOINT),ppw(1:2,1:NPPOINT),ppd(1:2,1:NPPOINT),ppr(1:2,1:NPPOINT),ppall(1:2,1:NPPOINT))

!omeI0=omeI-k0*U*cos(beta)
if(abs(omeI0)>1e-10)then
	if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
		temp1=k0*g0/(omeI0*cosh(k0*H))
		temp2=g0/(omeI0*cosh(k0*H))
	elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
		temp1=k0*g0/omeI0
		temp2=g0/omeI0
	endif
endif

do i=1,NPPoint
!--------------------------------------------------------------------------
	pps(1,I)=-rou*g0*(motion(1,3)+PressPoint(2,i)*motion(1,4)-PressPoint(1,i)*motion(1,5))
	pps(2,I)=-rou*g0*(motion(2,3)+PressPoint(2,i)*motion(2,4)-PressPoint(1,i)*motion(2,5))
!----入射波压力---------------------------------------------------------------------------
	temp3=cos(k0*(PressPoint(1,i)*cos(beta)-PressPoint(2,i)*sin(beta)))
	temp4=sin(k0*(PressPoint(1,i)*cos(beta)-PressPoint(2,i)*sin(beta)))
	if((trim(adjustl(WATERDEPTH))=="FINITE"))then  !有限水深
		tphi(1,1)=temp3*cosh(k0*(PressPoint(3,i)+H))*(-cos(beta)) !导数实部
		tphi(1,2)=temp4*cosh(k0*(PressPoint(3,i)+H))*(-cos(beta)) !导数虚部
		tphi(2,1)=-cosh(k0*(PressPoint(3,i)+H))*temp4 !入射波势实部
		tphi(2,2)=cosh(k0*(PressPoint(3,i)+H))*temp3  !入射波势虚部
	elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
		tphi(1,1)=temp3*exp(k0*PressPoint(3,i))*(-cos(beta))	!导数实部
		tphi(1,2)=temp4*exp(k0*PressPoint(3,i))*(-cos(beta))    !导数虚部
		tphi(2,1)=-exp(k0*PressPoint(3,i))*temp4 !入射波势实部
		tphi(2,2)=exp(k0*PressPoint(3,i))*temp3  !入射波势虚部
	endif

	tphi(1,:)=tphi(1,:)*temp1
	tphi(2,:)=tphi(2,:)*temp2
	ppw(1,I)=rou*ampI*(omeI*tphi(2,2)+U*tphi(1,1))
	ppw(2,I)=-rou*ampI*(omeI*tphi(2,1)-U*tphi(1,2))
!----绕射压力------------------------------------------------------------------------------
	ppd(1,I)=rou*ampI*(omeI*phi_press(2,1,7,i)+U*phi_press(1,2,7,i))
	ppd(2,I)=-rou*ampI*(omeI*phi_press(1,1,7,i)-U*phi_press(2,2,7,i))
!----辐射压力-----------------------------------------------------------------------------------
!有航速待修改
	ppr(:,I)=0.0
	do j=1,6
		ppr(1,I)=ppr(1,I)+omeI*omeI*motion(1,J)*phi_press(1,1,j,i)-omeI*omeI*motion(2,J)*phi_press(2,1,j,i)&
		               -U*omeI*motion(1,J)*phi_press(2,2,j,i)-U*omeI*motion(2,J)*phi_press(1,2,j,i)
		ppr(2,I)=ppr(2,I)+omeI*omeI*motion(1,J)*phi_press(2,1,j,i)+omeI*omeI*motion(2,J)*phi_press(1,1,j,i)&
		               +U*omeI*motion(1,J)*phi_press(1,2,j,i)-U*omeI*motion(2,J)*phi_press(2,2,j,i)
		if(j==5)then
			ppr(1,I)=ppr(1,I)+omeI*U*motion(1,J)*phi_press(2,1,3,i)+omeI*U*motion(2,J)*phi_press(1,1,3,i)&
			               +U*U*motion(1,J)*phi_press(1,2,3,i)-U*U*motion(2,J)*phi_press(2,2,3,i)
			ppr(2,I)=ppr(2,I)-omeI*U*motion(1,J)*phi_press(1,1,3,i)+omeI*U*motion(2,J)*phi_press(2,1,3,i)&
			               +U*U*motion(1,J)*phi_press(2,2,3,i)+U*U*motion(2,J)*phi_press(1,2,3,i)
		endif
		if(j==6)then
			ppr(1,I)=ppr(1,I)-omeI*U*motion(1,J)*phi_press(2,1,2,i)-omeI*U*motion(2,J)*phi_press(1,1,2,i)&
			               -U*U*motion(1,J)*phi_press(1,2,2,i)+U*U*motion(2,J)*phi_press(2,2,2,i)
			ppr(2,I)=ppr(2,I)+omeI*U*motion(1,J)*phi_press(1,1,2,i)-omeI*U*motion(2,J)*phi_press(2,1,2,i)&
			               -U*U*motion(1,J)*phi_press(2,2,2,i)-U*U*motion(2,J)*phi_press(1,2,2,i)
		endif
	enddo
	ppr(:,I)=rou*ppr(:,I)
	ppall(:,I)=pps(:,I)+ppw(:,I)+ppd(:,I)+ppr(:,I)
enddo

if(IP>0)   then
    CALL pp_putouts(BETA,k0,omeI0,omeI,omeI1,ratioPTS,pps,ppw,ppd,ppr,ppall,WavePhase)
endif

DEALLOCATE(pps,ppw,ppd,ppr,ppall)

end subroutine shippointpress
!----------------------------------------------------------------------------------------------------------



!------------------------------------------------------------------
!输出压力结果
!------------------------------------------------------------------
subroutine pp_putouts(BETA,k0,omeI0,omeI,omeI1,ratioPTS,pps,ppw,ppd,ppr,ppall,WavePhase)
implicit none

real,intent(in)::ratioPTS
! integer(4),dimension(1:NPoint)::ElmNum
! integer(4),dimension(1:1000000)	::	WetEleID
real(8),intent(in)::BETA,k0,WavePhase,omeI,omeI0,omeI1
real(8),intent(in),dimension(1:2,1:NPPoint)::pps,ppw,ppd,ppr,ppall
real(8),dimension(1:2,1:5)::temp
integer(4)::i,j,k,jj
real(8)::pa100
! PI=4.0*atan(1.0)
if(trim(adjustl(ANSMETHOD))=="NO")then !谱分析值
! 	do j=1,NP
		write(133,'(a,F5.1,a,f6.3,a,f6.3,a,F6.3)')"#HEAD= ",BETA*180.0/PI,"  wl/sl=",1.0/(k0*sl/(2*PI)),"  Ome=",omeI,"  k0=",k0
		write(133,'(a)')"# ome0  ome"
	    write(133,'(2f7.3)')omei0,omeI1
		write(133,'(a,5x,a,8x,a,9x,a,7x,a,7x,a,6x,a,7x,a,6x,a,7x,a,6x,a,7x,a,7x,a,7x,a,1x)')"#    Pt","x","y","z","Psa","PHA","Pwa","PHA","Pda","PHA","Pra","PHA","Pa","PHA"
		do i=1,NPPoint
			call comp_to_ampha(pps(1,I),pps(2,I),temp(1,1),temp(2,1))
			call comp_to_ampha(ppw(1,I),ppw(2,I),temp(1,2),temp(2,2))
			call comp_to_ampha(ppd(1,I),ppd(2,I),temp(1,3),temp(2,3))
			call comp_to_ampha(ppr(1,I),ppr(2,I),temp(1,4),temp(2,4))
			call comp_to_ampha(ppall(1,I),ppall(2,I),temp(1,5),temp(2,5))
			if(PressPoint(3,i)>0.0)temp(:,:)=0.0
			temp(1,:)=temp(1,:)*1000/(ratioPTS*ratioPTS)
			write(133,'(i7,3f9.3,5(e12.4,f7.2))')I,(PressPoint(jj,i),jj=1,3),((temp(K,JJ),k=1,2),jj=1,5)
		enddo
		WRITE(133,*)
! 	enddo
elseif(trim(adjustl(ANSMETHOD))=="YES") then !确定值 !总压力
! 	do j=1,NP
		write(18,'(a,F5.1,a,f6.3,a,f6.3,a,F6.3)')"#HEAD= ",BETA*180.0/PI,"  wl/sl=",1.0/(k0*sl/(2*PI)),"  Ome=",omeI1,"  k0=",k0
		write(18,'(a,5x,a,8x,a,9x,a,7x,a,7x,a,6x,a,7x,a,6x,a,7x,a,6x,a,7x,a,7x,a,7x,a,1x)')"#    Pt","x","y","z","Psa","PHA","Pwa","PHA","Pda","PHA","Pra","PHA","Pa","PHA","Pa100"
		do i=1,NPPoint
			call comp_to_ampha(pps(1,I),pps(2,I),temp(1,1),temp(2,1))
			call comp_to_ampha(ppw(1,I),ppw(2,I),temp(1,2),temp(2,2))
			call comp_to_ampha(ppd(1,I),ppd(2,I),temp(1,3),temp(2,3))
			call comp_to_ampha(ppr(1,I),ppr(2,I),temp(1,4),temp(2,4))
			call comp_to_ampha(ppall(1,I),ppall(2,I),temp(1,5),temp(2,5))
			pa100=temp(1,5)*cos((temp(2,5)+WavePhase)/180.0*PI)
			pa100=pa100*1000/(ratioPTS*ratioPTS)
			if(PressPoint(3,i)>0.0)pa100=0.0
			write(18,'(i7,3f9.3,5(e12.4,f7.2))')I,(PressPoint(jj,i),jj=1,3),((temp(jj,k),JJ=1,2),K=1,5)
		enddo
		WRITE(18,*)
! 	enddo
! elseif(inpr8==2) then !确定值 !总压力!增加了静水压力部分
! 	do j=1,NP
! 		write(FileNum,'(a,i6,a,f7.3,a,f7.3,a,e9.2)')" IP= ",j,"  wl/sl= ",1.0/(k0(j)*sl/(2*PI))," Ome0= ",ome0(j)," k0= ",k0(j)
! 		write(FileNum,*)"NPoint"
! 		write(FileNum,'(i10)')NPoint
! 		write(FileNum,'(4x,a,4x,a,3x)')"Pt","Pa100"
! 		do i=1,NPoint
! 			call comp_to_ampha(ppall(j,i,1),ppall(j,i,2),temp(5,1),temp(5,2))
! 			pa100=amp(j)*temp(5,1)*cos((temp(5,2)+WavePhase(j))/180.0*PI)-rou*g0*PressPoint(3,i)
! 			pa100=pa100*1000/(ratio*ratio)
! 			if(PressPoint(3,i)>0.0)pa100=0.0
! 			write(FileNum,'(i10,e12.4)')ElmNum(i),pa100
! 		enddo
! 	enddo
! elseif(inpr8==3) then
! 	do j=1,NP
! 		write(FileNum,'(a,i6,a,f7.3,a,f7.3,a,e9.2)')" IP= ",j,"  wl/sl= ",1.0/(k0(j)*sl/(2*PI))," Ome0= ",ome0(j)," k0= ",k0(j)
! 		write(FileNum,*)"NPoint"
! 		write(FileNum,'(i6)')NPoint
! 		write(FileNum,'(4x,a,4x,a,3x)')"Pt","Pa100"
! 		do i=1,NPoint
! 			call comp_to_ampha(ppall(j,i,1),ppall(j,i,2),temp(5,1),temp(5,2))
! 			pa100=amp(j)*temp(5,1)*cos((temp(5,2)+WavePhase(j))/180.0*PI)-rou*g0*PressPoint(3,i)
! 			pa100=pa100*1000/(ratio*ratio)
! 			if(PressPoint(3,i)>0.0)pa100=0.0
! 			write(FileNum,'(i10,a,e12.4)')ElmNum(i),", P, ",pa100
! 		enddo
! 	enddo
endif

end subroutine pp_putouts
!-----------------------------------------------------------------------------------------






!----------------------------------------------------------
END MODULE POINTPRESS_MOD