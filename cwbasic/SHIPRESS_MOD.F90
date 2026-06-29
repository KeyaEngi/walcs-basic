MODULE SHIPPRESS_MOD

USE ENVIRONMENT_MOD,ONLY:ROU,G0,H,PI
USE CAL_MOD,ONLY:WATERDEPTH
USE MESH_MOD,ONLY:NB,SL,Factor,TarNum,Space
USE INOUTACCESS_MOD


IMPLICIT NONE

CONTAINS
!--------------------------------------------------


!******************************************************************************************
!程序功能：求解船舶表面压力
!输出参数：ps,pw,pd,pr,pall分别为静水压力,入射波压力,绕射压力,辐射压力和总压力的实部和虚部       
!******************************************************************************************
subroutine shippress(omeI,omeI1,UI,Beta,ampI,k0,xav0,phi,motion,ps,pw,pd,pr,pall,IP)
implicit none

real(8),intent(in)::omeI,omeI1,UI,Beta,k0,ampI
real(8),intent(in),dimension(1:3,1:NB)::xav0
real(8),intent(in),dimension(1:2,1:2,1:7,1:Factor,1:NB)::phi
real(8),intent(in),dimension(1:2,1:6)::motion
real(8),intent(out),dimension(1:2,1:Factor,1:NB)::ps,pw,pd,pr,pall
integer(4)::IP
integer(4)::i,j,L,count
real(8)::omeI0,temp1,temp2,temp3,temp4,x00(2)
real(8),dimension(1:2,1:2)::tphi !tphi(1,1:2)是入射波势对x导数的实虚部，tphi(2,1:2)是入射波势的实虚部
L=0

omeI0=omeI-k0*UI*cos(beta)
if(abs(omeI0)>1e-10)then
	if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
		temp1=k0*g0/(omeI0*cosh(k0*H))
		temp2=g0/(omeI0*cosh(k0*H))
	elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
		temp1=k0*g0/omeI0
		temp2=g0/omeI0
	endif
endif
do count=1,TarNum,Space
   L=L+1
   do i=1,NB
     select case(count)
       case(1)
              x00(1)=xav0(1,i)
              x00(2)=xav0(2,i)
        case(2)
              x00(1)=-xav0(1,i)
              x00(2)=xav0(2,i)       
       case(3)
             x00(1)=-xav0(1,i)
             x00(2)=-xav0(2,i)
        case(4)
            x00(1)=xav0(1,i)
            x00(2)=-xav0(2,i)
 end select
!--------------------------------------------------------------------------
	ps(1,L,I)=-rou*g0*(motion(1,3)+x00(2)*motion(1,4)-x00(1)*motion(1,5))
	ps(2,L,I)=-rou*g0*(motion(2,3)+x00(2)*motion(2,4)-x00(1)*motion(2,5))

!----入射波压力---------------------------------------------------------------------------
	temp3=cos(k0*(x00(1)*cos(beta)-x00(2)*sin(beta)))
	temp4=sin(k0*(x00(1)*cos(beta)-x00(2)*sin(beta)))
	if((trim(adjustl(WATERDEPTH))=="FINITE"))then  !有限水深
		tphi(1,1)=temp3*cosh(k0*(xav0(3,i)+H))*(-cos(beta)) !导数实部
		tphi(1,2)=temp4*cosh(k0*(xav0(3,i)+H))*(-cos(beta)) !导数虚部
	elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
		tphi(1,1)=temp3*exp(k0*xav0(3,i))*(-cos(beta))	!导数实部
		tphi(1,2)=temp4*exp(k0*xav0(3,i))*(-cos(beta))    !导数虚部
	endif
	if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
		tphi(2,1)=-cosh(k0*(xav0(3,i)+H))*temp4 !入射波势实部
		tphi(2,2)=cosh(k0*(xav0(3,i)+H))*temp3  !入射波势虚部
	elseif((trim(adjustl(WATERDEPTH))=="INFINITE"))then !无限水深
		tphi(2,1)=-exp(k0*xav0(3,i))*temp4 !入射波势实部
		tphi(2,2)=exp(k0*xav0(3,i))*temp3  !入射波势虚部
	endif
	tphi(1,:)=tphi(1,:)*temp1
	tphi(2,:)=tphi(2,:)*temp2
	pw(1,L,I)=rou*ampI*(omeI*tphi(2,2)+UI*tphi(1,1))
	pw(2,L,I)=-rou*ampI*(omeI*tphi(2,1)-UI*tphi(1,2))
!----绕射压力------------------------------------------------------------------------------
	pd(1,L,I)=rou*ampI*(omeI*phi(2,1,7,L,i)+UI*phi(1,2,7,L,i))
	pd(2,L,I)=-rou*ampI*(omeI*phi(1,1,7,L,i)-UI*phi(2,2,7,L,i))
!----辐射压力-----------------------------------------------------------------------------------
!有航速待修改
	pr(:,L,I)=0.0
	do j=1,6
		pr(1,L,I)=pr(1,L,I)+omeI*omeI*motion(1,J)*phi(1,1,j,L,i)-omeI*omeI*motion(2,J)*phi(2,1,j,L,i)&
		               -UI*omeI*motion(1,J)*phi(2,2,j,L,i)-UI*omeI*motion(2,J)*phi(1,2,j,L,i)
		pr(2,L,I)=pr(2,L,I)+omeI*omeI*motion(1,J)*phi(2,1,j,L,i)+omeI*omeI*motion(2,J)*phi(1,1,j,L,i)&
		               +UI*omeI*motion(1,J)*phi(1,2,j,L,i)-UI*omeI*motion(2,J)*phi(2,2,j,L,i)

	    if(j==5)then	
	         pr(1,L,I)=pr(1,L,I)+omeI*UI*motion(1,5)*phi(2,1,3,L,i)+omeI*UI*motion(2,5)*phi(1,1,3,L,i)&
			                +UI*UI*motion(1,5)*phi(1,2,3,L,i)-UI*UI*motion(2,5)*phi(2,2,3,L,i)
	         pr(2,L,I)=pr(2,L,I)-omeI*UI*motion(1,5)*phi(1,1,3,L,i)+omeI*UI*motion(2,5)*phi(2,1,3,L,i)&
			                +UI*UI*motion(1,5)*phi(2,2,3,L,i)+UI*UI*motion(2,5)*phi(1,2,3,L,i)
		ENDIF

		if(j==6)then	
	         pr(1,L,I)=pr(1,L,I)-omeI*UI*motion(1,6)*phi(2,1,2,L,i)-omeI*UI*motion(2,6)*phi(1,1,2,L,i)&
			                -UI*UI*motion(1,6)*phi(1,2,2,L,i)+UI*UI*motion(2,6)*phi(2,2,2,L,i)
	         pr(2,L,I)=pr(2,L,I)+omeI*UI*motion(1,6)*phi(1,1,2,L,i)-omeI*UI*motion(2,6)*phi(2,1,2,L,i)&
			                -UI*UI*motion(1,6)*phi(2,2,2,L,i)-UI*UI*motion(2,6)*phi(1,2,2,L,i)
		ENDIF
	enddo		
	pr(:,L,I)=rou*pr(:,L,I)
	pall(:,L,I)=ps(:,L,I)+pw(:,L,I)+pd(:,L,I)+pr(:,L,I)
!	write(10,*)L,I
   enddo
end do

if(IP>0)   then
    CALL p_putouts(BETA,omeI,omeI1,xav0,ps,pw,pd,pr,pall)
endif

end subroutine shippress
!---------------------------------------------------------------


!------------------------
!输出压力结果
!------------------------
subroutine p_putouts(BETA,omeI,omeI1,xav0,ps,pw,pd,pr,pall)
implicit none

real(8),intent(in)::omeI,BETA,omeI1
real(8),intent(in),dimension(1:3,1:NB)::xav0
real(8),intent(in),dimension(1:2,1:Factor,1:NB)::ps,pw,pd,pr,pall
real(8),dimension(1:5,1:2)::temp

real(8),dimension(1:3,1:NB)::x00
integer(4)::i,j,k,L,count
L=0
!PREOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.PRE'
!OPEN (3,FILE=PREOUTPUT)
!! real(8)::PI
!! PI=4.0*atan(1.0)
!write(3,"(a)")"#              [ WALCS Version 2.0 ]"
!write(3,*)""
!write(3,"(a)")"#              [ 压力响应计算结果输出文件 ]"
!write(3,*)""
write(3,'(a,f5.1,a,f7.3)')"HEAD= ",BETA*180/PI,"     Ome= ",omeI1
write(3,'(a,4x,a,5x,a,3x,2(5x,a,3x),5x,a,7x,a,6x,a,7x,a,6x,a,7x,a,6x,a,7x,a,7x,a,7x,a,1x)')"#","I","x","y","z","Psa","PHA","Pwa","PHA","Pda","PHA","Pra","PHA","Pa","PHA"
do count=1,TarNum,Space
L=L+1
select case(count)
       case(1)
                   x00=xav0
        case(2)
              x00(1,:)=-xav0(1,:)
              x00(2,:)=xav0(2,:)   
       case(3)
             x00(1,:)=-xav0(1,:)
             x00(2,:)=-xav0(2,:)
        case(4)
            x00(1,:)=xav0(1,:)
            x00(2,:)=-xav0(2,:)
 end select
do i=1,NB
	call comp_to_ampha(ps(1,L,I),ps(2,L,I),temp(1,1),temp(1,2))
	call comp_to_ampha(pw(1,L,I),pw(2,L,I),temp(2,1),temp(2,2))
	call comp_to_ampha(pd(1,L,I),pd(2,L,I),temp(3,1),temp(3,2))
	call comp_to_ampha(pr(1,L,I),pr(2,L,I),temp(4,1),temp(4,2))
	call comp_to_ampha(pall(1,L,I),pall(2,L,I),temp(5,1),temp(5,2))
	temp(:,1)=temp(:,1)*1000	
	write(3,'(i5,3f9.3,5(e12.4,f7.2))')i+(L-1)*NB,(X00(j,i),j=1,3),((temp(j,k),k=1,2),j=1,5)
enddo
end do
end subroutine p_putouts
!---------------------------------------------------------------------------------

!--------------------------------------------------
END MODULE SHIPPRESS_MOD