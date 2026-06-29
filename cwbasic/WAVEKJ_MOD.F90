!********************************************************************************
!模块功能：计算浅水色散方程的实数根(k0)和虚数根(kj)
!输出参数：wkj(1)为实数根(k0)，wkj(2:wkjNum)为虚数根(kj)
!********************************************************************************
MODULE WAVEKJ_MOD

USE ENVIRONMENT_MOD,ONLY:PI

IMPLICIT NONE
PRIVATE
PUBLIC::WAVEKJ

CONTAINS 
!----------------------------------------------------------------------
subroutine wavekj(wkjNum,wkj,hd,pnu)
implicit none
integer(4),intent(in)::wkjNum
real(8),intent(in)::hd,pnu
real(8),intent(out),dimension(1:wkjNum)::wkj
real(8)::ph,eps,aa,bb,cc,c,wor
integer(4)::i
! PI=4.0*atan(1.0)
eps=1e-10
ph=pnu*hd
if(ph>150)then
	wkj(1)=pnu
else
	i=0
	do 
		i=i+1
		wor=fmo(i*pnu,hd,pnu,1)
		if(wor>0.0)then
			aa=(i-1)*pnu
			bb=aa+pnu
			cc=(aa+bb)/2
			exit
		 endif
	enddo
	call ws(aa,bb,eps,wkj(1),hd,pnu,1)
endif
do i=2,wkjNum
	c=abs(wkj(i-1)*hd-(i-2)*PI)
	if(c>eps)then
		aa=((0.5+eps)*PI+(i-2)*PI)/hd
		bb=aa+0.5*PI/hd
		cc=(aa+bb)/2
		call ws(aa,bb,eps,wkj(i),hd,pnu,2)
	else
		wkj(i)=(i-1)*PI/hd
	endif
enddo
end subroutine wavekj

!------------------------------------------------------
subroutine ws(aa,bb,eps,x,hd,pnu,flag)
implicit none
real(8),intent(in)::aa,bb,eps,hd,pnu
integer(4),intent(in)::flag
real(8),intent(out)::x
real(8)::x1,x2,f,fa,fb
fa=fmo(aa,hd,pnu,flag)
fb=fmo(bb,hd,pnu,flag)
if(fa<=0.0.and.fb>0.0)then
	x1=aa
	x2=bb
elseif(fa>0.0.and.fb<0.0)then
	x1=bb
	x2=aa
else
	stop '根隔离失败!'
endif
do 
	x=(x1+x2)/2
	f=fmo(x,hd,pnu,flag)
	if(f<0.0)then
		x1=x
	else
		x2=x
	endif
	if(abs(x1*hd-x2*hd)<eps)then
		x=(x1+x2)/2
		exit
	endif
enddo
end subroutine ws				

!---------------------------------------------------
!浅水色散方程
function fmo(x,hd,pnu,flag) result(fmo_result)
implicit none
real(8),intent(in)::x,hd,pnu
integer(4),intent(in)::flag
real(8)::fmo_result,t
if(flag==1)then
	t=tanh(hd*x)
	fmo_result=x*tanh(hd*x)-pnu
elseif(flag==2)then
	fmo_result=x*tan(hd*x)+pnu
endif
end function fmo
!---------------------------------------------------

END MODULE WAVEKJ_MOD