!--------------------------------------------------------------------------------------------
!模块功能：借助于浮体结构的两个几何对称面(XOZ，YOZ)的关系，计算与面元有关的几何量以及基本源在面元上的单层势及诱导速度
!输出参数：SRPHI(1:NPanel,1:Npanel,1:4,1:Factor)为第一下标(场点)在第二下标(源点)下的想x,y,z诱导速度及单层势，最后的1:Factor表示在不同象限的映射点值。
!                       DRPHI1:NPanel,1:Npanel,1:Factor)为上述诱导速度在场点面元的法向投影
!模块修改：王川 2012.07.26
!--------------------------------------------------------------------------------------------
MODULE BASICSOURCE_MOD
USE INOUTACCESS_MOD
USE MESH_MOD,ONLY:XAV,E,NPANEL,XQ,EA,NB,Factor,TarNum,Space
USE ENVIRONMENT_MOD,ONLY:H,PI
USE CAL_MOD,ONLY:WATERDEPTH
USE PVPOINT_MOD,ONLY:NPPOINT,PRESSPOINT

IMPLICIT NONE 
PRIVATE
PUBLIC::GR,GR_P
CONTAINS
!--------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!计算基本源在面元上的单层势及诱导速度
!--------------------------------------------------------------------------------------------
subroutine GR(SRPHI,DRPHI)
implicit none
real(8),intent(out),dimension(1:NPanel,1:NPanel,1:4,1:Factor)::SRPHI	!诱导速度及单层势
real(8),intent(out),dimension(1:NPanel,1:NPanel,1:Factor)::DRPHI	  !诱导速度法向投影
real(8),dimension(1:NPanel,1:NPanel)::AIJ
real(8),dimension(1:3)::xp
real(8),dimension(1:2)::SR
real(8),dimension(1:3,1:2)::SS
real(8),allocatable::x00(:,:),ee(:,:,:),xqq(:,:,:)
CHARACTER(LEN=500)::BASS
real(8)::s0,s1,s2,s3,Gx,Gy,Gz
real(8)::eps
integer(4)::i,j,k,jj,ii,L,count
allocate(x00(1:3,1:Npanel),ee(1:3,1:3,1:Npanel),xqq(1:2,1:4,1:Npanel))
eps=1.0e-6
SRPHI=0. ;DRPHI=0.
L=0

do count=1,TarNum,Space
        select case(count)
               case(1)
                      x00=xav
                      ee=e
                      xqq=xq
                case(2)
                      x00(1,:)=-xav(1,:)
                      x00(2,:)=xav(2,:)   
                      ee(1,3,:)=-e(1,3,:)
                      ee(2,3,:)=e(2,3,:)
                      ee(3,3,:)=e(3,3,:)    
                      ee(1,1,:)=-e(1,1,:)
                      ee(2,1,:)=e(2,1,:)
                      ee(3,1,:)=e(3,1,:)
                      ee(1,2,:)=e(1,2,:)
                      ee(2,2,:)=-e(2,2,:)
                      ee(3,2,:)=-e(3,2,:)
                      xqq(1,:,:)=xq(1,:,:)
                      xqq(2,:,:)=-xq(2,:,:)
               case(3)
                      x00(1,:)=-xav(1,:)
                      x00(2,:)=-xav(2,:)
                      ee(1,3,:)=-e(1,3,:)
                      ee(2,3,:)=-e(2,3,:)
                      ee(3,3,:)=e(3,3,:)
                      ee(1,1,:)=-e(1,1,:)
                      ee(2,1,:)=-e(2,1,:)
                      ee(3,1,:)=e(3,1,:)
                      ee(1,2,:)=-e(1,2,:)
                      ee(2,2,:)=-e(2,2,:)
                      ee(3,2,:)=e(3,2,:)
                      xqq(1,:,:)=xq(1,:,:)
                      xqq(2,:,:)=xq(2,:,:)
                case(4)
                      x00(1,:)=xav(1,:)
                      x00(2,:)=-xav(2,:)
                      ee(1,3,:)=e(1,3,:)
                      ee(2,3,:)=-e(2,3,:)
                      ee(3,3,:)=e(3,3,:)
                      ee(1,1,:)=e(1,1,:)
                      ee(2,1,:)=-e(2,1,:)
                      ee(3,1,:)=e(3,1,:)
                      ee(1,2,:)=-e(1,2,:)
                      ee(2,2,:)=e(2,2,:)
                      ee(3,2,:)=-e(3,2,:)
                      xqq(1,:,:)=xq(1,:,:)
                      xqq(2,:,:)=-xq(2,:,:)
         end select
                      L=L+1
do i=1,NPanel  !p
      do j=1,NPanel  !q
	    	do jj=1,2
	 	    	call pan5(x00(:,j),xav(:,i),jj,ee(:,:,j),xp)
			    call induce(s0,s1,s2,s3,xqq(:,:,j),xp,1)
                     if(count==2.or.count==4)then
                             s0=-s0
                             s1=-s1
                             s2=-s2
                             s3=-s3
                     end if
		             if((i==j).and.(jj==1).and.(count==1))then
			                  s1=0.0
			                  s2=0.0
			                  s3=2.0*PI
		             endif
		             if((i==j).and.(jj==2).and.(count==1).and.(i>NB).and.(trim(adjustl(WATERDEPTH))=="INFINITE"))then
			                  s1=0.0
			                  s2=0.0
			                  s3=0.0
		              endif
!--------------------------------------------------------------------------------------------
!（矢量转换）变换为全局坐标
!--------------------------------------------------------------------------------------------
				    SS(1,jj)=s1*ee(1,1,j)+s2*ee(1,2,j)+s3*ee(1,3,j)
				    SS(2,jj)=s1*ee(2,1,j)+s2*ee(2,2,j)+s3*ee(2,3,j)
				    SS(3,jj)=s1*ee(3,1,j)+s2*ee(3,2,j)+s3*ee(3,3,j)
!--------------------------------------------------------------------------------------------			
!标量无需转换坐标
!--------------------------------------------------------------------------------------------
			        SR(jj)=s0
		    enddo
SRPHI(i,j,4,L)=SR(1)+SR(2)
if(i==j.and.count==1)then
			Gx=SS(1,2)
			Gy=SS(2,2)
			Gz=-SS(3,2)
			DRPHI(i,j,L)=e(1,3,i)*Gx+e(2,3,i)*Gy+e(3,3,i)*Gz
			DRPHI(i,j,L)=2*pi+DRPHI(i,j,L)
else
			Gx=SS(1,1)+SS(1,2)
			Gy=SS(2,1)+SS(2,2)
			Gz=SS(3,1)-SS(3,2)
!--------------------------------------------------------------------------------------------
!求诱导速度法向投影
!--------------------------------------------------------------------------------------------
		    DRPHI(i,j,L)=e(1,3,i)*Gx+e(2,3,i)*Gy+e(3,3,i)*Gz
end if
		    SRPHI(i,j,1,L)=Gx
	     	SRPHI(i,j,2,L)=Gy
		    SRPHI(i,j,3,L)=Gz
		enddo
enddo
end do	
!    BASS=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.BASS'
!    OPEN (30021,FILE=BASS)
!do i=1,1  !p
!if(XAV(2,i)>=0)then
!      do j=1,NPanel  !q
!      if(XAV(2,J)>=0)then
!WRITE(30021,"(9F15.5)")XAV(1,I),XAV(2,I),XAV(3,I),XAV(1,J),XAV(2,J),XAV(3,J),SRPHI(i,j,1,1),SRPHI(i,j,2,1),SRPHI(i,j,3,1)
!end if
!END DO
!end if
!END DO
!
!do i=1,1 !p
!if(XAV(2,i)>=0)then
!      do j=1,NPanel  !q
!      if(XAV(2,J)<0)then
!WRITE(30021,"(9F15.5)")XAV(1,I),XAV(2,I),XAV(3,I),XAV(1,J),-XAV(2,J),XAV(3,J),SRPHI(i,j,1,1),SRPHI(i,j,2,1),SRPHI(i,j,3,1)
!end if
!END DO
!end if
!END DO
deallocate(x00,ee,xqq)
end subroutine GR
!--------------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------------
!程序功能：计算基本源对任意空间点(压力计算点)的单层势及诱导速度
!输出参数：s_press为压力计算点在各面元下的诱导速度及单层势
!调用子程序：pan5,induce
!程序编制：张海彬    时间：2002年4月16日        
!--------------------------------------------------------------------------------------------
subroutine GR_P(s_press)
implicit none
real(8),intent(out),dimension(1:NPPoint,1:NPanel,1:4,1:4)::s_press	!诱导速度及单层势
real(8),dimension(1:3)::xp
real(8),dimension(1:2)::SR
real(8),dimension(1:3,1:2)::SS
real(8)::s0,s1,s2,s3,Gx,Gy,Gz
real(8)::eps,dist
real(8),allocatable::x00(:,:),ee(:,:,:),xqq(:,:,:)

integer(4)::i,j,k,jj,ii,L,count
allocate(x00(1:3,1:Npanel),ee(1:3,1:3,Npanel),xqq(1:2,1:4,1:Npanel))
L=0
eps=1.0e-6

do count=1,TarNum,Space
        select case(count)
               case(1)
                      x00=xav
                      ee=e
                      xqq=xq
                case(2)
                      x00(1,:)=-xav(1,:)
                      x00(2,:)=xav(2,:)   
                      ee(1,3,:)=-e(1,3,:)
                      ee(2,3,:)=e(2,3,:)
                      ee(3,3,:)=e(3,3,:)    
                      ee(1,1,:)=-e(1,1,:)
                      ee(2,1,:)=e(2,1,:)
                      ee(3,1,:)=e(3,1,:)
                      ee(1,2,:)=e(1,2,:)
                      ee(2,2,:)=-e(2,2,:)
                      ee(3,2,:)=-e(3,2,:)
                      xqq(1,:,:)=xq(1,:,:)
                      xqq(2,:,:)=-xq(2,:,:)
               case(3)
                      x00(1,:)=-xav(1,:)
                      x00(2,:)=-xav(2,:)
                      ee(1,3,:)=-e(1,3,:)
                      ee(2,3,:)=-e(2,3,:)
                      ee(3,3,:)=e(3,3,:)
                      ee(1,1,:)=-e(1,1,:)
                      ee(2,1,:)=-e(2,1,:)
                      ee(3,1,:)=e(3,1,:)
                      ee(1,2,:)=-e(1,2,:)
                      ee(2,2,:)=-e(2,2,:)
                      ee(3,2,:)=e(3,2,:)
                      xqq(1,:,:)=xq(1,:,:)
                      xqq(2,:,:)=xq(2,:,:)
                case(4)
                      x00(1,:)=xav(1,:)
                      x00(2,:)=-xav(2,:)
                      ee(1,3,:)=e(1,3,:)
                      ee(2,3,:)=-e(2,3,:)
                      ee(3,3,:)=e(3,3,:)
                      ee(1,1,:)=e(1,1,:)
                      ee(2,1,:)=-e(2,1,:)
                      ee(3,1,:)=e(3,1,:)
                      ee(1,2,:)=-e(1,2,:)
                      ee(2,2,:)=e(2,2,:)
                      ee(3,2,:)=-e(3,2,:)
                      xqq(1,:,:)=xq(1,:,:)
                      xqq(2,:,:)=-xq(2,:,:)
         end select
                      L=L+1
do i=1,NPPoint  !p
!采用对称性的坐标变换   edited by 王川 2012.07.28
	do j=1,NPanel  !q

		dist=sqrt(sum((Presspoint(:,i)-x00(:,j))**2))
		do jj=1,2
			call pan5(x00(:,j),PressPoint(:,i),jj,ee(:,:,j),xp)
			call induce(s0,s1,s2,s3,xqq(:,:,j),xp,1)
         if(count==2.or.count==4)then
         s0=-s0
         s1=-s1
         s2=-s2
         s3=-s3
         end if
			if((dist<eps).and.(jj==1))then
				s1=0.0
				s2=0.0
				s3=2.0*PI
			endif
			do k=1,3
!--------------------------------------------------------------------------------------------			
!变换为全局坐标
!--------------------------------------------------------------------------------------------
			SS(k,jj)=s1*ee(k,1,j)+s2*ee(k,2,j)+s3*ee(k,3,j)
			enddo
			SR(jj)=s0
		enddo
			s_press(i,j,4,L)=SR(1)+SR(2)
			Gx=SS(1,1)+SS(1,2)
			Gy=SS(2,1)+SS(2,2)
			Gz=SS(3,1)-SS(3,2)			
		
			s_press(i,j,1,L)=Gx
		    s_press(i,j,2,L)=Gy
		    s_press(i,j,3,L)=Gz		
	enddo
enddo
end do
deallocate(x00,ee,xqq)
end subroutine GR_P
!--------------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------------
!pan5A计算场点p(一象限)与源点q(一象限)(flagq=1)或映像点(一象限)(flagq=2)间的距离在源点q局部坐标系下的投影
!--------------------------------------------------------------------------------------------
subroutine pan5(qcoor,pcoor,flagq,eq,xp)
implicit none
real(8),intent(in),dimension(1:3)::qcoor,pcoor
integer(4),intent(in)::flagq
real(8),intent(in),dimension(1:3,1:3)::eq
real(8),intent(out),dimension(1:3)::xp
real(8),dimension(1:3)::PQR
integer(4)::i
                     PQR(1)=pcoor(1)-qcoor(1)
                     if(flagq==1)then
                        	PQR(2:3)=pcoor(2:3)-qcoor(2:3)
                     elseif(flagq==2)then
	                 if((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
	                     	PQR(2)=pcoor(2)-qcoor(2)
	                    	PQR(3)=-pcoor(3)-qcoor(3)
	                 elseif((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
		                    PQR(2)=pcoor(2)-qcoor(2)
                     		PQR(3)=-2*H-pcoor(3)-qcoor(3)
                 	endif
                    endif
!坐标变换
do i=1,3
	xp(i)=dot_product(eq(:,i),PQR)
enddo
end subroutine pan5
!--------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------
!计算1/rpq在面元上的单层势s及诱导速度sx,sy,sz(flagd=1)或双层势诱导速度dx,dy,dz(flagd=2)
!--------------------------------------------------------------------------------------------
subroutine induce(s,s1,s2,s3,xqq,xp,flagd)
implicit none
real(8),intent(in),dimension(1:2,1:4)::xqq !源点
real(8),intent(in),dimension(1:3)::xp	  !场点
integer(4),intent(in)::flagd			  !控制参数
real(8),intent(out)::s,s1,s2,s3
real(8),dimension(1:2,1:4)::w
real(8),dimension(1:4)::L,r,AF,BT,M,N
real(8)::D,D1,D2,D3,eps,A,B
real(8)::w1,w2,w3,w4,w5,ww,x1,x2,x3,x4
integer(4)::i,i1
s=0.0;s1=0.0;s2=0.0;s3=0.0;D=0.0;D1=0.0;D2=0.0;D3=0.0
AF=0.0;BT=0.0;M=0.0;N=0.0
eps=1.0e-6
do i=1,4
	i1=i+1
	if(i==4)i1=1
	w(:,i)=xqq(:,i1)-xqq(:,i)
	L(i)=sqrt(sum(w(:,i)**2))
	w1=xp(1)-xqq(1,i)
	w2=xp(2)-xqq(2,i)
	w3=xp(3)
	r(i)=sqrt(w1*w1+w2*w2+w3*w3)
	if(abs(L(i))<eps) cycle
	AF(i)=w(1,i)/L(i)
	BT(i)=w(2,i)/L(i)
	M(i)=-AF(i)*w1-BT(i)*w2
	N(i)=AF(i)*w2-BT(i)*w1
enddo
do i=1,4
	if(abs(L(i))<eps) cycle
	i1=i+1
	if(i==4) i1=1
	w4=r(i)+r(i1)-L(i)
	w5=r(i)+r(i1)+L(i)
	if(abs(w4)<eps) cycle !stop '!场点在面元一边上' 不计该边贡献
	s=s+log(w5/w4)*N(i)
	s1=s1+log(w5/w4)*(-BT(i))
	s2=s2+log(w5/w4)*AF(i)
	ww=L(i)*(r(i)+r(i1))/(r(i)*r(i1)*w4*w5)
	D1=D1-2*BT(i)*xp(3)*ww
	D2=D2+2*AF(i)*xp(3)*ww
	D3=D3-2*N(i)*ww
enddo
if(abs(xp(3))>eps)then !如果场点与四边形共面且不在某边上，则法向速度贡献为0
	do i=1,4
		i1=i+1
		if(i==4) i1=1
		A=xqq(1,i1)-xqq(1,i)
		if(abs(A)<eps) cycle
		B=(xqq(2,i1)-xqq(2,i))/A
		x1=xqq(1,i)-xp(1)
		x2=xqq(2,i)-xp(2)
		x3=xqq(1,i1)-xp(1)
		x4=xqq(2,i1)-xp(2)
		s3=s3+atan((B*(x1*x1+xp(3)**2)-x1*x2)/(xp(3)*r(i)))&
		     -atan((B*(x3*x3+xp(3)**2)-x3*x4)/(xp(3)*r(i1)))
	    enddo
	s=s+xp(3)*s3
endif
if(flagd==2) then
	s1=D1
	s2=D2
	s3=D3
	s=D
endif
end subroutine induce
!--------------------------------------------------------------------------------------------



!--------------------------------------------------------------------------------------------
END MODULE BASICSOURCE_MOD