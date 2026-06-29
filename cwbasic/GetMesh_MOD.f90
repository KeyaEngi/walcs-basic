MODULE GETMESH_MOD

USE MESH_MOD
USE INOUTACCESS_MOD
USE CAL_MOD,ONLY:POINTPRESS,FLOATATION
USE ENVIRONMENT_MOD,ONLY:PI

IMPLICIT NONE
PRIVATE
PUBLIC::READMESH


CONTAINS

!---------------------------读取水动力网格(ProjName.gdf)----------------------------------------
SUBROUTINE READMESH()

IMPLICIT NONE

INTEGER::NStern,NBow
real(8)::TrimAng,temp !纵倾角
CHARACTER(LEN=300)::GDFINPUT, MEDINPUT
INTEGER::I,J,K
INTEGER,ALLOCATABLE,DIMENSION(:)::IEC !面元类型

MEDINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.MED' !建模方法格式文件 .MED
open(32,file=MEDINPUT,STATUS="OLD")
call IUTMP(32)
READ(32,"(A)") MESHTYPE
CLOSE(32)

GDFINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.FSA'
open(40,file=GDFINPUT,STATUS="OLD")
call IUTMP(40)

if(trim(adjustl(MESHTYPE))=="PAM")	then !PAM=参数化建模
    read(40,*) SYMTYPE		!对称性symmetry类型  SYMTYPE=NO or X or XY						
	read(40,*)    TrimAng !纵倾角
	read(40,*)	NPanel,NB,NStern,NBow,ArcNum !算上内外网格的总网格数、湿表面网格数、尾斜线数据点数、首斜线数据点数、首尾水线等分弧长数目
	POINTPRESS=0

else if(trim(adjustl(MESHTYPE))=="BDF")then !bdf有限元模型导入
    read(40,*) SYMTYPE  !读入对称信息

	read(40,*)    TrimAng !读入纵倾角
	read(40,*)	NPanel,NB	!浮体网格数目，湿表面网格数目
	POINTPRESS=1		    			
endif
if(trim(adjustl(SYMTYPE))=="NO")then
Factor=1
TarNum=1
Space=1
else if(trim(adjustl(SYMTYPE))=="X")then
Factor=2
TarNum=4
Space=3
ArcNum=2*ArcNum
else if(trim(adjustl(SYMTYPE))=="XY")then
Factor=4
TarNum=4
Space=1
end if

ALLOCATE(Xn(1:3,1:4,1:NPanel))
! ALLOCATE(TEMPXn(1:3,1:4,1:NPanel))

do	i=1,	NPanel
	read(40,*)	((xn(j,k,i),j=1,3),k=1,4)				!读入网格节点坐标(外看顺时针排列)
end	do
read(40,*)	SL,SB,temp,temp,ST !垂线间长SL，型宽SB，吃水ST
!---------------------------check file--------------------------------
write(12,'(a)')'#[ 水动力网格信息文件 ]'
write(12,*)
write(12,'(a)')"[浮体网格总数目,浮体湿表面网格数目]"
write(12,*)
write(12,'(2i5)')NPanel,NB
write(12,*)
write(12,'(a)')"[垂线间长(m),型宽(m),吃水(m)]"
write(12,*)
write(12,'(3f8.3)')sl, sb, ST
write(12,*)
!----------------------------------------------------------------------
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$改（Hu）$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------纵倾坐标变换----------
Ty=0.0
!TrimAng=TrimAng/180.0*PI   重复弧度变化 by shiyuyun
Ty(1,1)=cos(TrimAng);	Ty(1,3)=sin(TrimAng);	Ty(2,2)=1.0
Ty(3,1)=-sin(TrimAng);	Ty(3,3)=cos(TrimAng)
!do	i=1,	NPanel
!	do	j=1,	4
!		xn(1:3,j,i)=Matmul(Ty,xn(1:3,j,i))
!	end	do
!end	do
!----------------------------------------------------------------------
CALL GridCheck() !用于略去网格面积小于0.01的单元
ALLOCATE(IEC(1:NPANEL))

CALL PAN1(IEC)!判断面元类型(三角形或四边形)及节点排序规格化处理
ALLOCATE(XQ(1:2,1:4,1:NPanel),xav(1:3,1:NPanel),e(1:3,1:3,1:NPanel),EA(1:NPanel))

CALL PAN2(IEC)!计算面元坐标系基底E、面元中心坐标XAV和面元顶点坐标投影
CALL PAN3()!计算面元投影四边形的面积EA及各顶点相对面元中心之矢径在面元局部坐标系x,y轴的投影XQ
CALL getVol() !程序功能：获得船舶的排水体积和浮心坐标

deltax=XB
do	
	if(abs(XB)>=1e-2)	then
		do	i=1,	Npanel
			xn(1,:,i)=xn(1,:,i)-XB			!x坐标变换为以浮心坐标XB为参照点
		end	do

! 		CALL PAN1(IEC)
        CALL PAN2(IEC)
        CALL PAN3()
        CALL getVol()
		deltax=deltax+XB
	else
		exit
	endif
end	do

!---------------------------check file--------------------------------
write(12,'(a)')"[浮体排水体积(m^3)]"
write(12,*)
write(12,'(a,f12.3)')" 沿X方向积分:",vol_x
write(12,'(a,f12.3)')" 沿Y方向积分:",vol_y
write(12,'(a,f12.3)')" 沿Z方向积分:",vol_z
write(12,'(a,f12.3)')" 平均排水体积: ",vol
write(12,*)
write(12,'(a)')"[浮体浮心坐标(m)]"
write(12,*)
write(12,'(2(a,f8.3))')"XB= ",XB,"    deltax=",deltax
write(12,'(a,f8.3)')"YB= ",YB
write(12,'(a,f8.3)')"ZB= ",ZB
write(12,*)
write(12,'(a)')"#------------------------------------------------------------"
write(12,*)

!----------------------------------------------------------------------
DEALLOCATE(IEC)

END SUBROUTINE READMESH

!------------------------------------------------------------------------------
!用于略去网格面积小于0.01的单元
!------------------------------------------------------------------------------
SUBROUTINE GridCheck()

IMPLICIT NONE

integer(4)	::	Npanel0,NB0
real(8),ALLOCATABLE,DIMENSION(:,:,:)::xn1
real(8)	::	s(1:5),x1(1:3,1:6)
real(8)	::	p1,p2,W,w1,w2
integer(4)	::	ID(1:NPanel)
integer	::	i,j

ALLOCATE(xn1(1:3,1:4,1:NPanel))

xn1(:,:,:)=xn(:,:,:)
Npanel0=Npanel
NB0=NB
ID=0
do i=1,NPanel
	do j=1,4
		x1(:,j)=xn(:,j,i)
	enddo
	x1(:,5)=xn(:,1,i)
	x1(:,6)=xn(:,3,i)
	!计算面元面积
	do j=1,5
		s(j)=sqrt(sum((x1(:,j)-x1(:,j+1))**2))
	enddo
	p1=0.5*(s(1)+s(2)+s(5))
	p2=0.5*(s(3)+s(4)+s(5))
	w1=sqrt(p1*(p1-s(1))*(p1-s(2))*(p1-s(5)))
	w2=sqrt(p2*(p2-s(3))*(p2-s(4))*(p2-s(5)))
	w=w1+w2
	if(w<=1e-10)	ID(i)=1	
enddo
do	i=1,	NPanel
	if(ID(i)==1)	then
		if(i<=NB)	NB0=NB0-1
		NPanel0=Npanel0-1
		do	j=i,	NPanel-1
			xn1(:,:,j)=xn1(:,:,j+1)
		end	do
	end if
end	do

deallocate(xn)
Npanel=Npanel0
NB=NB0
allocate(xn(1:3,1:4,1:Npanel))
xn(1:3,1:4,1:Npanel)=xn1(1:3,1:4,1:Npanel)

deallocate(xn1)

end	subroutine	gridcheck

!--------------------------------------------------------------------------------------
!判断面元类型(三角形或四边形)及节点排序规格化处理
!--------------------------------------------------------------------------------------
subroutine pan1(IEC)
implicit none

integer,intent(out),dimension(1:NPanel)::IEC !面元类型
integer::i,i1,i0,k
real(8),dimension(1:4)::si
do k=1,NPanel
	IEC(k)=4
	i0=0
	do i=1,4
		i1=i+1
		if(i==4)i1=1
		si(i)=sqrt(sum((xn(:,i1,k)-xn(:,i,k))**2))
		if(si(i)<1e-6)then
			IEC(k)=3
			i0=i
			exit
		endif
	enddo
	if(IEC(k)/=4)then
		select case(i0)
		case(1)
			xn(:,2,k)=xn(:,3,k)
			xn(:,3,k)=xn(:,4,k)
		case(2)
			xn(:,3,k)=xn(:,4,k)
		case(4)
			xn(:,4,k)=xn(:,3,k)
		end select
	endif
enddo
end subroutine pan1

!--------------------------------------------------------------------------------------
!计算面元坐标系基底E、面元中心坐标XAV和面元顶点坐标投影
!--------------------------------------------------------------------------------------
subroutine pan2(IEC)
implicit none

integer,intent(in),dimension(1:NPanel)::IEC	 !面元类型
real(8),dimension(1:3)::r13,r24
real(8)::t
integer::i,k

do k=1,NPanel
	xav(:,k)=sum(xn(:,1:IEC(k),k),dim=2)/IEC(k) !计算面元中心
	r13=xn(:,3,k)-xn(:,1,k)
	r24=xn(:,4,k)-xn(:,2,k)
	!计算r13叉乘r24
	e(1,3,k)=r13(2)*r24(3)-r13(3)*r24(2)
	e(2,3,k)=r13(3)*r24(1)-r13(1)*r24(3)
	e(3,3,k)=r13(1)*r24(2)-r13(2)*r24(1)
	e(:,3,k)=e(:,3,k)/(sqrt(sum(e(:,3,k)**2))) !单位化
	do i=1,4
		t=dot_product(e(:,3,k),xn(:,i,k)-xav(:,k))   !形成以e(1:3,3)为法向量且通过中心点xav(1:3)的四边形平面单元
		xn(:,i,k)=xn(:,i,k)-t*e(:,3,k)
	enddo
	e(:,1,k)=xn(:,1,k)-xn(:,3,k)  !形成单元表面向量，适用于含三角形以及四边形单元——2012.07.25  王川
	e(:,1,k)=e(:,1,k)/(sqrt(sum(e(:,1,k)**2))) !单位化
	e(1,2,k)=e(2,3,k)*e(3,1,k)-e(3,3,k)*e(2,1,k)
	e(2,2,k)=e(3,3,k)*e(1,1,k)-e(1,3,k)*e(3,1,k)
	e(3,2,k)=e(1,3,k)*e(2,1,k)-e(2,3,k)*e(1,1,k)
enddo
end subroutine pan2

!--------------------------------------------------------------------------------------
!计算面元投影四边形的面积EA及各顶点相对面元中心之矢径在面元局部坐标系x,y轴的投影XQ
!--------------------------------------------------------------------------------------
subroutine pan3()
implicit none

integer::i,j,k
real(8),dimension(1:2,1:6)::x1
real(8),dimension(1:5)::ss
real(8)::p1,p2,w1,w2,ea1,ea2,eps
eps=1.0e-6
EA=0.0
do k=1,NPanel
    !计算面元顶点局部坐标
	do i=1,4
		do j=1,2
			xq(j,i,k)=dot_product(e(:,j,k),xn(:,i,k)-xav(:,k))
		enddo
		x1(:,i)=xq(:,i,k)
	enddo
	x1(:,5)=xq(:,1,k)
	x1(:,6)=xq(:,3,k)
	!计算面元面积
	do i=1,5
		ss(i)=sqrt(sum((x1(:,i)-x1(:,i+1))**2))
	enddo
	p1=0.5*(ss(1)+ss(2)+ss(5))
	p2=0.5*(ss(3)+ss(4)+ss(5))
	w1=p1*(p1-ss(1))*(p1-ss(2))*(p1-ss(5))
	w2=p2*(p2-ss(3))*(p2-ss(4))*(p2-ss(5))
	ea1=0.0;ea2=0.0
	if(w1>eps)ea1=sqrt(w1)
	if(w2>eps)ea2=sqrt(w2)
	ea(k)=ea1+ea2
enddo
end subroutine pan3

!************************************************************************************
!程序功能：获得船舶的排水体积和浮心坐标
!输入参数：NB为船体网格数目
!          xav(1:3,1:NB)面元中心点坐标, e(1:3,1:3,1:NB)面元基底向量，ea(1:NB)面元面积     
!************************************************************************************
subroutine getVol()
implicit none

integer::i

vol_x=0.0;vol_y=0.0;vol_z=0.0
do i=1,NB
	vol_x=vol_x-Factor*e(1,3,i)*xav(1,i)*ea(i)                                                !计算全船的排水体积，所以在四分之一的模型基础上*4——2012.07.25  王川
	vol_y=vol_y-Factor*e(2,3,i)*xav(2,i)*ea(i)
	vol_z=vol_z-Factor*e(3,3,i)*xav(3,i)*ea(i)
enddo
vol=(vol_x+vol_y+vol_z)/3.0
if(vol<=0.0)then
write(*,*)"水动力模型排水体积为负，请调整水动力网格法向。"
stop
end if
if(abs((vol-vol_x)/vol)>=0.3)then
write(*,*)"水动力模型不闭合，请检查模型。"
stop
end if
XB=0.0;YB=0.0;ZB=0.0
!没有对称面关系
if(Factor==1)then
do i=1,NB
	XB=XB-e(1,3,i)*xav(1,i)*xav(1,i)*ea(i)
	YB=YB-e(2,3,i)*xav(2,i)*xav(2,i)*ea(i)
	ZB=ZB-e(3,3,i)*xav(3,i)*xav(3,i)*ea(i)
enddo
!结构物关于X轴对称
elseif(Factor==2)then
do i=1,NB
	XB=XB-2*e(1,3,i)*xav(1,i)*xav(1,i)*ea(i)
	ZB=ZB-2*e(3,3,i)*xav(3,i)*xav(3,i)*ea(i)
enddo
!结构物关于X轴和Y轴对称
elseif(Factor==4)then
do i=1,NB
	ZB=ZB-4*e(3,3,i)*xav(3,i)*xav(3,i)*ea(i)
enddo
END IF
XB=XB/vol/2.0                                                                                                   
YB=YB/vol/2.0                                                                                                 
ZB=ZB/vol/2.0

end subroutine getVol




END MODULE GETMESH_MOD