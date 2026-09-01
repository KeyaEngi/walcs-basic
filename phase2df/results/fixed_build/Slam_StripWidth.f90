!输入：SlamNumLine（典型砰击节点数）；SlamAngle（考虑船体纵倾影响的砰击角）；SlamNumLP（半横剖线上的细化节点数）
!SlamNode(用户坐标系下，砰击节点坐标) ；SlamIntNumP（砰击剖线典型节点）；SlamStripTyp（切片厚度典型节点）   
!SlamWidthMcoef（典型砰击节点位置处对应的带宽修正系数）    
!输出：SlamStripBexist（切片边界划分与否状态0未划分，1划分）指的是切片厚度方向
!SlamStripBAngle（边界倾角，随船平动坐标系，与SlamAngle类似）指的是切片厚度方向    
!temStripWidth（剖线节点对应带宽,为后续积分点提供插值）；SlamStripBNode（切片左右两侧边界节点坐标）    

!功能：计算剖线左右两侧边界倾角（指的是切片厚度方向）；剖线各层节点（竖直方向）对应带宽（切片厚度方向）；切片左右两侧边界节点坐标（切片厚度方向）  
    
    
    
!--------计算剖线切片带宽
subroutine Slam_StripWidth(Numline,temStripWidth,temStripBNode )

  use Slamming,only:SlamNumLine,SlamAngle,SlamStripBexist,SlamStripBAngle,SlamNumLP,&
                  & SlamNode,SlamIntNumP,SlamStripTyp,SlamWidthMcoef
  use Constant,only:Pi
  implicit none

  integer(4)::Numline
  real(8),dimension(SlamIntNumP  )::temStripWidth
  real(8),dimension(SlamIntNumP,2,3 )::temStripBNode

  integer(4)::Num0
  real(8),dimension(SlamIntNumP,3)::temSlamNode


  real(8)::temAngle                 !---砰击计算剖线倾角
  real(8)::temAngle1,temAngle2      !---切片左右边界倾角

  real(8),dimension(3)::x0,x1,x2
  real(8)::k0,b0,k1,b1,k2,b2

  real(8),dimension(3)::xx1,xx2    !----交点

  real(8)::s,t,s1,s2
  integer(4)::i,j,k


  !**********先判断剖线所在左右边界是否已划分
  !----砰击计算剖线倾角
  temAngle=SlamAngle(Numline )

  if( Numline==1 ) then
        !-----左边界
        if(SlamStripBexist(Numline,1)==0  ) then
            SlamStripBexist(Numline,1)=1;
            SlamStripBAngle(Numline,1)=temAngle;
        end if

        !-----右边界
        if(SlamStripBexist(Numline,2)==0  ) then
            !----先判断相邻曲线在该边界上有无划分
            if( SlamStripBexist(Numline+1,1)==1 ) then
                SlamStripBexist(Numline,2)=1;
                SlamStripBAngle(Numline,2)=SlamStripBAngle(Numline+1,1)
            elseif( SlamStripBexist(Numline+1,1)==0 ) then
                SlamStripBexist(Numline,2)=1;
                SlamStripBAngle(Numline,2)=temAngle;
            end if
        end if

  elseif(Numline==SlamNumLine ) then
    
        !-----左边界
        if(SlamStripBexist(Numline,1)==0  ) then
            !----先判断相邻曲线在该边界上有无划分
            if( SlamStripBexist(Numline-1,2)==1 ) then
                SlamStripBexist(Numline,1)=1;
                SlamStripBAngle(Numline,1)=SlamStripBAngle(Numline-1,2)

            elseif( SlamStripBexist(Numline-1,2)==0 ) then
                SlamStripBexist(Numline,1)=1;
                SlamStripBAngle(Numline,1)=temAngle;

            end if
        end if

        !-----右边界
        if(SlamStripBexist(Numline,2)==0  ) then
            SlamStripBexist(Numline,2)=1;
            SlamStripBAngle(Numline,2)=temAngle;
        end if


  else
    
        !-----左边界
        if(SlamStripBexist(Numline,1)==0  ) then
            !----先判断相邻曲线在该边界上有无划分
            if( SlamStripBexist(Numline-1,2)==1 ) then
                SlamStripBexist(Numline,1)=1;
                SlamStripBAngle(Numline,1)=SlamStripBAngle(Numline-1,2)

            elseif( SlamStripBexist(Numline-1,2)==0 ) then
                SlamStripBexist(Numline,1)=1;
                SlamStripBAngle(Numline,1)=temAngle;

            end if
        end if
    
        !-----右边界
        if(SlamStripBexist(Numline,2)==0  ) then
            !----先判断相邻曲线在该边界上有无划分
            if( SlamStripBexist(Numline+1,1)==1 ) then
                SlamStripBexist(Numline,2)=1;
                SlamStripBAngle(Numline,2)=SlamStripBAngle(Numline+1,1)
            elseif( SlamStripBexist(Numline+1,1)==0 ) then
                SlamStripBexist(Numline,2)=1;
                SlamStripBAngle(Numline,2)=temAngle;
            end if
        end if

  end if
  !---------------经上述判断，可以确定本曲线的左右边界的倾角
  
  temAngle1=SlamStripBAngle(Numline,1)     !第Numline条剖线左边界倾角
  temAngle2=SlamStripBAngle(Numline,2)     !第Numline条剖线右边界倾角
  
  !******************************************
  !------开始计算剖线带宽
  !-----计算压力的剖线节点在随船平动坐标系下的三维坐标
  temSlamNode=0.0;
  Num0=SlamNumLP(Numline )
  temSlamNode(1:Num0,1:3)=SlamNode(Numline,1:Num0,1:3 )

  !-----左右界限典型点
  x1(:)=SlamStripTyp(Numline,1,: )
  x2(:)=SlamStripTyp(Numline,2,: )

  if( abs(temAngle1-Pi/2.0 )>1.0e-5 ) then  !---计算左右界限直线参数
    k1=tan(temAngle1 )
    b1=x1(3)-k1*x1(1)
  end if

  if( abs(temAngle2-Pi/2.0 )>1.0e-5 ) then
    k2=tan(temAngle2 )
    b2=x2(3)-k2*x2(1)
  end if

  !***********开始计算带宽
  temStripWidth=0.0;   temStripBNode=0.0;

  do i=1,Num0    !对剖先节点数进行循环

     x0(:)=temSlamNode(i,1:3);

     k0=tan( temAngle-Pi/2.0  )
     b0=x0(3)-k0*x0(1)

     !-----左舷交点
     if( abs(temAngle1-Pi/2.0 )<=1.0e-5 ) then
       xx1(1)=x1(1);
       xx1(2)=x0(2);
       xx1(3)=k0*xx1(1)+b0;
     else
       xx1(1)=(b0-b1)/(k1-k0)
       xx1(2)=x0(2);
       xx1(3)=k0*xx1(1)+b0;
     end if

     !-----右舷交点
     if( abs(temAngle2-Pi/2.0 )<=1.0e-5 ) then
       xx2(1)=x2(1);
       xx2(2)=x0(2);
       xx2(3)=k0*xx2(1)+b0;
     else
       xx2(1)=(b0-b2)/(k2-k0)
       xx2(2)=x0(2);
       xx2(3)=k0*xx2(1)+b0;
     end if

     !------判断计算的交点是否合理
     if( xx1(1)>x0(1) ) then
       xx1(:)=x0(:)
     end if
     if( xx2(1)<x0(1) ) then
       xx2(:)=x0(:)
     end if

     !----计算带宽
     s=sqrt( (xx1(1)-xx2(1))**2.0+(xx1(3)-xx2(3))**2.0 )
     temStripWidth(i)=s

     !-----计算左右侧边界上的节点坐标

     temStripBNode(i,1,1:3)=xx1(1:3);

     temStripBNode(i,2,1:3)=xx2(1:3);

!     s1=sqrt( (xx1(1)-x0(1) )**2.0+(xx1(3)-x0(3) )**2.0 )   !---左侧带宽
!     s2=sqrt( (xx2(1)-x0(1) )**2.0+(xx2(3)-x0(3) )**2.0 )   !---右侧带宽
!
!     if( s<=1.0e-6 ) then
!         temStripWidthRatio(i,1)=0.5;
!         temStripWidthRatio(i,2)=0.5;
!     else
!         temStripWidthRatio(i,1)=s1/s;
!         temStripWidthRatio(i,2)=1.0-temStripWidthRatio(i,1);
!     end if

  end do

!!  !************带宽修改
!!  temStripWidth=temStripWidth*0.6

  temStripWidth=temStripWidth*SlamWidthMcoef( Numline )




  return
end subroutine Slam_StripWidth