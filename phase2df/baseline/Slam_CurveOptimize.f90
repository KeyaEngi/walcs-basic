!输入：num0（粗等分弧长节点进行质量检查并进行圆角过度后的节点数目（大于30））；point0（对粗等分弧长节点进行质量检查并进行圆角过度后的节点坐标）
!m:每条砰击曲线（半横剖线）样条插值点数（理论最多）,lt中给定500;
    
!输出：NumBounPoint（每条曲线拥有的型值点数,500）；Bpoint2（局部坐标系中，每条半横剖线细等分弧长型值点坐标(y,z),用于计算压力）
!Bdx（局部坐标系中，加密剖线斜率，前一节点与后一节点组成线段斜率）；Numincline（每条砰击剖线中需要抛弃积分的区域个数）；inclineB（每个放弃的区域上下限z坐标(局部坐标系下)）
    
!功能：插值获取半横剖线细化节点坐标；获取局部坐标系中，加密剖线节点处斜率；获取每条砰击剖线中需要抛弃积分的区域个数
    !并给出每个放弃的区域上下限z坐标(局部坐标系下)

subroutine Slam_CurveOptimize(num0,point0,m,NumBounPoint,&
                              & Bpoint2,Bdx,Numincline,inclineB    )

   implicit none

   integer(4)::num0
   real(8),dimension(num0,2)::point0  !---原始插值节点

   integer(4)::Num3
   real(8),dimension(num0,2)::point3  !---对初始直线段的简易处理
   
   integer(4)::m    !----总插值点数
   integer(4)::NumBounLine   !---可用曲线段数目
   integer(4)::NumBounPoint   !---每条可用曲线对应插值点数
   real(8),dimension(m,2)::Bpoint2   !----插值节点坐标(平移曲线)
   real(8),dimension(m)::Bdx,Bddx   !---插值节点对应一阶导数，二阶导数
   !-------对于内凹区域在实际积分时，不进行积分,需要根据高度确定内凹区间
   integer(4)::Numincline     !------处理的内凹区间个数
   real(8),dimension( 100,2 )::inclineB   !---对应内凹区间的上下限(Z,局部坐标系下)



   integer(4),dimension(4,2)::BounLine   !----型线区间
   integer(4)::m1
   real(8),allocatable,dimension(:,:)::point1  !---加密插值节点
   real(8),allocatable,dimension(:)::dx,ddx   !---插值节点对应一阶导数，二阶导数
   real(8),allocatable,dimension(:)::length
   real(8),allocatable,dimension(:,:)::func_xt,func_yt
   real(8),allocatable,dimension(:,:)::tempP

   real(8)::ds1,dsn    !---边界条件，以一阶导数作为边界条件
   real(8),allocatable,dimension(:,:)::IntPoint
   integer(4)::mLine

   !***********************处理内部凹陷区域
   real(8)::temAngle,temAngle2
   integer(4)::Bdry1,Bdry2
   real(8)::Nsbdry1,Nsbdry2      !-----边界斜率
   real(8),dimension(3,2)::temBdryP   !---边界点

   real(8),dimension(2,2)::circleP
   real(8)::lk1,lk2,dlk1,dlk2
   real(8)::lb1,lb2,dlb1,dlb2
   real(8)::leng1,leng2

   real(8),dimension(m,2)::temNode
   real(8),dimension(m)::temdx,temddx

   !-------等弧长均分曲线
   integer(4)::m2
   real(8),dimension(m,2)::temNode2
   real(8),dimension(m)::temdx2,temddx2


   real(8)::Pi

   real(8)::a1,a2,a3,a4,b1,b2,b3,b4
   real(8)::xt,xtt,yt,ytt

   integer(4)::i,j,k,k1,k2,i1,i2,ii,jj,kk,iii,jjj,kkk
   real(8)::s,t,s1,s2,t1,t2,t3,s3
   real(8),dimension(3)::temx1,temx2
   real(8)::Kexi

!*********************
  Pi=atan(1.0)*4.0

!!  open(unit=21,file='check1.txt')
!!  open(unit=22,file='check2.txt')

 mLine=10

 !---处理初始直线段
 point3=0.0;

 s=(point0(2,2)-point0(1,2) )/(point0(2,1)-point0(1,1) )     !前两个节点构成直线的斜率
 k=0;
 do i=3,num0
   s1=(point0(i,2)-point0(1,2) )/(point0(i,1)-point0(1,1) )   !第i（i>3）节点与第1个节点构成直线的斜率

   if( abs(s1-s)<=0.000001 ) then
     k=k+1;
   else
     exit
   end if
 end do
 k=2+k    !---曲线首段出现直线的节点数目(前两个节点组成的线段默认为直线)

 Num3=1;
 point3(Num3,:)=point0(1,:)    !将剔除直线节点之外粗等分弧长节点赋给point3，但保留了第一节点   
 do i=k,num0
   Num3=Num3+1;
   point3(Num3,:)=point0(i,:)     
 end do

 !----尝试
 num0=Num3;    !剔除直线节点数目之后的粗选节点数目
 point0(1:num0,1:2)=point3(1:Num3,1:2)    !再将point3赋给point0，现在point0已经剔除直线节点

!! do i=1,num0
!!   write(22,"(i8,2(f20.9,1x))") i,point0(i,1),point0(i,2)
!! end do



 !----强制后面的点比前一个点高
   s=point0(1,2)    !第一个节点z坐标（局部坐标系中，坐标只含有y，z，2则代表z坐标）
   do i=2,num0
     if( point0(i,2)<=s ) then   !如果后一个节点z坐标比前一个节点z坐标小，则强制后一个节点z坐标等于前节点z坐标加0.01
       point0(i,2)=s+0.01
     end if
     s=point0(i,2)
   end do


!***************************此处开始执行曲线优化操作




  if( num0==2 ) then   !----整个剖面轮廓为直线段
      !----先限制底升角
      s=85.0/180.0*Pi      !低升角最大85度
      !----斜率
      t=(point0(2,2)-point0(1,2) )/(point0(2,1)-point0(1,1) )
      if( t>tan(s) ) then         !如果低升角大于85度，则默认为85度
        point0(2,2)=point0(1,2)+tan(s)*(  point0(2,1)-point0(1,1) )
      end if

      t1=(point0(2,2)-point0(1,2) )/(point0(2,1)-point0(1,1) )

      NumBounLine=1;       !可用曲线段数目
      NumBounPoint=m       !每条可用曲线对应插值点数 m=500 
      do i=1,m
        Kexi=real(i-1)/real(m-1)

        Bpoint2(i,:)=point0(1,:)*(1.0-Kexi)+point0(2,:)*Kexi    !细化后的节点坐标（线性插值）
        Bdx(i)=t1        !局部坐标系中，加密剖线斜率
        Bddx(i)=0.0
      end do

      Numincline=0;      !单条砰击剖线中需要抛弃积分的区域个数
      inclineB=0.0;      !每个放弃的区域上下限(局部坐标系下)

  else          !如果整个剖线轮廓不是直线段
  !------加入条件，外侧内倾或直壁起始点位
  s=point0(num0,1)      !最高点的y值
  kk=0
  do i=num0-1,1,-1
      t=point0(i,1)

      if( (t-s)>=0.0 ) then      !由高到底检查砰击点y值与最高点y值关系，判断船体舷板上部直壁或内倾的位置
        kk=kk+1
      else
        exit
      end if
      s=t;
  end do
  kk=num0-kk;   !----节点上限,后面的点位直壁或内倾壁（直壁没问题，内倾是不是有一点问题）

!*************************************
!---------开始插值
  !----第一段采用线性插值，为后面的曲线插值提供边界条件

  allocate( IntPoint(num0-1,2) )  !---构造曲线插值节点
  do i=1,num0-1
    IntPoint(i,:)=point0(i+1,:)    !将第二到第num0个节点的坐标赋给IntPoint
  end do

  !----第一段线段的斜率
  t2=(point0(2,2)-point0(1,2))/(point0(2,1)-point0(1,1))

  !----开始累加弦长插值加密
  allocate( length(num0-1) )
  length=0.0;
  s=0.0
  do i=2,num0-1
!!    s=s+sqrt( (point0(i,1)-point0(i-1,1))**2.0+(point0(i,2)-point0(i-1,2))**2.0 )
    s=s+sqrt( (IntPoint(i,1)-IntPoint(i-1,1))**2.0+(IntPoint(i,2)-IntPoint(i-1,2))**2.0 )    
    length(i)=s     !从第一个节点到num0-1个节点的累加弦长
  end do

  allocate(func_xt(num0-2,4),func_yt(num0-2,4)  )
  allocate( tempP(num0-1,2) )
  func_xt=0.0;  func_yt=0.0

!-----求解插值参数函数func_xt;func_yt  
  if( kk==num0 ) then   !----不存在直壁

    tempP(:,1)=length(:);
    tempP(:,2)=IntPoint(:,1);
    !----给定边界条件
    ds1=(IntPoint(2,1)-IntPoint(1,1) )/(length(2)-length(1))
    dsn=(IntPoint(num0-1,1)-IntPoint(num0-2,1) )/(length(num0-1)-length(num0-2))
    !-----插值参数函数xt
    call spline3(num0-1,tempP,ds1,dsn,func_xt )

    !----给定边界条件
    ds1=t2*ds1
    dsn=(IntPoint(num0-1,2)-IntPoint(num0-2,2) )/(length(num0-1)-length(num0-2))*1.5
    tempP(:,2)=IntPoint(:,2);
    !-----插值参数函数yt
    call spline3(num0-1,tempP,ds1,dsn,func_yt )

  else     !----存在直壁
    
    tempP(1:kk-1,1)=length(1:kk-1);
    tempP(1:kk-1,2)=IntPoint(1:kk-1,1);
    !----给定边界条件
    ds1=(IntPoint(2,1)-IntPoint(1,1) )/(length(2)-length(1))
    dsn=(IntPoint(kk-1,1)-IntPoint(kk-2,1) )/(length(kk-1)-length(kk-2))
    !-----插值参数函数xt
    call spline3(kk-1,tempP(1:kk-1,:),ds1,dsn,func_xt(1:kk-2,:) )

    !**********************插yt
    !----给定边界条件
    t3=(IntPoint(kk-1,2)-IntPoint(kk-2,2) )/(length(kk-1)-length(kk-2))

    ds1=t2*ds1
    dsn=t3*2.0
    tempP(1:kk-1,2)=IntPoint(1:kk-1,2);
    !-----插值参数函数yt
    call spline3(kk-1,tempP(1:kk-1,:),ds1,dsn,func_yt(1:kk-2,:) )

  end if


!*********************开始插值
  !-----节点加密,总共m个
  allocate( point1(m,2) )
  allocate( dx(m),ddx(m)  )

  s=length(kk-1)/real(m-mLine-1)
  ii=0;
  do i=mLine+1,m
    
    t=s*real(ii)
    if(i==m) t=length(kk-1)

    do j=1,num0-2   !---找区间
      s1=length(j);  s2=length(j+1);
      if( (t-s1)*(t-s2)<=0.0 ) then
          k=j
          exit
      end if
    end do

    a1=func_xt(k,1);  a2=func_xt(k,2);
    a3=func_xt(k,3);  a4=func_xt(k,4);

    b1=func_yt(k,1);  b2=func_yt(k,2);
    b3=func_yt(k,3);  b4=func_yt(k,4);


    point1(i,1)=func_xt(k,1)*t**3.0+func_xt(k,2)*t**2.0+func_xt(k,3)*t+func_xt(k,4)
    point1(i,2)=func_yt(k,1)*t**3.0+func_yt(k,2)*t**2.0+func_yt(k,3)*t+func_yt(k,4)

    xt=3.0*a1*t**2.0+2.0*a2*t+a3
    xtt=6.0*a1*t+2.0*a2
    yt=3.0*b1*t**2.0+2.0*b2*t+b3
    ytt=6.0*b1*t+2.0*b2

    dx(i)=yt/xt;
    ddx(i)=(ytt*xt-xtt*yt )/xt**3.0

    ii=ii+1

  end do

  s1=point0(1,1);  s2=point0(2,1);    !----直线段插值
  t1=point0(1,2);  t2=point0(2,2);
  do i=1,mLine
    t=real(i-1)/real(mLine+1)

    point1(i,1)=s1*(1.0-t)+s2*t;
    point1(i,2)=t1*(1.0-t)+t2*t;

    dx(i)=(t2-t1)/(s2-s1);
    ddx(i)=0.0
  end do

  deallocate( length,func_xt,func_yt,tempP )


    !-----处理舷外侧直壁或内倾壁
    t3=tan( 85.0/180.0*Pi )    !----斜率上限

    k=0;
    s=point1(m,1);
    t=dx(m)

    if( t<0.0.or.t>t3 ) then
      k=k+1
    end if

    do i=m-1,1,-1
      s1=point1(i,1);   !----点坐标
      t1=dx(i);         !----节点处斜率
      if(s1>=s.or.t1<0.0.or.t1>t3  ) then
        k=k+1;
      else
        exit
      end if
      s=s1;
    end do
    m1=m-k;    !-----确定抛开外侧直壁后的插值点数目

    !----先对各点高度关系做出限制
    !----处理平底部分的插值
    i=1
    do 
      
      if(i>=m1) exit

      t1=point1(i,2)
      t2=point1(i+1,2)

      if(t2<=t1) then   !---确定低于前一点的区间
        k=1
        do j=i+2,m1
          t2=point1(j,2)
          if( t2<=t1 ) then
            k=k+1
          else
            exit
          end if
        end do

        if( (i+k+1)>m1-1  ) then
          t2=t1+0.01
          do j=i+1,i+k
            s=(point1(j,1)-point1(i,1) )/(point1(m1,1)-point1(i,1) )
            point1(j,2)=t1*(1.0-s)+t2*s

            dx(j)=(point1(j,2)-point1(j-1,2) )/(point1(j,1)-point1(j-1,1) )
            ddx(j)=(dx(j)-dx(j-1) )/(point1(j,1)-point1(j-1,1) )
          end do
        else
          t2=point1(i+k+1,2)

          do j=i+1,i+k
            s=(point1(j,1)-point1(i,1) )/(point1(i+k+1,1)-point1(i,1) )
            point1(j,2)=t1*(1.0-s)+t2*s

            dx(j)=(point1(j,2)-point1(j-1,2) )/(point1(j,1)-point1(j-1,1) )
            ddx(j)=(dx(j)-dx(j-1) )/(point1(j,1)-point1(j-1,1) )
          end do

          point1(i+k+1,2)=point1(i+k,2)*0.7+point1(i+k+2,2)*0.3
          dx(i+k+1)=(point1(i+k+1,2)-point1(i+k,2) )/(point1(i+k+1,1)-point1(i+k,1) )

          dx(i+k+2)=(point1(i+k+2,2)-point1(i+k+1,2) )/(point1(i+k+2,1)-point1(i+k+1,1) )
        end if

        i=i+k
      end if

      i=i+1
    end do


!!  do i=1,m
!!    write(21,"(i8,4(f20.9,1x))") i,point1(i,1),point1(i,2),dx(i),ddx(i)
!!  end do
!!  stop

    !*******************************     
    !-------处理内部凹陷区域
    Numincline=0;
    inclineB=0.0;

    temAngle=85.0/180.0*Pi
    temAngle2=80.0/180.0*Pi

    temNode=0.0;  temdx=0.0;  temddx=0.0;

    kk=0;
    i=1;
    do 
       
       !---当前节点编号
       kk=kk+1
       temNode(kk,:)=point1(i,:);
       temdx(kk)=dx(i);
       temddx(kk)=ddx(i);        

       if( i>=m1-1 ) then
           kk=kk+1
           temNode(kk,:)=point1(i+1,:);
           temdx(kk)=dx(i+1);
           temddx(kk)=ddx(i+1); 
           exit           
       end if

       !----起始节点
       s=point1(i,1);
       t=point1(i,2);
       
       jj=0;
       do j=i+1,m1   !---考察节点i后面的一段是否存在凹/直壁曲面
         s1=point1(j,1);
         t1=point1(j,2);
         s2=dx(j)

         b1=t-tan(TemAngle2)*s    !---斜线在y轴的交点
         b2=(t1-b1 )/tan(TemAngle2)   !---按临界斜线算，等高下，节点相对位置

         if(s1<=b2.or.s2<0.0.or.s2>tan(TemAngle2 )) then
           jj=jj+1;
         else
           exit
         end if

       end do

       if(jj==0) then
         i=i+1
       else

          Numincline=Numincline+1;
          inclineB(Numincline,1)=t;   !----凹陷区域起始位置

         !-----先排除外侧节点
         if(i+jj==m1 )  then

           !----采用斜线代替
           s1=dx(i)*1.5;
           t1=t-s1*s;   !----与纵轴交点

           t2=point1(m1,2);   !----重新计算最后一个点
           s2=(t2-t1)/s1;

           do j=i+1,i+jj
               kk=kk+1;
               Kexi=real(j-i )/real(jj)

               temNode(kk,1)=s*(1.0-Kexi)+s2*Kexi;
               temNode(kk,2)=t*(1.0-Kexi)+t2*Kexi;
               temdx(kk)=s1;
               temddx(kk)=0.0;
           end do

           inclineB(Numincline,2)=t2   !---凹陷结束位置

           exit
         end if         

           jjj=i+jj+1;
           !------------构造直线连接圆弧的连接方式
           if(i==1) then
               if( dx(i)>tan(TemAngle) ) then
                   lk1=tan(TemAngle);
                   dx(i)=tan(TemAngle);
                   temdx(kk)=dx(i);
               elseif( dx(i)<tan(TemAngle2) ) then
                   dx(i)=tan( (TemAngle+TemAngle2)/2.0 )
                   temdx(kk)=dx(i);
                   lk1=dx(i);   !----第一条直线斜率
               else
                   lk1=dx(i);
               end if
           else
               if( dx(i)>tan(TemAngle) ) then
                   lk1=tan(TemAngle);
               elseif( dx(i)<tan(TemAngle2) ) then
                   lk1=tan( (TemAngle+TemAngle2)/2.0 )
               else
                   lk1=dx(i);
               end if
           end if
           
           dlk1=tan( atan(lk1)+Pi/2.0 )  !----垂线斜率
           lb1=t-lk1*s   !---line1与Y轴交点

           !---------------计算，确定圆弧信息
           circleP=0.0;
           circleP(2,:)=point1(jjj,:)   !---圆弧艉部点
           lk2=dx(jjj)  !---第二条直线斜率
           dlk2=tan( atan(lk2)+Pi/2.0 )   !----垂线斜率

           lb2=circleP(2,2)-lk2*circleP(2,1)  !---line2与Y轴交点
           dlb2=circleP(2,2)-dlk2*circleP(2,1)  !---line2垂线y轴交点
           !----计算line1与line2的交点
           s1=(lb1-lb2 )/(lk2-lk1)
           t1=lk1*s1+lb1
           !-----线段长度
           a1=sqrt( (s1-circleP(2,1) )**2.0+(t1-circleP(2,2) )**2.0 )
           a2=sqrt( (s1-s)**2.0+(t1-t)**2.0 )
           !-----计算比例系数
           Kexi=(a2-a1)/a2
           circleP(1,1)=s*(1.0-Kexi)+s1*Kexi;  !---圆弧首端点
           circleP(1,2)=t*(1.0-Kexi)+t1*Kexi;

           dlb1=circleP(1,2)-dlk1*circleP(1,1)  !---line1垂线y轴交点
           !-----计算圆心
           s2=(dlb1-dlb2 )/(dlk2-dlk1)
           t2=dlk1*s2+dlb1
           !-----计算半径
           b1=sqrt( (s2-circleP(1,1) )**2.0+(t2-circleP(1,2) )**2.0 )
           !-----计算圆弧端点对应弧度(第二象限)
           b2=atan(lk1)+Pi/2.0
           b3=atan(lk2)+Pi/2.0

           !------计算直线段与圆弧段的长度
           leng1=sqrt( (s-circleP(1,1) )**2.0+(t-circleP(1,2) )**2.0 )
           leng2=b1*abs(b3-b2)

           Kexi=leng1/(leng1+leng2)
           !-----开始插值
           iii=int(real(jj)*Kexi)   !----直线部分
           do j=1,iii
             Kexi=real(j)/real(iii)

             kk=kk+1
             temNode(kk,1)=s*(1.0-Kexi)+circleP(1,1)*Kexi
             temNode(kk,2)=t*(1.0-Kexi)+circleP(1,2)*Kexi
             temdx(kk)=lk1;
             temddx(kk)=0.0; 
           end do

           iii=jj-iii   !---圆弧部分
           do j=1,iii
             Kexi=real(j)/real(iii+1)

             b4=b2*(1.0-Kexi)+b3*Kexi

             kk=kk+1
             temNode(kk,1)=b1*cos(b4)+s2   !----点坐标
             temNode(kk,2)=b1*sin(b4)+t2
             
             temdx(kk)=tan( b4-Pi/2.0 )
             temddx(kk)=0.0;

           end do

           inclineB(Numincline,2)=temNode(kk,2)   !---凹陷结束位置

         i=i+jj+1
       end if

    end do

    m1=kk

    point1(1:m1,:)=temNode(1:m1,:);
    dx(1:m1)=temdx(1:m1);
    ddx(1:m1)=temddx(1:m1);

    !----先对各点高度关系做出限制
    !----处理平底部分的插值
    i=1
    do 
      
      if(i>=m1) exit

      t1=point1(i,2)
      t2=point1(i+1,2)

      if(t2<=t1) then   !---确定低于前一点的区间
        k=1
        do j=i+2,m1
          t2=point1(j,2)
          if( t2<=t1 ) then
            k=k+1
          else
            exit
          end if
        end do

        if( (i+k+1)>m1-1  ) then
          t2=t1+0.01
          do j=i+1,i+k
            s=(point1(j,1)-point1(i,1) )/(point1(m1,1)-point1(i,1) )
            point1(j,2)=t1*(1.0-s)+t2*s

            dx(j)=(point1(j,2)-point1(j-1,2) )/(point1(j,1)-point1(j-1,1) )
            ddx(j)=(dx(j)-dx(j-1) )/(point1(j,1)-point1(j-1,1) )
          end do
        else
          t2=point1(i+k+1,2)

          do j=i+1,i+k
            s=(point1(j,1)-point1(i,1) )/(point1(i+k+1,1)-point1(i,1) )
            point1(j,2)=t1*(1.0-s)+t2*s

            dx(j)=(point1(j,2)-point1(j-1,2) )/(point1(j,1)-point1(j-1,1) )
            ddx(j)=(dx(j)-dx(j-1) )/(point1(j,1)-point1(j-1,1) )
          end do

          point1(i+k+1,2)=point1(i+k,2)*0.7+point1(i+k+2,2)*0.3
          dx(i+k+1)=(point1(i+k+1,2)-point1(i+k,2) )/(point1(i+k+1,1)-point1(i+k,1) )

          dx(i+k+2)=(point1(i+k+2,2)-point1(i+k+1,2) )/(point1(i+k+2,1)-point1(i+k+1,1) )
        end if

        i=i+k
      end if

      i=i+1
    end do


    do i=2,m1
        if(dx(i)<=0.0) then
            dx(i)=dx(i-1)
            point1(i,1)=point1(i-1,1)+(point1(i,2)-point1(i-1,2) )/dx(i)
            ddx(i)=0.0;
        end if
    end do



    temNode=0.0;  temdx=0.0;  temddx=0.0;

    temNode(1:m1,:)=point1(1:m1,:);
    temdx(1:m1)=dx(1:m1);
    temddx(1:m1)=ddx(1:m1);


!!  do i=1,m1
!!    write(21,"(i8,4(f20.9,1x))") i,point1(i,1),point1(i,2),dx(i),ddx(i)
!!  end do
!!
!!  stop


!!  do i=1,m1
!!    write(21,"(i8,4(f20.9,1x))") i,point1(i,1),point1(i,2),dx(i),ddx(i)
!!  end do
!!  stop

      
      !-----20210406修改
      !-----增加最终判断机制
      i=1;
      do 

          if(i>=m1) exit
          s1=temNode(i,1);
          t1=temNode(i,2);

          s2=temNode(i+1,1);
          t2=temNode(i+1,2);

          if( s1>=s2.or.t1>=t2 ) then
              do j=i+1,m1-1
                  temNode(j,:)=temNode(j+1,:);
                  temdx(j)=temdx(j+1);
                  temddx(j)=temddx(j+1);
              end do
              if(i+1>m1-1) then
                  m1=m1-1
                  exit
              end if
              m1=m1-1;
              
              cycle
          end if
          i=i+1;
      end do

      do i=2,m1
          s1=temdx(i)
          if( s1<=0.0 ) then
              temdx(i)=temdx(i-1)
          end if
      end do

!!!*****************凹陷区域处理完毕
!!       NumBounLine=1;
!!       NumBounPoint=m1;
!!
!!       Bpoint2(1:m1,:)=temNode(1:m1,:);
!!       Bdx(1:m1)=temdx(1:m1);
!!       Bddx(1:m1)=temddx(1:m1);     


!****************************************
      !-----20210419
      !---------采用等弧长均布节点及其导数

      allocate( length( m1 ) )
      length=0.0;

      s=0.0;
      do i=2,m1
          temx1(1:2)=temNode(i-1,1:2);
          temx2(1:2)=temNode(i,1:2);

          t=sqrt( ( temx1(1)-temx2(1) )**2.0+( temx1(2)-temx2(2) )**2.0 )
          s=s+t;

          length(i)=s;
      end do

      !-----开始等弧长划分节点
      m2=m;
      temNode2=0.0;
      temdx2=0.0; temddx2=0.0;

      temNode2(1,:)=temNode(1,:);   !---起始点
      temdx2(1)=temdx(1);
      temddx2(1)=0.0;

      xt=length(m1)/real( m2-1 )
      do i=2,m2-1
          xtt=xt*real(i-1)    

          do ii=1,m1-1
              s1=length(ii);
              s2=length(ii+1);

              if( ( s1-xtt )*( s2-xtt )<=0.0 ) then
                  if( abs( s2-s1 )<=1.0e-8 ) then
                      Kexi=0.5
                  else
                      Kexi=( xtt-s1 )/(s2-s1)
                  end if

                  temNode2(i,:)=temNode(ii,:)*(1.0-Kexi)+temNode(ii+1,:)*Kexi
                  temdx2(i)=temdx(ii)*(1.0-Kexi)+temdx(ii+1)*Kexi
                  temddx2(i)=0.0
                  exit
              end if
          end do


      end do

      !----终止点
      temNode2(m2,:)=temNode(m1,:);   !---起始点
      temdx2(m2)=temdx(m1);
      temddx2(m2)=0.0;

!!      !-----判断各点能否积分
!!
!!      if( Numincline==0 ) then   !---整段都能积分
!!          temSlamItgcase(1:m2)=1;
!!      else
!!          temSlamItgcase(1:m2)=1;    !-----存在某段不能积分
!!          do j=1,m2
!!              xt=temNode2(j,2)
!!
!!              do i=1,Numincline
!!                  s1=inclineB(i,1);
!!                  s2=inclineB(i,2);
!!
!!                  if( ( xt-s1 )*( xt-s2 )<=0.0 ) then
!!                      temSlamItgcase(j)=0
!!                      exit
!!                  end if
!!              end do
!!          end do
!!      end if



!*****************凹陷区域处理完毕
       NumBounLine=1;       !可用曲线段数目
       NumBounPoint=m2;     !每条可用曲线对应插值点数m2=m       

       Bpoint2(1:m2,:)=temNode2(1:m2,:);     !局部坐标系中，每条曲线型值点坐标(y,z),细化后
       Bdx(1:m2)=temdx2(1:m2);     !局部坐标系中，加密剖线斜率（插值节点对应的一阶导数）
       Bddx(1:m2)=temddx2(1:m2);       !插值节点对应的二阶导数  
 
    end if



      return
end subroutine Slam_CurveOptimize

