!输入：Num0每条半横剖线细化后拥有的型值点数（最多500）；Point0：局部坐标系下每条半横剖线细化后拥有的型值点坐标
!Sdx局部坐标系下加密剖线的斜率；NumIntC剖线半宽等分Ct的点数100    
!输出：IntC，选定Ct插值点，不包括Ct=0；IntNC，Ct对应的侵入深度ht,不包括Ct=0
!DerNc，D(ht)/D(Ct)，包括Ct=0；temRise，Ct对应接触点位置处液面抬升高度,不包括Ct=0    
    
!功能：获得Ct对应的插值点IntC（y方向），进而计算插值点IntC位置对应的浸湿深度IntNC；
    !另外，还计算了z向浸湿深度在y方向Ct插值点的变化率D(ht)/D(Ct)，最后给出IntC位置对应的液面抬升高度（0自由液面以上）
    
    
    
!----根据Wagner条件，计算剖线浸入深度ht与接触点Ct的关系

subroutine Slam_WagCondition(Num0,Point0,Sdx,NumIntC,IntC,IntNC,DerNc,temRise)

  use Constant,only:Pi
  implicit none

  integer(4)::Num0    !----剖线节点数目
  real(8),dimension(1:Num0,1:2)::Point0   !----剖线节点坐标
  real(8),dimension(1:Num0)::Sdx     !----剖线节点对应一阶导数

  integer(4)::NumIntC
  real(8),dimension(NumIntC)::IntC
  real(8),dimension(NumIntC)::IntNC
  real(8),dimension(NumIntC+1)::DerNc

  real(8),dimension(NumIntC )::temRise

  real(8),dimension(1:Num0,1:2)::Point1   !---变换坐标系的剖线节点
  real(8)::width
  real(8)::temC

  integer(4)::RomK    !---Romberg积分阶次
!!  real(8),allocatable,dimension(:,:)::RomTkj


  real(8),allocatable,dimension(:)::length
  real(8),allocatable,dimension(:,:)::temp1,temp2
  real(8),allocatable,dimension(:,:)::funcx,funcy
  real(8)::ds1,dsn

  integer(4)::TemI,TemI1
  integer(4)::i,j,k,ii,jj,kk
  real(8)::s,t,s1,s2,t1,t2
  real(8)::a1,a2,a3,a4,b1,b2,b3,b4
  real(8)::a,b          !---Romberg积分上下限

!!  open(unit=24,file='ckeckNC.txt')

  !*************变换坐标原点
  s=Point0(1,1);
  t=Point0(1,2);
  do i=1,Num0 
    Point1(i,1)=Point0(i,1)-s;
    Point1(i,2)=Point0(i,2)-t;
  enddo

  width=Point1(Num0,1)   !----剖面半宽

  s=width/real(NumIntC)
  do i=1,NumIntC
    IntC(i)=s*real(i)
  end do

  !-----开始计算积分
  TemI=20
  RomK=4

  IntNC=0.0;

  do i=1,NumIntC
    
    temC=IntC(i)

    !-----积分上下界(总)
    s=0.0;
    t=Pi/2.0;

    TemI1=TemI*i   !----小的积分段数目

    s1=(t-s)/real(TemI1)
    a=s
    do j=1,TemI1
      b=a+s1
      if( j==TemI1 ) b=t

      !-----开始Romberg积分
      call Slam_Romberg( RomK,a,b,temC,Num0,Point1,t1  )
      !---阶次，下限，上限
      IntNC(i)=IntNC(i)+t1

      a=b
    end do
    !----即对应于每一个Ct,计算出相应的入侵深度Ht
    IntNC(i)=IntNC(i)*2.0/Pi

!!    write(20,"(2(f20.9,1x))") IntC(i),IntNC(i)
  end do

  !-----上述得到IntNC(i)与IntC(i)的关系
  !----即  h=G(c)
  !----下面计算DG/Dc

  !-----样条插值的方式计算节点斜率
  allocate( length(NumIntC+1 ) )   !----累加弦长
  allocate( temp1(NumIntC+1,2),temp2(NumIntC+1,2) )
  allocate( funcx(NumIntC,4),funcy(NumIntC,4) )

  length=0.0;
  s=0.0;
  do i=1,NumIntC
    if(i==1) then
      t=sqrt( (IntNC(i))**2.0+(IntC(i))**2.0 )
    else
      t=sqrt( (IntNC(i)-IntNC(i-1) )**2.0+(IntC(i)-IntC(i-1) )**2.0 )
    end if
    s=s+t
    length(i+1)=s
  end do

  temp1=0.0;  
  temp1(:,1)=length(:);
  temp1(2:NumIntC+1,2)=IntC(1:NumIntC);  !---插x
  ds1=(temp1(2,2)-temp1(1,2) )/( temp1(2,1)-temp1(1,1) )
  dsn=(temp1(NumIntC+1,2)-temp1(NumIntC,2) )/( temp1(NumIntC+1,1)-temp1(NumIntC,1) )
  !-----插值中间函数Xt
  call spline3(NumIntC+1,temp1,ds1,dsn,funcx )

  temp2=0.0;
  temp2(:,1)=length(:);
  temp2(2:NumIntC+1,2)=IntNC(1:NumIntC);  !---插y
  ds1=(temp2(2,2)-temp2(1,2) )/( temp2(2,1)-temp2(1,1) )
  dsn=(temp2(NumIntC+1,2)-temp2(NumIntC,2) )/( temp2(NumIntC+1,1)-temp2(NumIntC,1) )
  !-----插值中间函数yt
  call spline3(NumIntC+1,temp2,ds1,dsn,funcy )


  do i=1,NumIntC+1
    s=length(i)
    if(i==1) then
      a1=funcx(i,1);  a2=funcx(i,2); 
      a3=funcx(i,3);  a4=funcx(i,4); 
      b1=funcy(i,1);  b2=funcy(i,2); 
      b3=funcy(i,3);  b4=funcy(i,4);
    else
      a1=funcx(i-1,1);  a2=funcx(i-1,2); 
      a3=funcx(i-1,3);  a4=funcx(i-1,4); 
      b1=funcy(i-1,1);  b2=funcy(i-1,2); 
      b3=funcy(i-1,3);  b4=funcy(i-1,4);
    end if

    s1=3.0*a1*s**2.0+2.0*a2*s+a3   !---dx/dt
    s2=3.0*b1*s**2.0+2.0*b2*s+b3   !---dy/dt

    DerNc(i)=s2/s1    !----dy/dx
  end do

  temRise=0.0;    !----选定Ct处对应液面抬升(从0开始)
  do i=1,NumIntC
    call SlamPoint_Interpolation(Num0,Point1(:,1),Point1(:,2),IntC(i),temRise(i)  )
    temRise(i)=temRise(i)-IntNC(i)
  end do


!!  do i=1,NumIntC
!!    write(24,"(4(f20.9,1x))") IntNC(i),IntC(i),temRise(i),DerNc(i+1)
!!  enddo

  



  return
end subroutine Slam_WagCondition



