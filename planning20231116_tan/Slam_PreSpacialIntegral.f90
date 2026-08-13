!输入：Dis船体瞬时摇荡位移；Numline：对砰击剖线循环的变量；NumPreNode：计算压力节点数
!PreNode:计算压力节点坐标; SlamPre:计算点压力
    
!SlamIntNumP:每条原始曲线理论最多插值节点（500）；SlamNumLP：每条剖线上砰击节点数
!SlamPoint:局部坐标系，砰击节点坐标；SlamNumIntegration：计算砰击压力的点数（200）
!SlamUr:左舷剖线节点模态位移；SlamRUr：左舷剖线节点模态位移梯度
!SlamWidth:剖线节点对应带宽；SlamNode；随船平动坐标系下的砰击节点坐标   
    
    
!输出：SForce：节点压力积分得到的剖线9个方向压力；SlineForce：单个砰击剖线上的砰击压力积分(在平动坐标系下的三个方向的力)    
    

!-----砰击压力在轮廓线上的空间积分
!-----包括耦合进模态位移
subroutine Slam_PreSpacialIntegral(Dis,Numline,NumPreNode,PreNode,SlamPre,SForce,SlineForce )

  use Slamming,only:SlamIntNumP,SlamNumLP,SlamPoint,SlamNumIntegration,SlamUr,&
                  &  SlamRUr,SlamAngle,SlamWidth,SlamNode
  
  use ShipHullVar,only:Cog
  use Constant,only:rou,NR
  use ArrayOperations
  implicit none

  real(8),dimension(6)::Dis    !---time时刻摇荡运动位移
  integer(4)::Numline
  integer(4)::NumPreNode   !----计算压力的节点数量
  real(8),dimension(0:NumPreNode,1:2 )::PreNode  !---压力节点坐标(局部坐标系下)
  real(8),dimension(0:NumPreNode )::SlamPre   !---压力

  real(8),dimension(NR)::SForce    !----砰击载荷(乘以了带宽，乘以了rou)
  
  real(8),dimension(3)::SlineForce  !---单个砰击剖线压力积分(未乘以带宽，乘以了rou)

  real(8),dimension(SlamNumIntegration,NR,3 )::NpreUr    !---压力积分点对应模态位移
  real(8),dimension(SlamNumIntegration,NR,3 )::NpreRUr   !---积分点模态位移，右舷
  real(8),dimension(SlamNumIntegration)::NperSwidth      !---积分点处对应切片带宽



  !---------------------------剖面几何信息
  integer(4)::Num0
  real(8),dimension(SlamIntNumP,2)::Point0
  real(8)::temAngle     !----平动坐标系下的剖面倾角
  real(8),dimension(SlamIntNumP )::temStripWidth    !---计算剖线节点对应切片带宽
  real(8),dimension(SlamIntNumP,3 )::temSlamNode

  !------剖线对应模态信息
  real(8),dimension(SlamIntNumP,NR,3 )::temSlamUr
  real(8),dimension(SlamIntNumP,NR,3 )::temSlamRUr

  real(8),dimension(SlamIntNumP,3 )::temPreSlamNode,temPreSlamRENode
!!  real(8),dimension(SlamIntNumP,NR,3,3 )::temSlamDUr
!!  real(8),dimension(SlamIntNumP,NR,3,3 )::temSlamRDUr

  real(8),dimension(3)::SlamNor    !---剖线方向向量
  real(8),dimension(3)::SegForceL,SegForceR
  
  integer(4)::i,j,k,ii,jj,kk
  real(8),dimension(2)::x1,x2,x3
  real(8)::s,t,s1,t1
  real(8)::P2,P3
  real(8),dimension(3)::vel,vel2

!***********************************************

  Num0=SlamNumLP(Numline )    !每条剖线节点数
  Point0=0.0;
  Point0(1:Num0,1:2)=SlamPoint(NumLine,1:Num0,1:2 );       !局部坐标系剖线节点坐标

  s=Point0(1,1);    !转换坐标系（将局部坐标系坐标原点移至第一个砰击节点位置）
  t=Point0(1,2);
  do i=1,Num0
    Point0(i,1)=Point0(i,1)-s
    Point0(i,2)=Point0(i,2)-t
  end do
  
  !------注意：此时的倾角还是平动坐标系下初始位置对应的倾角
  temAngle=SlamAngle(NumLine)
  !-----需要考虑船体摇荡运动的影响
  temAngle=temAngle-Dis(5)   !---考虑船体纵倾的影响(平动坐标系下，真实的剖线倾角)
!!  temAngle=SlamAngle(NumLine)


  temStripWidth=0.0;
  temStripWidth(1:Num0)=SlamWidth(NumLine,1:Num0 )    !----带宽

  temSlamNode=0.0;
  temSlamNode(1:Num0,1:3 )=SlamNode(NumLine,1:Num0,1:3  )   !----剖线对应平动坐标系下的节点坐标


  SlamNor=0.0                 !----砰击力的投影方向
  SlamNor(1)=cos(temAngle )
  SlamNor(2)=0.0;
  SlamNor(3)=sin(temAngle )


  !-------插值积分线段中心点对应的模态信息
  temSlamUr=0.0;   temSlamRUr=0.0;
!!  temSlamDUr=0.0;  temSlamRDUr=0.0;

  temSlamUr(1:Num0,:,: )=SlamUr(Numline,1:Num0,:,: )
  temSlamRUr(1:Num0,:,: )=SlamRUr(Numline,1:Num0,:,: )


  !-------插值积分线段中心点对应的模态信息
    NpreUr=0.0;      !压力积分对应的模态位移
    NpreRUr=0.0;      !压力积分对应的模态位移梯度
    NperSwidth=0.0;   !积分点处对应的切片带宽
    temPreSlamNode=0.0;     !压力计算点在随船平动坐标系下的位置
    temPreSlamRENode=0.0;    !右舷压力计算点在随船平动坐标系下的位置

    do i=1,NumPreNode    !对压力积分计算点数进行循环
      x2(:)=PreNode(i-1,:)     !压力积分节点坐标
      x3(:)=PreNode(i,:)
      !-----积分线段中心点(用于插值模态信息，以及插值切片宽度)
      x1(:)=(x2(:)+x3(:))/2.0

      do ii=7,NR    !----插值积分点位置处的模态信息(只需要插弹性体模态的)
        do jj=1,3
          call SlamPoint_Interpolation(Num0,Point0(1:Num0,1),temSlamUr(1:Num0,ii,jj ),x1(1),NpreUr(i,ii,jj)   )

          call SlamPoint_Interpolation(Num0,Point0(1:Num0,1),temSlamRUr(1:Num0,ii,jj ),x1(1),NpreRUr(i,ii,jj)   )
        end do
      end do

      !----插值积分点处的带宽信息
      call SlamPoint_Interpolation(Num0,Point0(1:Num0,1),temStripWidth(1:Num0),x1(1),NperSwidth(i)   )

      !----插值压力计算点对应随船平动坐标系下的三维节点位置(x,y,z)
      call SlamPoint_Interpolation( Num0,Point0(1:Num0,1),temSlamNode(1:Num0,1),x1(1),temPreSlamNode(i,1) )
      call SlamPoint_Interpolation( Num0,Point0(1:Num0,1),temSlamNode(1:Num0,2),x1(1),temPreSlamNode(i,2) )
      call SlamPoint_Interpolation( Num0,Point0(1:Num0,1),temSlamNode(1:Num0,3),x1(1),temPreSlamNode(i,3) )

      !----右舷对应节点
      temPreSlamRENode(i,1)=temPreSlamNode(i,1);
      temPreSlamRENode(i,2)=-temPreSlamNode(i,2);
      temPreSlamRENode(i,3)=temPreSlamNode(i,3);

    end do

    !-----注意：此时插值得到的积分点的各模态位移是平动坐标系下初始位置的
    !-----需要转换到瞬时位置：即考虑船体转动的影响(忽略横摇)
    vel=0.0;
    vel(1:3)=Dis(4:6);

    !------弹性体模态对应节点位移采用插值+旋转变换的方式得到
    do i=1,NumPreNode
      do ii=7,NR    !----插值积分点位置处的模态信息
            NpreUr(i,ii,1:3)=VectorL2G( NpreUr(i,ii,1:3),vel(1:3) )
            NpreRUr(i,ii,1:3)=VectorL2G( NpreRUr(i,ii,1:3),vel(1:3) )
      end do
    end do

    !------确定刚体运动的模态对应节点位移
    !------先将对应的平动坐标系下的节点变换到瞬时位置关于重心
    do i=1,NumPreNode
        vel(1:3)=temPreSlamNode(i,1:3)-Cog(1:3)
        vel(1:3)=VectorL2G( vel(1:3),Dis(4:6) )   !----瞬时位置关于瞬时重心

        !-----确定响应的刚体模态位移
        NpreUr(i,1,1)=1.0;   
        NpreUr(i,2,2)=1.0;
        NpreUr(i,3,3)=1.0;

        NpreUr(i,4,2)=-vel(3);  NpreUr(i,4,3)=vel(2);
        NpreUr(i,5,1)=vel(3);   NpreUr(i,5,3)=-vel(1);
        NpreUr(i,6,1)=-vel(2);  NpreUr(i,6,2)=vel(1);

        !-----右舷
        vel(1:3)=temPreSlamRENode(i,1:3)-Cog(1:3)
        vel(1:3)=VectorL2G( vel(1:3),Dis(4:6) )   !----瞬时位置关于瞬时重心

        !-----确定响应的刚体模态位移(右舷)
        NpreRUr(i,1,1)=1.0;   
        NpreRUr(i,2,2)=1.0;
        NpreRUr(i,3,3)=1.0;

        NpreRUr(i,4,2)=-vel(3);  NpreRUr(i,4,3)=vel(2);
        NpreRUr(i,5,1)=vel(3);   NpreRUr(i,5,3)=-vel(1);
        NpreRUr(i,6,1)=-vel(2);  NpreRUr(i,6,2)=vel(1);
    end do


    !****************************
    !-----开始砰击压力的积分
    SForce=0.0;
    SlineForce=0.0;
    do i=1,NumPreNode
      !----线单元端点坐标
      x2(:)=PreNode(i-1,:);
      x3(:)=PreNode(i,:);
      !----单元节点对应压力
      P2=SlamPre(i-1);
      P3=SlamPre(i);

      !-----压力积分（小单元）,此时，已经乘以了积分点处的带宽
      s=(P2+P3)*(x3(2)-x2(2) )/2.0*NperSwidth(i)   !----局部坐标系中，关于水平方向的压力积分
      t=(P2+P3)*(x3(1)-x2(1) )/2.0*NperSwidth(i)    !----局部坐标系中，关于竖直方向的压力积分

      !------------左舷
      SegForceL(1)=t*cos(temAngle )
      SegForceL(2)=-s
      SegForceL(3)=t*sin(temAngle )    !---注意符号
      do j=1,NR    
         SForce(j)=SForce(j)+dot_Product(NpreUr(i,j,1:3),SegForceL(1:3) )
      end do
      
      !------左舷：单根剖线压力积分
      SlineForce(1)=SlineForce(1)+(P2+P3)*(x3(1)-x2(1) )/2.0*cos(temAngle )
      SlineForce(2)=SlineForce(2)-(P2+P3)*(x3(2)-x2(2) )/2.0
      SlineForce(3)=SlineForce(3)+(P2+P3)*(x3(1)-x2(1) )/2.0*sin(temAngle )
      
      
      !------------右舷
      SegForceR(1)=t*cos(temAngle )
      SegForceR(2)=s
      SegForceR(3)=t*sin(temAngle )    !---注意符号
      do j=1,NR    
         SForce(j)=SForce(j)+dot_Product(NpreRUr(i,j,1:3),SegForceR(1:3) )
      end do

      !------左舷：单根剖线压力积分
      SlineForce(1)=SlineForce(1)+(P2+P3)*(x3(1)-x2(1) )/2.0*cos(temAngle )
      SlineForce(2)=SlineForce(2)+(P2+P3)*(x3(2)-x2(2) )/2.0
      SlineForce(3)=SlineForce(3)+(P2+P3)*(x3(1)-x2(1) )/2.0*sin(temAngle )     
      
    end do

    !-----此处乘以密度，得到完整的剖面切片砰击载荷
    SForce=SForce*rou
    
    !-----乘以密度，得到单根剖线上的完整的压力积分
    SlineForce=SlineForce*rou

  return
end subroutine Slam_PreSpacialIntegral