subroutine Shiphull
    use verification
    use ShipHullVar
    use Constant
    use Slamming
    use ArrayOperations
    use PanelGeometry
    
    use lift
    
    implicit none
    
    integer(4)::i,j,k,jj,j1,iii,jjj
    real(8)::deltax
    real(8),dimension(3)::dp1,dp2
    
    real(8)::TmpMij(1:6,1:6),TmpVec(1:3)          !构造质量矩阵的临时变量
    
    real(8),dimension(6)::nor                     !求解初始静浮力时的变量
    real(8),dimension(3)::vel
    
    real(8),dimension(NR,NR)::temHerm1,temHerm2,temHerm3
    
    real(8),dimension(1:NR,1:3)::eu,eu1,eu2,eu3,eu4
    real(8),dimension(1:NR,1:3,1:3)::deu,deu1,deu2,deu3,deu4
    
    
    real(8),dimension(Nwh,Nl,3)::Node0            !划分好的湿表面网格在用户坐标系下关于重心的位置 
    
    
    real(8),allocatable,dimension(:,:)::temSlamNode     !-----有关砰击计算的临时数据存储变量
    
    
    integer(4)::Nwh0,Nl0
    
    !******计算剖面载荷相关变量20230425******
    real(8),dimension(NR,6)::fdrm1,fdrmN      !----样条插值边界条件
    real(8),dimension(NBSECT,NR,6)::Ns,Nds,Ndds
    
    !******计算计算与梁振动运动相关的模态，用于考虑水动升力弹性20230426******
    real(8),dimension(NR,6)::Mfdrm1,MfdrmN
    real(8),dimension(SlamNumLine,NR,6)::MNs,MNds,MNdds
    
    !*******与机翼升力相关的变量20230506******
    real(8)::Stemp_angle(Swing_num),Mtemp_angle(Mwing_num)        !坐标变换临时变量
    
    
!**********前面是变量定义*********** 
    open(unit=11,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//'.che')
    
    
    Volx=0.0; Voly=0.0; Volz=0.0
    Cobx=0.0; Coby=0.0; Cobz=0.0
    
    CoB(1:3)=0.0;WetArea=0.0;
    
    do j=1,NL
        do i=1,Nwh-1
        if(j<Nl)then
        x1=Node(i,j,1);y1=Node(i,j,2);z1=Node(i,j,3)
        x2=Node(i+1,j,1);y2=Node(i+1,j,2);z2=Node(i+1,j,3)
        x3=Node(i+1,j+1,1);y3=Node(i+1,j+1,2);z3=Node(i+1,j+1,3)
        x4=Node(i,j+1,1);y4=Node(i,j+1,2);z4=Node(i,j+1,3)
        else
        x1=Node(i,j,1);y1=Node(i,j,2);z1=Node(i,j,3)
        x2=Node(i+1,j,1);y2=Node(i+1,j,2);z2=Node(i+1,j,3)
        x3=Node(i+1,1,1);y3=Node(i+1,1,2);z3=Node(i+1,1,3)
        x4=Node(i,1,1);y4=Node(i,1,2);z4=Node(i,1,3)      
        end if
        
        DxDxi=x2-x1;DyDxi=y2-y1;DzDxi=z2-z1
        DxDet=x4-x1;DyDet=y4-y1;DzDet=z4-z1
        hs1_x(1) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(1) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(1) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(1)=(hs1_x(1)**2+hs1_y(1)**2+hs1_z(1)**2)**0.5

        DxDxi=x3-x2;DyDxi=y3-y2;DzDxi=z3-z2
        DxDet=x1-x2;DyDet=y1-y2;DzDet=z1-z2
        hs1_x(2) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(2) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(2) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(2)=(hs1_x(2)**2+hs1_y(2)**2+hs1_z(2)**2)**0.5
  
        DxDxi=x4-x3;DyDxi=y4-y3;DzDxi=z4-z3
        DxDet=x2-x3;DyDet=y2-y3;DzDet=z2-z3
        hs1_x(3) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(3) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(3) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(3)=(hs1_x(3)**2+hs1_y(3)**2+hs1_z(3)**2)**0.5

        DxDxi=x1-x4;DyDxi=y1-y4;DzDxi=z1-z4
        DxDet=x3-x4;DyDet=y3-y4;DzDet=z3-z4
        hs1_x(4) = DyDxi*DzDet - DyDet*DzDxi
        hs1_y(4) = DzDxi*DxDet - DzDet*DxDxi
        hs1_z(4) = DxDxi*DyDet - DxDet*DyDxi            
        hs1(4)=(hs1_x(4)**2+hs1_y(4)**2+hs1_z(4)**2)**0.5
        
        do jj=1,2
            do j1=1,2
            ph1=(1.+cor2(j1))*(1.+cor2(jj))/4.
		    ph2=(1.-cor2(j1))*(1.+cor2(jj))/4.
		    ph3=(1.-cor2(j1))*(1.-cor2(jj))/4.
		    ph4=(1.+cor2(j1))*(1.-cor2(jj))/4.

            xq=ph1*x1+ph2*x2+ph3*x3+ph4*x4
            yq=ph1*y1+ph2*y2+ph3*y3+ph4*y4
            zq=ph1*z1+ph2*z2+ph3*z3+ph4*z4
            hs_x=(ph1*hs1_x(1)+ph2*hs1_x(2)+ph3*hs1_x(3)+ph4*hs1_x(4))/4.        
            hs_y=(ph1*hs1_y(1)+ph2*hs1_y(2)+ph3*hs1_y(3)+ph4*hs1_y(4))/4.         
            hs_z=(ph1*hs1_z(1)+ph2*hs1_z(2)+ph3*hs1_z(3)+ph4*hs1_z(4))/4.
            hs=(ph1*hs1(1)+ph2*hs1(2)+ph3*hs1(3)+ph4*hs1(4))/4.

            hs_x=-hs_x
            hs_y=-hs_y       
            hs_z=-hs_z       
            Volx=Volx-ww(j1)*ww(jj)*hs_x*xq
            Voly=Voly-ww(j1)*ww(jj)*hs_y*yq
            Volz=Volz-ww(j1)*ww(jj)*hs_z*zq
        
            Cobx=Cobx-ww(j1)*ww(jj)*hs_x*xq*xq
            Coby=Coby-ww(j1)*ww(jj)*hs_y*yq*yq
            Cobz=Cobz-ww(j1)*ww(jj)*hs_z*zq*zq
        
            WetArea=WetArea+ww(j1)*ww(jj)*hs
                
                
            end do    
        end do
        
  
        end do
    end do
    
    !默认对称
    Coby=0.0
    !湿表面积
    WetArea=WetArea
    !排水体积
    vol=(volX+volY+volZ)/3.0
    !浮心
    COB(1)=CobX/vol/2.0      
    COB(2)=CobY/vol/2.0
    COB(3)=CobZ/vol/2.0
    deltax=COB(1)
    
    !**部分变量开始转换到随船平动坐标系中
    !全船网格，湿表面网格只需x方向平移，前两步在waterline模块里完成
    Node(:,:,1)=Node(:,:,1)-deltax
    Nodeb(:,:,1)=Nodeb(:,:,1)-deltax
    InstBreakNode(:,1)=InstBreakNode(:,1)-deltax
    
    !浮心坐标转换到随船平动坐标系下
    COB(1)=0.0
    !-------得到平动坐标下下的重心坐标
    Cog(3)=Cog(3)-Ta
    Cog(1:3)=VectorL2G(Cog(1:3),d_mc(4:6,1)) 
    Cog(1)=Cog(1)-deltax
    Cog(2)=0.0;    !---默认对称
    
    !将分段质心转换到随船平动坐标系下
    PointCor2(3,:)=PointCor2(3,:)-Ta
    do i=1,NITEM
        PointCor2(1:3,i)=VectorL2G(PointCor2(1:3,i),d_mc(4:6,1) )
    end do
    PointCor2(1,:)=PointCor2(1,:)-deltax        !----平动坐标系下，分段质心
    !-----随船平动坐标系下，计算分段质心相对于整体重心位置
    PointCor2(1,:)=PointCor2(1,:)-cog(1)
    PointCor2(2,:)=PointCor2(2,:)-cog(2)
    PointCor2(3,:)=PointCor2(3,:)-cog(3)      !----平动坐标系，关于质心

    !------将剖面形心坐标SecVec2转化到随船平动坐标系下
    SecVec2(3,:)=SecVec2(3,:)-Ta
    do i=1,NBSECT
        SecVec2(1:3,i)=VectorL2G(SecVec2(1:3,i),d_mc(4:6,1) )
    end do
    SecVec2(1,:)=SecVec2(1,:)-deltax
    
    !*******将读取的bdf文件网格节点转换到平动坐标系下
    !-------此时认为选择非线性计算方法，并且采用读取bdf文件的方式
    if( trim(adjustl(Non_Linear))=='NL'.and.nonlinearCtrl==1 ) then
        
        bdfnode(:,3)=bdfnode(:,3)-Ta;
        do i=1,bdfnum_node
            bdfnode(i,1:3)=VectorL2G(bdfnode(i,1:3),d_mc(4:6,1) )
        end do
        bdfnode(:,1)=bdfnode(:,1)-deltax
        
    end if
    
    
    !***弹性网格转换***这里先不处理
    Node0(:,:,1)=Node(:,:,1)-Cog(1)
    Node0(:,:,2)=Node(:,:,2)-Cog(2)
    Node0(:,:,3)=Node(:,:,3)-Cog(3)
    
    !-----开始计算弹性体模态中，网格节点的位移，位移梯度
    !----- 得到各模态对应网格节点的位移，位移梯度
    Nwh0=Nwh
    
    Nl0=NL
    
    call Emesh(Nwh0,Nl0,Node0)
    
   

    !*******************机翼升力相关变量处理***
    if( trim(adjustl(Airlift))=='YES' )then
    !尾翼
    !计算中弧线节点坐标（随动坐标系）    
    do i=1,Swing_num
        do j=1,Spoint_num(i)
            Swing_point(i,j,1)=( SUpoint(i,j,1)+SDpoint(i,j,1) )/2.0
            Swing_point(i,j,2)=( SUpoint(i,j,2)+SDpoint(i,j,2) )/2.0
            Swing_point(i,j,3)=( SUpoint(i,j,3)+SDpoint(i,j,3) )/2.0            
        end do
    end do
    
    !将中弧线上节点坐标转到局部坐标系
    
    !先将坐标原点平移到第一个节点位置
    do i=1,Swing_num
        do j=1,Spoint_num(i)
            
            Stemp_point(i,j,1)=Swing_point(i,j,1)-Swing_point(i,1,1)
            Stemp_point(i,j,2)=Swing_point(i,j,2)-Swing_point(i,1,2)
            Stemp_point(i,j,3)=Swing_point(i,j,3)-Swing_point(i,1,3)            
        end do        
    end do
    
    !计算夹角
    do i=1,Swing_num
        Stemp_angle(i)=atan( Stemp_point(i,Spoint_num(i),3)/Stemp_point(i,Spoint_num(i),1) )        
    end do
    
    !最后旋转到局部坐标系
    do i=1,Swing_num
        do j=1,Spoint_num(i)
            Swing_node(i,j,1)=Stemp_point(i,j,1)*cos( Stemp_angle(i) )+Stemp_point(i,j,3)*sin( Stemp_angle(i) )
            
            Swing_node(i,j,2)=Stemp_point(i,j,3)*cos( Stemp_angle(i) )-Stemp_point(i,j,1)*sin( Stemp_angle(i) )
        end do
    end do
    
    !计算弦长
    allocate( Sca(Swing_num) )
    Sca=0.0;
    do i=1,Swing_num
        Sca(i)=Swing_node(i,Spoint_num(i),1)-Swing_node(i,1,1)        
    end do
    
    !中部机翼
    !计算中弧线节点坐标（随动坐标系） 
    do i=1,Mwing_num
        do j=1,Mpoint_num(i)
            Mwing_point(i,j,1)=( MUpoint(i,j,1)+MDpoint(i,j,1) )/2.0
            Mwing_point(i,j,2)=( MUpoint(i,j,2)+MDpoint(i,j,2) )/2.0
            Mwing_point(i,j,3)=( MUpoint(i,j,3)+MDpoint(i,j,3) )/2.0
        end do     
    end do
    
    !将中弧线上节点坐标转到局部坐标系
    !先将坐标原点平移到第一个节点位置
    do i=1,Mwing_num
        do j=1,Mpoint_num(i)
            Mtemp_point(i,j,1)=Mwing_point(i,j,1)-Mwing_point(i,1,1)
            Mtemp_point(i,j,2)=Mwing_point(i,j,2)-Mwing_point(i,1,2)
            Mtemp_point(i,j,3)=Mwing_point(i,j,3)-Mwing_point(i,1,3)
        end do
    end do
    
    !计算夹角
    do i=1,Mwing_num
        Mtemp_angle(i)=atan( Mtemp_point(i,Mpoint_num(i),3)/Mtemp_point(i,Mpoint_num(i),1) )        
    end do
    
    !最后旋转到局部坐标系
    do i=1,Mwing_num
        do j=1,Mpoint_num(i)
            Mwing_node(i,j,1)=Mtemp_point(i,j,1)*cos( Mtemp_angle(i) )+Mtemp_point(i,j,3)*sin( Mtemp_angle(i) )
        
            Mwing_node(i,j,2)=Mtemp_point(i,j,3)*cos( Mtemp_angle(i) )-Mtemp_point(i,j,1)*sin( Mtemp_angle(i) )
        end do
    end do
        
    !计算弦长
    allocate( Mca(Mwing_num) )
    Mca=0.0;
    
    do i=1,Mwing_num
        Mca(i)=Mwing_node(i,Mpoint_num(i),1)-Mwing_node(i,1,1)        
    end do
    
 
    end if
    
    
    
    !----水动升力计算节点转换到随船平动坐标系下(这里先不考虑弹性)
    if( trim(adjustl(Liftctrl))=='MLM' )then
        
        !********以下是与模态节点位移有关的全局变量
        allocate( SlamUr(SlamNumLine,SlamIntNumP,NR,3 ) )     !---左舷剖线节点模态位移
        allocate( SlamRUr(SlamNumLine,SlamIntNumP,NR,3 ) )     !---右舷剖线节点模态位移
        allocate( SlamDur(SlamNumLine,SlamIntNumP,NR,3,3 ) )  !---左舷剖线节点模态位移梯度
        allocate( SlamRDur(SlamNumLine,SlamIntNumP,NR,3,3 ) )   !---右舷剖线节点模态位移梯度
        
        
        !先将剖线SlamNode转换为随船平动坐标系下
        do i=1,SlamNumLine      !对典型砰击节点数进行循环（对应位置有可能砰击剖线未生成）
            if( SlamLineCase(i)==0 ) cycle   !-----进行相关计算时，必须判断Slamcase是否允许进行计算

            !-----开始进行坐标转换(计算初始节点已处于用户坐标系下，关于重心)
            k=SlamNumLP(i)
            do j=1,k      !将细化后的砰击节点坐标（用户坐标系相对于重心），转到随船平动坐标系中
                SlamNode(i,j,1:3)=VectorL2G( SlamNode(i,j,1:3),d_mc(4:6,1)  )

                SlamNode(i,j,1)=SlamNode(i,j,1)+Cog(1)
                SlamNode(i,j,2)=SlamNode(i,j,2)+Cog(2)
                SlamNode(i,j,3)=SlamNode(i,j,3)+Cog(3)
            end do

            !-----典型节点(2,3号)也转到随船平动坐标系中
            do j=2,3
                SlamIniType(i,j,1:3)=VectorL2G( SlamIniType(i,j,1:3),d_mc(4:6,1)  )
                SlamIniType(i,j,1)=SlamIniType(i,j,1)+Cog(1)
                SlamIniType(i,j,2)=SlamIniType(i,j,2)+Cog(2)
                SlamIniType(i,j,3)=SlamIniType(i,j,3)+Cog(3)
            end do
        end do

        !-----将典型节点转换到随船平动坐标系下(1号，原本是在用户坐标系下)
        !-----注意：一号典型节点在后续中会用来判断是否发生砰击，因此一经确认不会再改变
        !           并且一直存在
        do i=1,SlamNumLine
            j=1;
            SlamIniType(i,j,3)=SlamIniType(i,j,3)-Ta
            SlamIniType(i,j,1:3)=VectorL2G( SlamIniType(i,j,1:3),d_mc(4:6,1)  )
            SlamIniType(i,j,1)=SlamIniType(i,j,1)-deltax
        end do
               
        
        
        
        !20231115修改，添加输出船中纵剖线各点与波面的相对位置
        if( trim(adjustl(SlamRelaMctr))=='YES' )then
            allocate( SlamShipRMP(SlamNumShipRMP,3),SlamWaveRMP(SlamNumWaveRMP,3) )
            SlamShipRMP=0.0;  SlamWaveRMP=0.0;
            
            !开始坐标转换(船体纵剖线上的点)
            do i=1,SlamNumShipRMP
                dp1=0.0;
            
                dp1(1:3)=SlamShipIniRMP(i,1:3);
                dp1(3)=dp1(3)-Ta;
                dp1(1:3)=VectorL2G( dp1(1:3),d_mc(4:6,1)  );   !----调好了纵倾角
            
                SlamShipIniRMP(i,1:3)=dp1(1:3);    !----确定第一种节点坐标
            
                dp1(1)=dp1(1)-deltax    !---减去浮心坐标
                SlamShipRMP(i,1:3)=dp1(1:3);      !----确定第二种节点坐标
            end do
            
            !波面上的点进行坐标转换
            do i=1,SlamNumWaveRMP
                dp1=0.0;
            
                dp1(1:3)=SlamWaveIniRMP(i,1:3);
                dp1(3)=dp1(3)-Ta;
                dp1(1:3)=VectorL2G( dp1(1:3),d_mc(4:6,1)  );   !----调好了纵倾角            
                dp1(3)=0.0;
            
                SlamWaveIniRMP(i,1:3)=dp1(1:3);   !----确定第一种节点的坐标
            
                dp1(1)=dp1(1)-deltax    !---减去浮心坐标
                SlamWaveRMP(i,1:3)=dp1(1:3);   !----确定第二种节点坐标
            end do
 
        end if
        
        
        
        
        !*******************确定砰击计算切片有关变量
        !--------计算砰击曲线切片边界点
        !-----由数据读取/当前计算出来
        allocate( SlamStripTyp(SlamNumLine,2,3 ) )   !----切片边界点(典型点，定位)，指的是厚度左右边界
        allocate( SlamStripBexist(SlamNumLine,2) )   !----切片边界划分与否状态0未划分，1划分
        allocate( SlamStripBAngle(SlamNumLine,2) )   !----边界倾角，随船平动坐标系，与SlamAngle类似
        allocate( SlamWidth(SlamNumLine,SlamIntNumP ) )  !---剖线节点对应带宽,为后续积分点提供插值
        allocate( SlamStripBNode(SlamNumLine,SlamIntNumP,2,3) )  !---切片左右两侧边界节点坐标
        
        SlamStripTyp=0.0;    !----一经确认，不会再改变
        do i=1,SlamNumLine
          if(i==1) then
            dp1(1:3)=SlamIniType(i,1,1:3)
            dp2(1:3)=SlamIniType(i+1,1,1:3)
          elseif(i==SlamNumLine) then
            dp1(1:3)=SlamIniType(i-1,1,1:3)
            dp2(1:3)=SlamIniType(i,1,1:3)
          else
            dp1(1:3)=SlamIniType(i-1,1,1:3)
            dp2(1:3)=SlamIniType(i+1,1,1:3)
          end if

          SlamStripTyp(i,1,:)=(dp1(:)+SlamIniType(i,1,:) )/2.0      !确定切片边界点（典型点）坐标（左）
          SlamStripTyp(i,2,:)=(dp2(:)+SlamIniType(i,1,:) )/2.0      !右

        end do
        
        !----确认边界倾角
        SlamStripBexist=0;
        SlamStripBAngle=0.0;
        SlamWidth=0.0;
        SlamStripBNode=0.0;
        
        do i=1,SlamNumLine    !对典型砰击节点进行循环
      
          if( SlamLineCase(i)==0 ) cycle
          !*******************
          !----注：在后续实时修改中，需要根据是否发生砰击的情况来确认

          !为什么全局变量还可以传参，这没有问题吗
          call Slam_StripWidth( i,SlamWidth(i,:),SlamStripBNode(i,:,:,:) )

        end do
        
        
        allocate( temSlamNode(SlamIntNumP,3) )
        !-----注意，实时计算时，一定要初始化
        SlamUr=0.0;    SlamRUr=0.0;
        SlamDur=0.0;   SlamRDur=0.0;
        
        do i=1,SlamNumLine     !对典型砰击节点进行循环
            if( SlamLineCase(i)==0 ) cycle
            
            k=SlamNumLP(i)
            temSlamNode=0.0    !----先将临时节点转换到关于重心的位置（随船平动坐标系下关于重心）
            
            do j=1,k
                temSlamNode(j,1)=SlamNode(i,j,1)-Cog(1)
                temSlamNode(j,2)=SlamNode(i,j,2)-Cog(2)
                temSlamNode(j,3)=SlamNode(i,j,3)-Cog(3) 
            end do
            
            !-----注意，实时计算时，一定要初始化
            SlamUr(i,:,:,:)=0.0;    SlamRUr(i,:,:,:)=0.0;    !左右舷剖线节点模态位移
            SlamDur(i,:,:,:,:)=0.0; SlamRDur(i,:,:,:,:)=0.0;  !左右舷剖线节点模态位移梯度
            
            !k:剖线上的细化节点数；temSlamNode：随船平动坐标系下关于重心的砰击节点坐标
            !-----再进行弹性体模态(包含刚体)节点位移的计算以及位移梯度的计算
            
            call SlamEmesh(k,temSlamNode(1:k,1:3),SlamUr(i,1:k,:,:),SlamDur(i,1:k,:,:,:),&
                   &  SlamRUr(i,1:k,:,:),SlamRDur(i,1:k,:,:,:) )
            
         
        end do
 
    end if
        
    !!***check***    
    write(11,"(A)") '# 重心坐标(随船平动坐标系)'
    write(11,"(3(f12.3,1x))") Cog(1),Cog(2),Cog(3)
    write(11,"(A)") '# 浮心坐标(随船平动坐标系)'
    write(11,"(3(f12.3,1x))") Cob(1),Cob(2),Cob(3)
    write(11,"(A)") '# 排水体积(m^3))'
    write(11,"(e15.6)") vol
    
    
    !*******20230420重心浮心排水体积检查问题*******
    
    
    
    !构造结构质量矩阵
    if(trim(adjustl(Masssolve))=='WHOLE')then
        TmpMij=0.0;
        Mass=1000*Mass
        TmpMij(1,1)=MASS; TmpMij(2,2)=MASS; TmpMij(3,3)=MASS
        TmpMij(4,4)=MASS*TotalI11**2 ; TmpMij(5,5)=Mass*TotalI22**2;  TmpMij(6,6)=Mass*TotalI33**2
        TmpMij(4,6)=Mass*TotalI13**2 ; TmpMij(6,4)=Mass*TotalI13**2
  
    else

        TmpMij=0.0
        TmpMij(1,1)=sum(MM(:))
        TmpMij(2,2)=TmpMij(1,1)
        TmpMij(3,3)=TmpMij(1,1)
        
        do i=1,NITEM
            
            !-----质量模型坐标系(用户坐标系下计算得到整船质量矩阵)
            TmpVec(1)=PointCor(1,I)
            TmpVec(2)=PointCor(2,I)
            TmpVec(3)=PointCor(3,I)                
  
            !----20200325  本次修改参考basic
            TmpMij(4,4)=TmpMij(4,4)+MM(i)*IX_R(i)**2.0+MM(i)*(TmpVec(2)**2.0+TmpVec(3)**2.0)   !---横摇转动惯量
            TmpMij(5,5)=TmpMij(5,5)+MM(i)*(x12(i,2)-x12(i,1))**2.0/12.0+MM(i)*( TmpVec(1)**2+TmpVec(3)**2 )
            TmpMij(6,6)=TmpMij(6,6)+MM(i)*(x12(i,2)-x12(i,1))**2.0/12.0+MM(i)*( TmpVec(1)**2+TmpVec(2)**2 )

            TmpMij(4,5)=TmpMij(4,5)-MM(i)*(  TmpVec(1)* TmpVec(2) )
            TmpMij(5,4)=TmpMij(4,5)
            TmpMij(4,6)=TmpMij(4,6)-MM(i)*(  TmpVec(1)* TmpVec(3) )
            TmpMij(6,4)=TmpMij(4,6)
            TmpMij(5,6)=TmpMij(5,6)-MM(i)*(  TmpVec(2)* TmpVec(3) )
            TmpMij(6,5)=TmpMij(5,6)
        end do

    end if

   allocate( Mij( MonDim,MonDim ) ) 
   Mij=0.0;
   
   Mij(1:6,1:6)=tmpMij(1:6,1:6)      !6阶刚体结构质量矩阵  
   
   StruM(1:6,1:6)=tmpMij(1:6,1:6)    !9阶弹性体结构质量矩阵
   
   do i=7,NR
        StruC(i,i)=Ihome(i)**2.0*StruM(i,i)   !----结构刚度矩阵
        StruB(i,i)=2.0*StruM(i,i)*Ihome(i)*ceb(i)   !---结构阻尼矩阵(类似临界阻尼系数)
   end do
   
 
   
   
   !open(unit=1001,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.txt')
   !call iutmp(1001)
   !do i=1,NR
   !     read(1001,*) ( StruB(i,j),j=1,NR )
   !end do
   !do i=1,NR
   !    read(1001,*)  ( StruC(i,j),j=1,NR )
   !end do
   !
   !close(1001)
   
   
   
   
   
   
   !****检查****
    write(11,"(A)") '#   结构质量矩阵'
    do i=1,NR
    do j=1,NR
        write(11,"(E15.4,1x,\)",advance='NO') StruM(i,j)
    end do
    write(11,"(/)") 
    end do

    write(11,"(A)") '#   结构阻尼矩阵'
    do i=1,NR
    do j=1,NR
        write(11,"(E15.4,1x,\)",advance='NO') StruB(i,j)
    end do
    write(11,"(/)") 
    end do

    write(11,"(A)") '#   结构刚度矩阵'
    do i=1,NR
    do j=1,NR
        write(11,"(E15.4,1x,\)",advance='NO') StruC(i,j)
    end do
    write(11,"(/)") 
    end do
   
   !********20230420结构质量矩阵、刚度矩阵、阻尼矩阵检查没问题*******
   
     
   
   !****如果是弹性体，质量矩阵应该是9阶****
   !****刚体时不考虑剖面载荷，因为达朗贝尔原理比较麻烦，等到弹性体构造完成之后再考虑剖面载荷****

   
   !**为计算原始重力与浮力的差值，这里需要计算原始重力及原始原始浮力,另外还搭建静水恢复力矩阵
   
   !原始重力
   Initial_Mg=0.0;
   
   allocate( Initial2_Mg(NR) )
   Initial2_Mg=0.0;      !20240109

   do i=1,NITEM
       Initial_Mg(3)=Initial_Mg(3)-MM(i)*g0;
       Initial_Mg(4)=Initial_Mg(4)-MM(i)*g0*PointCor2(2,i)
       Initial_Mg(5)=Initial_Mg(5)+MM(i)*g0*PointCor2(1,i)       
   end do
   
   Initial2_Mg(3)=Initial_Mg(3)       !20240109
   
   
   !--------Xg,Yg,Zg  也是平动坐标系下的重心坐标
    Zg=Cog(3)
    Xg=Cog(1)
    Yg=Cog(2)
    
   
   
   !初始浮力
   Initial_Fs=0.0;
   
   temHerm1=0.0; temHerm2=0.0;  temHerm3=0.0;
   
    HRM=0.0   !----静水恢复力矩阵(6自由度刚体运动的)
    HERM=0.0  !----包含刚体运动模态在内的弹性体模态
  
   
   do j=1,NL
       do i=1,Nwh-1
           if (j<Nl) then
                x1=Node(i,j,1);y1=Node(i,j,2);z1=Node(i,j,3)
                x2=Node(i+1,j,1);y2=Node(i+1,j,2);z2=Node(i+1,j,3)
                x3=Node(i+1,j+1,1);y3=Node(i+1,j+1,2);z3=Node(i+1,j+1,3)
                x4=Node(i,j+1,1);y4=Node(i,j+1,2);z4=Node(i,j+1,3)
                
                !----各模态位移
                eu1(1:NR,1:3)=Ur(i,j,1:NR,1:3);     eu2(1:NR,1:3)=Ur(i+1,j,1:NR,1:3); 
                eu3(1:NR,1:3)=Ur(i+1,j+1,1:NR,1:3); eu4(1:NR,1:3)=Ur(i,j+1,1:NR,1:3);
                !----各模态位移梯度
                deu1(1:NR,1:3,1:3)=dur(i,j,1:NR,1:3,1:3);
                deu2(1:NR,1:3,1:3)=dur(i+1,j,1:NR,1:3,1:3);
                deu3(1:NR,1:3,1:3)=dur(i+1,j+1,1:NR,1:3,1:3);
                deu4(1:NR,1:3,1:3)=dur(i,j+1,1:NR,1:3,1:3);
        
 
           else
                x1=Node(i,j,1);y1=Node(i,j,2);z1=Node(i,j,3)
                x2=Node(i+1,j,1);y2=Node(i+1,j,2);z2=Node(i+1,j,3)
                x3=Node(i+1,1,1);y3=Node(i+1,1,2);z3=Node(i+1,1,3)
                x4=Node(i,1,1);y4=Node(i,1,2);z4=Node(i,1,3)
                
                !----各模态位移
                eu1(1:NR,1:3)=Ur(i,j,1:NR,1:3);     eu2(1:NR,1:3)=Ur(i+1,j,1:NR,1:3); 
                eu3(1:NR,1:3)=Ur(i+1,1,1:NR,1:3);   eu4(1:NR,1:3)=Ur(i,1,1:NR,1:3);               
                !----各模态位移梯度
                deu1(1:NR,1:3,1:3)=dur(i,j,1:NR,1:3,1:3);
                deu2(1:NR,1:3,1:3)=dur(i+1,j,1:NR,1:3,1:3);
                deu3(1:NR,1:3,1:3)=dur(i+1,1,1:NR,1:3,1:3);
                deu4(1:NR,1:3,1:3)=dur(i,1,1:NR,1:3,1:3);      

               
           end if
           
            DxDxi=x2-x1;DyDxi=y2-y1;DzDxi=z2-z1
            DxDet=x4-x1;DyDet=y4-y1;DzDet=z4-z1
            hs1_x(1) = DyDxi*DzDet - DyDet*DzDxi
            hs1_y(1) = DzDxi*DxDet - DzDet*DxDxi
            hs1_z(1) = DxDxi*DyDet - DxDet*DyDxi            
            hs1(1)=(hs1_x(1)**2+hs1_y(1)**2+hs1_z(1)**2)**0.5
            
            DxDxi=x3-x2;DyDxi=y3-y2;DzDxi=z3-z2
            DxDet=x1-x2;DyDet=y1-y2;DzDet=z1-z2
            hs1_x(2) = DyDxi*DzDet - DyDet*DzDxi
            hs1_y(2) = DzDxi*DxDet - DzDet*DxDxi
            hs1_z(2) = DxDxi*DyDet - DxDet*DyDxi            
            hs1(2)=(hs1_x(2)**2+hs1_y(2)**2+hs1_z(2)**2)**0.5
            
            DxDxi=x4-x3;DyDxi=y4-y3;DzDxi=z4-z3
            DxDet=x2-x3;DyDet=y2-y3;DzDet=z2-z3
            hs1_x(3) = DyDxi*DzDet - DyDet*DzDxi
            hs1_y(3) = DzDxi*DxDet - DzDet*DxDxi
            hs1_z(3) = DxDxi*DyDet - DxDet*DyDxi            
            hs1(3)=(hs1_x(3)**2+hs1_y(3)**2+hs1_z(3)**2)**0.5
            
            DxDxi=x1-x4;DyDxi=y1-y4;DzDxi=z1-z4
            DxDet=x3-x4;DyDet=y3-y4;DzDet=z3-z4
            hs1_x(4) = DyDxi*DzDet - DyDet*DzDxi
            hs1_y(4) = DzDxi*DxDet - DzDet*DxDxi
            hs1_z(4) = DxDxi*DyDet - DxDet*DyDxi            
            hs1(4)=(hs1_x(4)**2+hs1_y(4)**2+hs1_z(4)**2)**0.5
            
            do jj=1,2
                do j1=1,2
                    
                   ph1=(1.+cor2(j1))*(1.+cor2(jj))/4.
                   ph2=(1.-cor2(j1))*(1.+cor2(jj))/4.
		           ph3=(1.-cor2(j1))*(1.-cor2(jj))/4.
		           ph4=(1.+cor2(j1))*(1.-cor2(jj))/4.
                    
                   xq=ph1*x1+ph2*x2+ph3*x3+ph4*x4
	               yq=ph1*y1+ph2*y2+ph3*y3+ph4*y4
                   zq=ph1*z1+ph2*z2+ph3*z3+ph4*z4
                   hs_x=(ph1*hs1_x(1)+ph2*hs1_x(2)+ph3*hs1_x(3)+ph4*hs1_x(4))/4.        
                   hs_y=(ph1*hs1_y(1)+ph2*hs1_y(2)+ph3*hs1_y(3)+ph4*hs1_y(4))/4.         
                   hs_z=(ph1*hs1_z(1)+ph2*hs1_z(2)+ph3*hs1_z(3)+ph4*hs1_z(4))/4.
                   hs=(ph1*hs1(1)+ph2*hs1(2)+ph3*hs1(3)+ph4*hs1(4))/4.
                   
                   eu(1:NR,1:3)=ph1*eu1(1:NR,1:3)+ph2*eu2(1:NR,1:3)+ph3*eu3(1:NR,1:3)+ph4*eu4(1:NR,1:3)
                   deu(:,:,:)=ph1*deu1(:,:,:)+ph2*deu2(:,:,:)+ph3*deu3(:,:,:)+ph4*deu4(:,:,:)
                                      
                   
                   hs_x=-hs_x
                   hs_y=-hs_y       
                   hs_z=-hs_z
                   
                   nor(1)=hs_x
                   nor(2)=hs_y       
                   nor(3)=hs_z
                   vel(1)=xq-Xg;vel(2)=yq-Yg;vel(3)=zq-Zg         !----从船体重心延伸出的矢径
                   nor(4:6)=R_Rcross( vel(1:3),nor(1:3) )
                   
                   do k=1,6
                       Initial_Fs(k)=Initial_Fs(k)-ww(j1)*ww(jj)*nor(k)*zq    !---船体初始静浮力  
                   end do
                   
                   !----静水力回复矩阵(刚体的，用于检验)
                    HRM(3,3)=HRM(3,3)+ww(j1)*ww(jj)*hs_z
                    HRM(3,5)=HRM(3,5)-ww(j1)*ww(jj)*hs_z*xq
                    HRM(4,4)=HRM(4,4)+ww(j1)*ww(jj)*hs_z*yq*yq
                    HRM(5,5)=HRM(5,5)+ww(j1)*ww(jj)*hs_z*xq*xq
                    
                    do iii=1,NR
                        
                        do jjj=1,NR
                            
                            !-----因节点垂向位移引起的静水压力的变化
                            temHerm1(iii,jjj)=temHerm1(iii,jjj)+ww(j1)*ww(jj)&
                                            &*(nor(1)*eu(iii,1)+nor(2)*eu(iii,2)+nor(3)*eu(iii,3) )*eu(jjj,3)
                            !-----因节点法向量变化引起的静水压力的变化
                            temHerm2(iii,jjj)=temHerm2(iii,jjj)+ww(j1)*ww(jj)*zq*&
                                            &(  ((deu(jjj,2,2)+deu(jjj,3,3) )*eu(iii,1)-deu(jjj,1,2)*eu(iii,2)-deu(jjj,1,3)*eu(iii,3) )*nor(1)&
                                            &+( -deu(jjj,2,1)*eu(iii,1)+(deu(jjj,1,1)+deu(jjj,3,3) )*eu(iii,2)-deu(jjj,2,3)*eu(iii,3) )*nor(2)&
                                            &+( -deu(jjj,3,1)*eu(iii,1)-deu(jjj,3,2)*eu(iii,2)+(deu(jjj,1,1)+deu(jjj,2,2) )*eu(iii,3) )*nor(3)   )
                    
                    
                            temHerm3(iii,jjj)=temHerm3(iii,jjj)+ww(j1)*ww(jj)*zq*&
                                            &( (eu(jjj,1)*deu(iii,1,1)+eu(jjj,2)*deu(iii,1,2)+eu(jjj,3)*deu(iii,1,3))*nor(1)&
                                            &+(eu(jjj,1)*deu(iii,2,1)+eu(jjj,2)*deu(iii,2,2)+eu(jjj,3)*deu(iii,2,3) )*nor(2)&
                                            &+(eu(jjj,1)*deu(iii,3,1)+eu(jjj,2)*deu(iii,3,2)+eu(jjj,3)*deu(iii,3,3) )*nor(3)  )
                            
                    
                        end do
                    
                    end do
                   
                   
               
                    
                end do
            end do
 
       end do  
   end do
   
   
   Initial_Fs(:)=Initial_Fs(:)*rou*g0       !初始静浮力
   
   
   HRM(5,5)=HRM(5,5)+vol*(COB(3)-COG(3))      !此行代码参考basic
   
   HRM(5,3)=HRM(3,5)
   
   HRM=HRM*rou*g0
   
   !-----包涵弹性体模态在内的静水恢复力矩阵
   HERM(:,:)=(temHerm1(:,:)+temHerm2(:,:)+temHerm3(:,:) )*rou*g0
   
   HERM(1:6,1:6)=HRM(1:6,1:6)
   
   !精简矩阵
   do i=1,NR
       do j=1,NR
           if( abs(HERM(i,j))<=1.0  ) then
               HERM(i,j)=0.0  
            end if
        end do    
   end do
   
   
   
    write(11,"(A)") '#   静水恢复力矩阵（刚体）'
    do i=1,6
    do j=1,6
        write(11,"(E15.4,1x,\)",advance='NO') HRM(i,j)
    end do
    write(11,"(/)") 
    end do
   
   
   
   
    write(11,"(A)") '#   静水恢复力矩阵（弹性）'
    do i=1,NR
    do j=1,NR
        write(11,"(E15.4,1x,\)",advance='NO') HERM(i,j)
    end do
    write(11,"(/)") 
    end do
   
 
   
   !至此与结构相关的质量矩阵，刚度矩阵，阻尼矩阵求解完毕
   !静水恢复力矩阵计算完毕
    
    
    !**计算与剖面载荷相关的模态20230425**
    allocate( Eloadr(NBSECT,NR,6 )  )
    Eloadr=0.0;
    
    fdrm1=0.0;   fdrmN=0.0;
    Ns=0.0; Nds=0.0; Ndds=0.0;
    
    do i=7,NR
        !----先付边界条件
        fdrm1(i,:)=( frm(i,2,:)-frm(i,1,:) )/( StruXN(2)-StruXN(1) )
        fdrmN(i,:)=( frm(i,StruSect,:)-frm(i,StruSect-1,:) )/( StruXN(StruSect)-StruXN(StruSect-1) )
        
        do j=1,6
            
            call trispline( StruSect,StruXN(:),frm(i,:,j),fdrm1(i,j),fdrmN(i,j),NBSECT,SecVec(1,:),Ns(:,i,j),Nds(:,i,j),Ndds(:,i,j) )
            
            Eloadr(:,i,j)=Ns(:,i,j)    !----注意单位
            
        end do
    
    end do
    
    !**计算与梁振动运动相关的模态，用于考虑水动升力弹性**
    allocate( MEloadr(SlamNumLine,NR,6) )
    MEloadr=0.0;
    Mfdrm1=0.0;   MfdrmN=0.0;
    MNs=0.0; MNds=0.0; MNdds=0.0;

    do i=7,NR
        !边界条件
        Mfdrm1(i,:)=( drm(i,2,:)-drm(i,1,:) )/( StruXN(2)-StruXN(1) )
        MfdrmN(i,:)=( drm(i,StruSect,:)-drm(i,StruSect-1,:) )/( StruXN(StruSect)-StruXN(StruSect-1) )
        
        do j=1,6
            
            call trispline( StruSect,StruXN(:),drm(i,:,j),Mfdrm1(i,j),MfdrmN(i,j),SlamNumLine,MSecVec(:),MNs(:,i,j),MNds(:,i,j),MNdds(:,i,j) )
            
            MEloadr(:,i,j)=MNs(:,i,j)      !----注意单位
            
        end do
   
    end do
    
    
    
 
   
   !***check***
    write(11,'(A)')adjustl('# 初始重力及静浮力(关于重心)')
    write(11,"(A8,2(e15.4,1x))") 'Fx',Initial_Mg(1),Initial_Fs(1)
    write(11,"(A8,2(e15.4,1x))") 'Fy',Initial_Mg(2),Initial_Fs(2)
    write(11,"(A8,2(e15.4,1x))") 'Fz',Initial_Mg(3),Initial_Fs(3)
    write(11,"(A8,2(e15.4,1x))") 'Mx',Initial_Mg(4),Initial_Fs(4)
    write(11,"(A8,2(e15.4,1x))") 'My',Initial_Mg(5),Initial_Fs(5)
    write(11,"(A8,2(e15.4,1x))") 'Mz',Initial_Mg(6),Initial_Fs(6)
   

   
    close(11)
    
end subroutine Shiphull