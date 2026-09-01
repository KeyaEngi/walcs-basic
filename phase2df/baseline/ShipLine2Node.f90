




subroutine ShipLine2Node(Numline,NumPoint,Point,Tf,Ta)

   use PanelGeometry,only:Nwh,Nl,Node,Nodeb,InstBreakNode,breakKexi
   use Constant,only:Lpp,TrimAng,d_mc
   use ArrayOperations

    implicit none

  integer(4)::Numline
  integer(4),dimension(Numline)::NumPoint
  real(8),dimension(Numline,100,3 )::Point

  real(8)::Tf,Ta     !----首尾吃水

  !******************************全船网格划分所需变量
  real(8),dimension(Numline,100 )::LineArc
  real(8),dimension(Numline,Nwh,3  )::LineEqArcP   !---各原始型线等分弧长点插值
  !-------
  integer(4),dimension(Numline)::LineCase     !----横剖线和首尾纵剖线相交情况
  integer(4)::TemNumIntNode                     
  real(8),dimension(Numline,3)::TemIntNode    !----用于纵向插值的节点
  integer(4),dimension(Numline )::TransNumInt  !----横剖线等分弧长份数
  !------全船网格节点
  real(8),dimension(Nwh,Nl,3)::TemNodeH


  !*****************************************
  !-----湿表面网格划分所需变量
  real(8)::App,Fpp    !----艉/艏吃水高度对应的纵向艉垂线/艏垂线位置

  real(8)::TemTa,TemTf   !---对应与首尾纵剖线水线交点的吃水深度

  real(8),dimension(Numline,101,3)::WetPoint    !----经由静水面截断后的型线
  integer(4),dimension(Numline)::NumwetP

  real(8),dimension(Nwh,Nl,3)::TemNode   !----湿表面网格


  real(8),dimension(Numline)::Length
  integer(4)::temNl

  real(8),dimension(3)::x1,x2,x3
  integer(4)::i,j,k,ii,jj,kk
  real(8)::s,t,s1,s2,t1,t2
  real(8)::k0,b0

  real(8)::Kexi
  !------------------------------
  temNl=Nl/2+1;    !----半船纵向划分节点数

  !------先把最高点对齐
  s=Point(1,NumPoint(1),3)
  do i=2,Numline
    s=max( s,Point(i,NumPoint(i),3) )
  end do

  do i=1,Numline
    Point(i,NumPoint(i),3)=s
  end do


  !******************************************
  !-----先插值全船网格

  !--------先根据Nwh寻找型线上等分弧长点
  LineArc=0.0
  do i=1,Numline

    k=NumPoint(i)
    t=0.0
    do j=2,k
      x1(:)=Point(i,j-1,:);
      x2(:)=Point(i,j,:);

      s=sqrt( (x2(1)-x1(1) )**2.0+(x2(2)-x1(2) )**2.0+(x2(3)-x1(3) )**2.0 )
      t=t+s
      LineArc(i,j)=t    !----各条线累加弦长
    end do
  end do

  !-----开始寻找等分弧长点
  LineEqArcP=0.0;

  do i=1,Numline
    
    k=NumPoint(i)
    !-----等分弧长距离
    s=LineArc(i,NumPoint(i) )/real( Nwh-1 )

    do j=2,Nwh-1
      t=s*real(j-1)   !----插值节点累加弦长

      do jj=2,k
        s1=LineArc(i,jj-1);
        s2=LineArc(i,jj);
        if( (t-s1)*(t-s2)<=0.0 ) then
          Kexi=(t-s1)/(s2-s1);
          exit
        end if
      end do
      !--------------找到等分弧长插值节点
      LineEqArcP(i,j,:)=Point(i,jj-1,:)*(1.0-Kexi)+Point(i,jj,:)*Kexi
    end do
    !----首尾两端点
    LineEqArcP(i,1,:)=Point(i,1,:);
    LineEqArcP(i,Nwh,:)=Point(i,k,:);

  end do

  !******************************************
  !------判别横剖线与首尾纵剖线的相对位置关系
    k=NumPoint(1)
    s1=Point(1,1,1);
    do j=2,k
    s1=max( s1,Point(1,j,1) )   !----艉纵剖线最前面的点
    end do

    k=NumPoint(Numline)
    s2=Point(Numline,1,1);
    do j=2,k
    s2=min( s2,Point(Numline,j,1) )  !---艏纵剖线最后面的点
    end do

  LineCase=0;
  do i=2,Numline-1

    s=Point(i,1,1);
    
    if(s<=s1 ) then    !----和尾部纵剖线由交点
      LineCase(i)=1
    elseif(s>=s2 ) then   !---和艏部纵剖线有交点
      LineCase(i)=2
    else
      LineCase(i)=0     !----无交点，夹在中间
    end if

  end do

  !********************************
  !-----此处重新定义已相交的横剖线的等分弧长数

  TransNumInt=0;

  do i=2,Numline-1
      x3(:)=Point(i,1,:)

      if(LineCase(i)==1  ) then    !----与艉部纵剖线有交点的横剖线的划分
          
          !---------------先判断曲线需要等分弧长的点数
          do j=Nwh,2,-1
            x2(:)=LineEqArcP(1,j,:)
            x1(:)=LineEqArcP(1,j-1,:)

            if( ( (x3(3)-x1(3))*(x3(3)-x2(3))<=0.0  ).and.((x3(1)-x1(1))*(x3(1)-x2(1))<=0.0  ) ) then
              TransNumInt(i)=(Nwh-j+1)
              exit
            end if
          end do

          if( TransNumInt(i)==0 ) then
              do j=Nwh,1,-1
                x2(:)=LineEqArcP(1,j,:)

                if(  x2(3)<=x3(3)  ) then
                  TransNumInt(i)=(Nwh-j)
                  exit
                end if
              end do            
          end if

          !**************************************************************
          !**********开始等分弧长
          LineEqArcP(i,:,:)=0.0;     !----重新定义第i条剖线的等分弧长点

              k=NumPoint(i)     !-----第i条剖线原有节点数
              !-----等分弧长距离
              s=LineArc(i,NumPoint(i) )/real( TransNumInt(i) )    !---根据第i条剖析啊需要划分的份数确定等分弧长           

              do ii=1,TransNumInt(i)-1
                t=s*real(ii)     !----待等分点累积弦长

                do jj=2,k
                  s1=LineArc(i,jj-1);
                  s2=LineArc(i,jj);
                  if( (t-s1)*(t-s2)<=0.0 ) then
                    Kexi=(t-s1)/(s2-s1);
                    exit
                  end if
                end do

                LineEqArcP(i,Nwh-(TransNumInt(i)-ii ),:)=Point(i,jj-1,:)*(1.0-Kexi)+Point(i,jj,:)*Kexi
              end do
              LineEqArcP(i,Nwh,:)=Point(i,k,:)

      elseif( LineCase(i)==2 ) then   !----与艏部纵剖线有交点的横剖线的划分
           
          do j=Nwh,2,-1
            x2(:)=LineEqArcP(Numline,j,:)
            x1(:)=LineEqArcP(Numline,j-1,:)

            if( ( (x3(3)-x1(3))*(x3(3)-x2(3))<=0.0  ).and.((x3(1)-x1(1))*(x3(1)-x2(1))<=0.0  ) ) then
              TransNumInt(i)=(Nwh-j+1)
              exit
            end if
          end do

          if( TransNumInt(i)==0 ) then
              do j=Nwh,1,-1
                x2(:)=LineEqArcP(Numline,j,:)

                if(  x2(3)<=x3(3)  ) then
                  TransNumInt(i)=(Nwh-j)
                  exit
                end if
              end do            
          end if

          !1*************************开始划分
              LineEqArcP(i,:,:)=0.0;     !----重新定义第i条剖线的等分弧长点

              k=NumPoint(i)     !-----第i条剖线原有节点数
              !-----等分弧长距离
              s=LineArc(i,NumPoint(i) )/real( TransNumInt(i) )    !---根据第i条剖析啊需要划分的份数确定等分弧长           

              do ii=1,TransNumInt(i)-1
                t=s*real(ii)     !----待等分点累积弦长

                do jj=2,k
                  s1=LineArc(i,jj-1);
                  s2=LineArc(i,jj);
                  if( (t-s1)*(t-s2)<=0.0 ) then
                    Kexi=(t-s1)/(s2-s1);
                    exit
                  end if
                end do

                LineEqArcP(i,Nwh-(TransNumInt(i)-ii ),:)=Point(i,jj-1,:)*(1.0-Kexi)+Point(i,jj,:)*Kexi
              end do
              LineEqArcP(i,Nwh,:)=Point(i,k,:)

      end if
  end do



  !****************************************
  !------开始插值
  !------------目前为止没有问题
  !-----开始纵向进行等弧长划分
  TemNodeH=0.0

  do j=1,Nwh
    
    !*******************************
    !----先构建待插值节点数组
    TemNumIntNode=0;
    TemIntNode=0.0;

    !-----艉部节点
    TemNumIntNode=TemNumIntNode+1;
    TemIntNode(TemNumIntNode,:)=LineEqArcP(1,j,: );

    do i=2,Numline-1
      
      !-----先判断曲线位置状态
      if( LineCase(i)==0 ) then
        TemNumIntNode=TemNumIntNode+1;

        TemIntNode(TemNumIntNode,: )=LineEqArcP(i,j,: )
      elseif( LineCase(i)==1 ) then
        
        k=TransNumInt(i)     !-----与艉纵剖线有交点的横剖线的等分弧长数
        if( j>=Nwh-(k-1) ) then
            TemNumIntNode=TemNumIntNode+1;

            TemIntNode(TemNumIntNode,: )=LineEqArcP(i,j,:)
        end if

      elseif( LineCase(i)==2 ) then

        k=TransNumInt(i)     !-----与艉纵剖线有交点的横剖线的等分弧长数
        if( j>=Nwh-(k-1) ) then
            TemNumIntNode=TemNumIntNode+1;

            TemIntNode(TemNumIntNode,: )=LineEqArcP(i,j,:)
        end if        
      end if
      
    end do
    !-----艏部剖线节点
    TemNumIntNode=TemNumIntNode+1;
    TemIntNode(TemNumIntNode,:)=LineEqArcP(Numline,j,: );

    
    !--------------先计算纵向累加弦长
    Length=0.0;
    t=0.0;
    do i=2,TemNumIntNode
      x2(:)=TemIntNode(i,:);
      x1(:)=TemIntNode(i-1,:);

      s=sqrt( (x2(1)-x1(1) )**2.0+(x2(2)-x1(2) )**2.0+(x2(3)-x1(3) )**2.0 )
      t=t+s;

      Length(i)=t;
    end do

    !--------------开始等分弧长划分

    s=Length(TemNumIntNode )/real( temNl-1 )
    do i=2,temNl-1
      t=s*real(i-1)   !----插值节点累加弦长

      do jj=2,TemNumIntNode
        s1=Length(jj-1)
        s2=Length(jj)
        if( (t-s1)*(t-s2)<=0.0 ) then
          Kexi=(t-s1)/(s2-s1);
          exit
        end if
      end do
      !--------------找到等分弧长插值节点
      
      TemNodeH(j,i,:)=TemIntNode(jj-1,:)*(1.0-Kexi)+TemIntNode(jj,:)*Kexi
      TemNodeH(j,Nl+2-i,1)=TemNodeH(j,i,1)
      TemNodeH(j,Nl+2-i,2)=-TemNodeH(j,i,2)
      TemNodeH(j,Nl+2-i,3)=TemNodeH(j,i,3)
    end do
    !----首尾两端点

    TemNodeH(j,1,:)=TemIntNode(1,:);
    TemNodeH(j,Nl/2+1,:)=TemIntNode(TemNumIntNode,:);
  end do
  !---------------至此：全船网格划分完毕
  !************************************************

  !------------------先根据首尾吃水寻找真正的首尾垂线位置
  if( Ta<Point(1,1,3) ) then
    write(*,"(A)") '艉吃水小于艉纵剖线最低点位置'
    stop
  elseif( Tf<Point(Numline,1,3) ) then
    write(*,"(A)") '艏吃水小于艏纵剖线最低点位置'
    stop
  end if

  !********************寻找静水面截断的型线
  NumwetP=0;
  WetPoint=0.0;


  kk=0;
  do i=1,NumPoint(1)-1    !----寻找尾部垂线纵向位置
    x1(:)=Point(1,i,:);
    x2(:)=Point(1,i+1,:);

    if( (Ta-x1(3) )*(Ta-x2(3) )<=0.0  ) then
        Kexi=(Ta-x1(3) )/(x2(3)-x1(3) )
        App=x1(1)*(1.0-Kexi)+x2(1)*Kexi
        kk=1;
        exit
    end if
  end do

  if(kk==0 ) then
     App=Point(1,NumPoint(1),1)
  end if

  !--------截断后的艉纵剖线
  NumwetP(1)=i+1;
  WetPoint(1,1:i,:)=Point(1,1:i,:)
  WetPoint(1,i+1,1)=App;
  WetPoint(1,i+1,3)=Ta;


  kk=0;
  do i=1,NumPoint(Numline)-1    !----寻找艏部垂线纵向位置
    x1(:)=Point(Numline,i,:);
    x2(:)=Point(Numline,i+1,:);

    if( (Tf-x1(3) )*(Tf-x2(3) )<=0.0  ) then
        Kexi=(Tf-x1(3) )/(x2(3)-x1(3) )
        Fpp=x1(1)*(1.0-Kexi)+x2(1)*Kexi
        kk=1;
        exit
    end if
  end do

  if(kk==0 ) then
     Fpp=Point(Numline,NumPoint(Numline),1)
  end if

  NumwetP(Numline)=i+1;
  WetPoint(Numline,1:i,:)=Point(Numline,1:i,:)
  WetPoint(Numline,i+1,1)=Fpp;
  WetPoint(Numline,i+1,3)=Tf;

  !*************************
  !--------------------此处重新定义垂线间长，纵倾角度
  Lpp=Fpp-App;

  TrimAng=(Tf-Ta)/Lpp
  d_mc=0.0;
  d_mc(5,1)=TrimAng

  TemTa=Ta
  TemTf=Tf

  !---------------------注意：在后续计算中，涉及到坐标转换的步骤里
  !---------都有一步是减去（X=0用户坐标系下）位置处对应的吃水(Ta)
  !---------因此，在此处将Ta改为X=0位置处的吃水
  Kexi=(0.0-App)/(Fpp-App)
  Ta=TemTa*(1.0-Kexi)+TemTf*Kexi

  !-------静水线斜率
  x1(:)=WetPoint(1,NumwetP(1),: );
  x2(:)=WetPoint(Numline,NumwetP(Numline),: );

  k0=(x2(3)-x1(3) )/(x2(1)-x1(1) )    !----水线所在平面在中纵平面上投影的直线斜率
  b0=x1(3)-k0*x1(1)

  !-----先对曲线进行截断处理，寻找水线交点

  do i=2,Numline-1
      
      k=NumPoint(i);
      kk=0;
      do j=1,k-1
        x1(:)=Point(i,j,:);
        x2(:)=Point(i,j+1,:);

        !-----计算x1,x2竖直对应在水线投影线上的点
        t1=k0*x1(1)+b0
        t2=k0*x2(1)+b0

        if( (t1-x1(3) )*(t2-x2(3)  )<=0.0 ) then
          Kexi=( t1-x1(3) )/(x2(3)-x1(3) )

          NumwetP(i)=j+1;
          WetPoint(i,j+1,:)=x1(:)*(1.0-Kexi)+x2(:)*Kexi;
          kk=1
          exit
        end if

      end do
      WetPoint(i,1:j,:)=Point(i,1:j,:)

      if(kk==0) then
        NumwetP(i)=j+1;

        x2(:)=Point(i,k,:);
        t2=k0*x2(1)+b0;

        WetPoint(i,j+1,1)=x2(1);
        WetPoint(i,j+1,2)=x2(2);
        WetPoint(i,j+1,3)=t2;
      end if

  end do
  !---------------------------------型线截断处理完毕

  !-----开始第一轮等分弧长
  !--------先根据Nwh寻找型线上等分弧长点
  LineArc=0.0
  do i=1,Numline

    k=NumwetP(i)
    t=0.0
    do j=2,k
      x1(:)=WetPoint(i,j-1,:);
      x2(:)=WetPoint(i,j,:);

      s=sqrt( (x2(1)-x1(1) )**2.0+(x2(2)-x1(2) )**2.0+(x2(3)-x1(3) )**2.0 )
      t=t+s
      LineArc(i,j)=t    !----各条线累加弦长
    end do
  end do


  !-----开始寻找等分弧长点
  LineEqArcP=0.0;

  do i=1,Numline
    
    k=NumwetP(i)
    !-----等分弧长距离
    s=LineArc(i,NumwetP(i) )/real( Nwh-1 )

    do j=2,Nwh-1
      t=s*real(j-1)   !----插值节点累加弦长

      do jj=2,k
        s1=LineArc(i,jj-1);
        s2=LineArc(i,jj);
        if( (t-s1)*(t-s2)<=0.0 ) then
          Kexi=(t-s1)/(s2-s1);
          exit
        end if
      end do
      !--------------找到等分弧长插值节点
      LineEqArcP(i,j,:)=WetPoint(i,jj-1,:)*(1.0-Kexi)+WetPoint(i,jj,:)*Kexi
    end do
    !----首尾两端点
    LineEqArcP(i,1,:)=WetPoint(i,1,:);
    LineEqArcP(i,Nwh,:)=WetPoint(i,k,:);

  end do
  !------------------第一轮等分弧长完毕


   !******************************************
   !------判别横剖线与首尾纵剖线的相对位置关系
    k=NumwetP(1)
    s1=WetPoint(1,1,1);
    do j=2,k
      s1=max( s1,WetPoint(1,j,1) )   !----艉纵剖线最前面的点
    end do

    k=NumwetP(Numline)
    s2=WetPoint(Numline,1,1);
    do j=2,k
      s2=min( s2,WetPoint(Numline,j,1) )  !---艏纵剖线最后面的点
    end do

  LineCase=0;
  do i=2,Numline-1

    s=WetPoint(i,1,1);
    
    if(s<=s1 ) then    !----和尾部纵剖线由交点
      LineCase(i)=1
    elseif(s>=s2 ) then   !---和艏部纵剖线有交点
      LineCase(i)=2
    else
      LineCase(i)=0     !----无交点，夹在中间
    end if

  end do

  !-------------------确定与纵剖线相交的横剖线的等分弧长点

  TransNumInt=0;
  do i=2,Numline-1
      x3(:)=WetPoint(i,1,:)

      if(LineCase(i)==1  ) then    !----与艉部纵剖线有交点的横剖线的划分
          

          if( x3(3)>LineEqArcP(1,Nwh,3).or.x3(1)<LineEqArcP(1,Nwh,1)   )  then
              LineEqArcP(i,:,:)=0.0;
              TransNumInt(i)=0;   !----此种情况表明该曲线无法划分网格
          else

              !---------------先判断曲线需要等分弧长的点数
              do j=Nwh,2,-1
                x2(:)=LineEqArcP(1,j,:)
                x1(:)=LineEqArcP(1,j-1,:)

                if( ( (x3(3)-x1(3))*(x3(3)-x2(3))<=0.0  ).and.((x3(1)-x1(1))*(x3(1)-x2(1))<=0.0  ) ) then
                  TransNumInt(i)=(Nwh-j+1)
                  exit
                end if
              end do

              if( TransNumInt(i)==0 ) then
                  do j=Nwh,1,-1
                    x2(:)=LineEqArcP(1,j,:)

                    if(  x2(3)<=x3(3)  ) then
                      TransNumInt(i)=(Nwh-j)
                      exit
                    end if
                  end do            
              end if
              !---------------第i条剖线的等分弧长数确定完毕
                  LineEqArcP(i,:,:)=0.0;     !----重新定义第i条剖线的等分弧长点

                  k=NumwetP(i)     !-----第i条剖线原有节点数
                  !-----等分弧长距离
                  s=LineArc(i,NumwetP(i) )/real( TransNumInt(i) )    !---根据第i条剖析啊需要划分的份数确定等分弧长           

                  do ii=1,TransNumInt(i)-1
                    t=s*real(ii)     !----待等分点累积弦长

                    do jj=2,k
                      s1=LineArc(i,jj-1);
                      s2=LineArc(i,jj);
                      if( (t-s1)*(t-s2)<=0.0 ) then
                        Kexi=(t-s1)/(s2-s1);
                        exit
                      end if
                    end do

                    LineEqArcP(i,Nwh-(TransNumInt(i)-ii ),:)=WetPoint(i,jj-1,:)*(1.0-Kexi)+WetPoint(i,jj,:)*Kexi
                  end do
                  LineEqArcP(i,Nwh,:)=WetPoint(i,k,:)

          end if


      elseif( LineCase(i)==2 ) then   !----与艏部纵剖线有交点的横剖线的划分
           

          if( x3(3)>LineEqArcP(Numline,Nwh,3).or.x3(1)>LineEqArcP(Numline,Nwh,1)   )  then
              LineEqArcP(i,:,:)=0.0;
              TransNumInt(i)=0;   !----此种情况表明该曲线无法划分网格              
          else

              do j=Nwh,2,-1
                x2(:)=LineEqArcP(Numline,j,:)
                x1(:)=LineEqArcP(Numline,j-1,:)

                if( ( (x3(3)-x1(3))*(x3(3)-x2(3))<=0.0  ).and.((x3(1)-x1(1))*(x3(1)-x2(1))<=0.0  ) ) then
                  TransNumInt(i)=(Nwh-j+1)
                  exit
                end if
              end do

              if( TransNumInt(i)==0 ) then
                  do j=Nwh,1,-1
                    x2(:)=LineEqArcP(Numline,j,:)

                    if(  x2(3)<=x3(3)  ) then
                      TransNumInt(i)=(Nwh-j)
                      exit
                    end if
                  end do            
              end if

              !1*************************开始划分
                  LineEqArcP(i,:,:)=0.0;     !----重新定义第i条剖线的等分弧长点

                  k=NumwetP(i)     !-----第i条剖线原有节点数
                  !-----等分弧长距离
                  s=LineArc(i,NumwetP(i) )/real( TransNumInt(i) )    !---根据第i条剖析啊需要划分的份数确定等分弧长           

                  do ii=1,TransNumInt(i)-1
                    t=s*real(ii)     !----待等分点累积弦长

                    do jj=2,k
                      s1=LineArc(i,jj-1);
                      s2=LineArc(i,jj);
                      if( (t-s1)*(t-s2)<=0.0 ) then
                        Kexi=(t-s1)/(s2-s1);
                        exit
                      end if
                    end do

                    LineEqArcP(i,Nwh-(TransNumInt(i)-ii ),:)=WetPoint(i,jj-1,:)*(1.0-Kexi)+WetPoint(i,jj,:)*Kexi
                  end do
                  LineEqArcP(i,Nwh,:)=WetPoint(i,k,:)
          end if

      end if
  end do

  !*************************************
  !----------开始纵向划分网格
  TemNode=0.0;


do j=1,Nwh
    
    !*******************************
    !----先构建待插值节点数组
    TemNumIntNode=0;
    TemIntNode=0.0;

    !-----艉部节点
    TemNumIntNode=TemNumIntNode+1;
    TemIntNode(TemNumIntNode,:)=LineEqArcP(1,j,: );

    do i=2,Numline-1
      
      !-----先判断曲线位置状态
      if( LineCase(i)==0 ) then
        TemNumIntNode=TemNumIntNode+1;

        TemIntNode(TemNumIntNode,: )=LineEqArcP(i,j,: )
      elseif( LineCase(i)==1 ) then
        
        if( TransNumInt(i)==0 ) cycle

        k=TransNumInt(i)     !-----与艉纵剖线有交点的横剖线的等分弧长数
        if( j>=Nwh-(k-1) ) then
            TemNumIntNode=TemNumIntNode+1;

            TemIntNode(TemNumIntNode,: )=LineEqArcP(i,j,:)
        end if

      elseif( LineCase(i)==2 ) then

        if( TransNumInt(i)==0 ) cycle

        k=TransNumInt(i)     !-----与艉纵剖线有交点的横剖线的等分弧长数
        if( j>=Nwh-(k-1) ) then
            TemNumIntNode=TemNumIntNode+1;

            TemIntNode(TemNumIntNode,: )=LineEqArcP(i,j,:)
        end if        
      end if
      
    end do
    !-----艏部剖线节点
    TemNumIntNode=TemNumIntNode+1;
    TemIntNode(TemNumIntNode,:)=LineEqArcP(Numline,j,: );

    
    !--------------先计算纵向累加弦长
    Length=0.0;
    t=0.0;
    do i=2,TemNumIntNode
      x2(:)=TemIntNode(i,:);
      x1(:)=TemIntNode(i-1,:);

      s=sqrt( (x2(1)-x1(1) )**2.0+(x2(2)-x1(2) )**2.0+(x2(3)-x1(3) )**2.0 )
      t=t+s;

      Length(i)=t;
    end do

    !--------------开始等分弧长划分

    s=Length(TemNumIntNode )/real( temNl-1 )
    do i=2,temNl-1
      t=s*real(i-1)   !----插值节点累加弦长

      do jj=2,TemNumIntNode
        s1=Length(jj-1)
        s2=Length(jj)
        if( (t-s1)*(t-s2)<=0.0 ) then
          Kexi=(t-s1)/(s2-s1);
          exit
        end if
      end do
      !--------------找到等分弧长插值节点
      
      TemNode(j,i,:)=TemIntNode(jj-1,:)*(1.0-Kexi)+TemIntNode(jj,:)*Kexi
      TemNode(j,Nl+2-i,1)=TemNode(j,i,1)
      TemNode(j,Nl+2-i,2)=-TemNode(j,i,2)
      TemNode(j,Nl+2-i,3)=TemNode(j,i,3)
    end do
    !----首尾两端点

    TemNode(j,1,:)=TemIntNode(1,:);
    TemNode(j,Nl/2+1,:)=TemIntNode(TemNumIntNode,:);
  end do

  !*****************************************************
  !-----------至此：全船网格与湿表面网格划分完毕
  !*****************************************************

  !-----目前仍在用户坐标系下
  Nodeb(1:Nwh,1:Nl,1:3)=TemNodeH(1:Nwh,1:Nl,1:3)

  Node(1:Nwh,1:Nl,1:3)=TemNode(1:Nwh,1:Nl,1:3)

  !------------------此处给定Nodeb的截取下限
  breakKexi=1.0/2.0
  allocate( InstBreakNode(Nl,3) )

  InstBreakNode(1:Nl,1:3)=Nodeb(2,1:Nl,1:3)
  do j=1,Nl
      x1(:)=Nodeb(1,j,1:3);
      x2(:)=Nodeb(2,j,1:3);
      InstBreakNode( j,1:3 )=x1(1:3)*(1.0-breakKexi)+x2(1:3)*breakKexi
  end do



  !----------------------------------------------------
  !------先在转化到一个临时坐标系，(还差一步平移就可以转换到随船平动坐标系中)
  !----坐标平移(Z方向的平移)
  Node(:,:,3)=Node(:,:,3)-Ta 
  Nodeb(:,:,3)=Nodeb(:,:,3)-Ta

  !----坐标旋转(尚需X方向上的平移)
  do	i=1,Nwh
    do	j=1,NL
        Node(i,j,1:3)= VectorL2G( Node(i,j,1:3),d_mc(4:6,1) ) 

        Nodeb(i,j,1:3)=VectorL2G( Nodeb(i,j,1:3),d_mc(4:6,1) )
    enddo
  end do 


  InstBreakNode(:,3)=InstBreakNode(:,3)-Ta
  do j=1,Nl
      InstBreakNode(j,1:3)=VectorL2G( InstBreakNode(j,1:3),d_mc(4:6,1) )
  end do


  return
end subroutine ShipLine2Node