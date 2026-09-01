!输入：NumLine：砰击剖线序号；IniPenetration；初始入侵距离；IniRise：初始入侵距离对应的压面抬升
!ht：相对入侵深度；Vel：相对运动速度；Ace：相对运动加速度
    
!输出：NumPreNode：计算压力节点数量；PreNode：计算压力节点坐标；SlamPre：计算节点压力
    
!SlamIntNumP：每条原始曲线样条插值点总数(理论最多，程序中给定500)；SlamNumLP：每条曲线拥有的型值点数
!SlamPoint：局部坐标系剖线节点坐标；Slamdx：局部坐标系，加密剖线斜率；    
!SlamNumIntC：：剖线半宽等分ct点数；SlamNumIntegration：每次进行砰击力积分时，半个剖面上计算了砰击压力的点数    
!SlamIntC：选定Ct插值点，不包括Ct=0;SlamNC:Ct对应的侵入深度ht,不包括Ct=0;SlamRiseC:Ct对应接触点位置处液面抬升高度,不包括Ct=0
!SlamDerNC:D(ht)/D(Ct)，包括Ct=0;SlamNumAbdent:每条砰击剖线需要抛弃积分的区域个数
!SlamAbdentBz:抛弃积分区域上下限；
!NumExpSlamCp：试验测得压力峰值系数点数；ExpSlamCp：压力峰值系数
!其余为给定点砰击压力计算
    
    

!-----剖面轮廓砰击压力计算
subroutine Slam_Pressure(NumLine,IniPenetration,IniRise,ht,Vel,Ace,NumPreNode,PreNode,SlamPre)
  
  use Slamming,only:SlamIntNumP,SlamNumLP,SlamPoint,Slamdx,SlamNumIntC,SlamNumIntegration,&
                      & SlamIntC,SlamNC,SlamRiseC,SlamDerNC,SlamNumAbdent,SlamAbdentBz,NumExpSlamCp,ExpSlamCp

  implicit none
  integer(4)::NumLine
  real(8)::IniPenetration    !---初始入侵距离
  real(8)::IniRise           !---前期数据中，初始入侵距离下对应的液面抬升高度
  real(8)::ht                !---瞬时入侵距离
  real(8)::Vel               !---瞬时速度
  real(8)::Ace               !---瞬时加速度
  
  !---------------------------剖面几何信息
  integer(4)::Num0
  real(8),dimension(SlamIntNumP,2)::Point0     !局部坐标系中，砰击节点坐标
  real(8),dimension(SlamIntNumP)::Dx           !局部坐标系中，加密剖线斜率

  !---------------------------计算砰击前期计算参数
  integer(4)::temNumIntC          !剖线半宽等分ct点数
  integer(4)::temNumIntegration    !每次进行砰击力积分时，半个剖面上计算了砰击压力的点数
  real(8),dimension(SlamNumIntC )::temIntC    !选定Ct插值点，不包括Ct=0
  real(8),dimension(SlamNumIntC )::temNC      !Ct对应的侵入深度ht,不包括Ct=0
  real(8),dimension(SlamNumIntC+1 )::temDerNc    !D(ht)/D(Ct)，包括Ct=0
  real(8),dimension(SlamNumIntC )::temRise      !Ct对应接触点位置处液面抬升高度,不包括Ct=0

  !------湿表面压力积分点相关变量
  real(8),dimension(SlamNumIntegration )::temxi   !---压力积分点X坐标
  real(8),dimension(SlamNumIntegration )::temDxi  !---压力积分点处对应斜率
  real(8),dimension(SlamNumIntegration )::temfx   !---压力积分点Y坐标

  real(8)::temDecCtNc
  real(8),dimension(SlamNumIntegration )::A1,A2,A3,A4

  integer(4)::NumPreNode   !----计算压力的节点数量
  real(8),dimension(0:SlamNumIntegration,1:2 )::PreNode  !---压力节点坐标
  real(8),dimension(0:SlamNumIntegration )::SlamPre   !---压力

  real(8)::TRise    !---t时刻前期数据对应抬升高度(理论)
  real(8)::RHt      !---实际入侵距离
  real(8)::Ct       !---实际对应的接触点位置
  real(8)::Force,Force2

  !---------------数据输出变量
  integer(4)::SpreCtr
  integer(4)::temNumAbdent
  !-------给定点砰击压力
  real(8)::temSPintPre  
  
  !---------------试验修正变量
  real(8)::temExpCtBeta     !----接触点底升角
  real(8)::temExpCp         !----对应底升角对应的砰击压力系数
  real(8)::temExpPreMax     !----砰击压力峰值
  character(len=15)::ExpCpModyCtr
  real(8)::temSimPreMax     !----压力积分点组中最大的砰击压力
  
  real(8)::temCtPreRatio
  real(8)::temModyB1    !----修正的下边界
  real(8),dimension( SlamNumIntegration )::temExpModyRatio
  

  real(8),allocatable,dimension(:,:)::temp1
  integer(4)::i,j,k,ii,jj,kk
  real(8)::s,t,s1,s2,t1,t2,s3,t3,tt1,tt2
  real(8)::ta1,ta2,ta3,ta4
  real(8)::tsfx,tsdx
  real(8)::Kexi
  real(8),dimension(2)::tempx
  
  
  
  !**************************************
  temExpModyRatio=1.0;

  
  !****************先选定待求曲线(几何信息)
  Point0=0.0;  Dx=0.0;
  
  Num0=SlamNumLP(NumLine )   !轮廓线包含几何点数
  Point0(1:Num0,1:2)=SlamPoint(NumLine,1:Num0,1:2 )    !局部坐标系中，剖线上砰击节点坐标
  Dx(1:Num0)=Slamdx(NumLine,1:Num0)      !局部坐标系中，细化剖线斜率

  s=Point0(1,1);    !转换坐标系（局部坐标系坐标原点移到第一个砰击节点）
  t=Point0(1,2);
  do i=1,Num0
    Point0(i,1)=Point0(i,1)-s
    Point0(i,2)=Point0(i,2)-t
  end do

  !*******************************压力计算前期储备数据
  !-----CT与ht的关系曲线
  temIntC=0.0;   temNC=0.0;  temDerNc=0.0;  temRise=0.0;
  temNumIntC=SlamNumIntC                       !----半剖线上插值Ct节点数
  temNumIntegration=SlamNumIntegration         !----半剖线砰击压力计算节点数
  temIntC(1:temNumIntC)=SlamIntC(NumLine,1:temNumIntC )
  temNC(1:temNumIntC)=SlamNC(NumLine,1:temNumIntC )
  temDerNc(1:temNumIntC+1)=SlamDerNC(NumLine,1:temNumIntC+1 )
  temRise(1:temNumIntC)=SlamRiseC(NumLine,1:temNumIntC )   !---实时液面抬升插值节点

  !------湿表面压力积分点位置xi
 
!!  write(*,"(A20,f12.3)") '浸入深度',ht
  !***********************************
  !***开始计算压力***
  SpreCtr=0

  if( ht<=0.01 ) then
     Force=0.0;
     Rht=ht;            !ht相对入侵深度；Rht实际入侵深度
     SpreCtr=0;
  else

  A1=0.0; A2=0.0; A3=0.0; A4=0.0;
  do
    
    !-----先确定真实的用于计算的入侵深度以及对应的接触点Ct
    if( IniRise==0.0 ) then       !初始入侵深度对应的压面抬升
        
        !----先根据ht插Ct(此时出现的情况是剖线完全浸没在水中)
        !----直接将载荷赋0是不是比较草率??
        if( ht>temNC(temNumIntC)  ) then
          Force=0.0;
          Rht=ht;
          SpreCtr=0;
          exit
        end if

        do i=1,temNumIntC      !对版剖线上的ct点数进行循环
          if(i==1) then
            s1=0.0;
            t1=0.0;
            tt1=0.0;
          else
            s1=temNC(i-1);      !----入侵深度
            t1=temIntC(i-1);    !----接触点x
            tt1=temRise(i-1);   !----液面抬升
          end if
          s2=temNC(i);
          t2=temIntC(i);
          tt2=temRise(i);

          if( (ht-s1)*(ht-s2)<=0.0 ) then
            Kexi=(ht-s1)/(s2-s1)
            
            Ct=t1*(1.0-Kexi)+t2*Kexi         !插值得到接触点位置Ct           
            TRise=tt1*(1.0-Kexi)+tt2*Kexi    !插值得到接触点位置处液面抬升高度
            exit
          end if
        end do

        Rht=ht;
    else
        
        do i=1,temNumIntC
          if(i==1) then
            s1=0.0;
            tt1=0.0;
          else
            s1=temNC(i-1);      !----入侵深度
            tt1=temRise(i-1);   !----液面抬升
          end if
          s2=temNC(i);
          tt2=temRise(i);

          if( ht>temNC(temNumIntC ) ) then     !如果ht大于最大的入侵深度，则取最大的入侵深度
            TRise=temRise(temNumIntC )
            exit
          end if

          if( (ht-s1)*(ht-s2)<=0.0 ) then
            Kexi=(ht-s1)/(s2-s1)
            
            TRise=tt1*(1.0-Kexi)+tt2*Kexi      !插值得到接触点位置处液面抬升高度
            exit
          end if
        end do        
        
        s=ht+TRise-IniRise     !----真实的，由顶点到接触点的高度(此时已经去除了初始入侵的影响)
        if(s>Point0(Num0,2)  ) then
          Force=0.0;
          Rht=ht;
          Ct=Point0(Num0,1);
          SpreCtr=0;
          exit
        end if

        !------再根据s插值Ct,然后计算ht
        !------对轮廓线，由y插x
        call SlamPoint_Interpolation(Num0,Point0(:,2),Point0(:,1),s,Ct   )
       
        !------根据Ct插Rht
        if( Ct>temIntC(temNumIntC) ) Ct=temIntC(temNumIntC)
        
        do i=1,temNumIntC
          if(i==1) then
            s1=0.0;
            t1=0.0;
          else
            t1=temNC(i-1);      !----入侵深度
            s1=temIntC(i-1);    !----接触点x
          end if
          t2=temNC(i);
          s2=temIntC(i);
          if( (Ct-s1)*(Ct-s2)<=0.0 ) then
            Kexi=(Ct-s1)/(s2-s1)
            
            Rht=t1*(1.0-Kexi)+t2*Kexi    !插值得到接触点位置Ct
            exit
          end if
        end do

    end if

    
    !-----开始计算接触区域内压力
    !----插值得到DerNc
    allocate( temp1(temNumIntC+1,2 ) )
    temp1=0.0;
    temp1(2:temNumIntC+1,1)=temIntC(1:temNumIntC)
    temp1(:,2)=temDerNc(:)
    call SlamPoint_Interpolation(temNumIntC+1,temp1(:,1),temp1(:,2),Ct,temDecCtNc   )

    deallocate(temp1 )

    s=Ct/real(temNumIntegration )
    do i=1,temNumIntegration        !-----计算压力节点在剖线上的坐标/斜率
      temxi(i)=s*real(i)
      if(i==temNumIntegration ) temxi(i)=Ct
      !----插值得到fx
      call SlamPoint_Interpolation(Num0,Point0(:,1),Point0(:,2),temxi(i),temfx(i)   )
      !----插值得到dx
      call SlamPoint_Interpolation(Num0,Point0(:,1),Dx(:),temxi(i),temDxi(i)   )
      
    end do
    
    !=====20210830
    !==========先根据Ct对应的底升角确定试验中的Cp
    temExpCtBeta= atan( temDxi(temNumIntegration) )   !----底升角(弧度)
    
    if( temExpCtBeta>ExpSlamCp(NumExpSlamCp,1) ) then
        
        ExpCpModyCtr='NO'
    else
        
        ExpCpModyCtr='YES'         !----需要修正
        !----先找出Ct处对应的实验砰击压力峰值系数
        call SlamPoint_Interpolation( NumExpSlamCp,ExpSlamCp(:,1),ExpSlamCp(:,2),temExpCtBeta,temExpCp   )
    
        temExpPreMax=temExpCp/2.0*Vel**2.0    !----不包括水密度(按试验测试砰击压力峰值系数计算出的砰击压力)
    end if
    
    
    do i=1,temNumIntegration-1
      
      A1(i)=1.0/temDecCtNc*Ct/sqrt( Ct**2.0-temxi(i)**2.0 )

      A2(i)=0.5*Ct**2.0/(Ct**2.0-temxi(i)**2.0 )/(1.0+temDxi(i)**2.0)

      A3(i)=0.5*temDxi(i)**2.0/(1.0+temDxi(i)**2.0)

      A4(i)=sqrt( Ct**2.0-temxi(i)**2.0 )+temfx(i)-Rht
    end do

    !------计算节点压力
    !---CT处的压力区域负无穷
    SlamPre=0.0;  PreNode=0.0;
    do i=1,temNumIntegration-1
      PreNode(i,1)=temxi(i);
      PreNode(i,2)=temfx(i);
      
      s1=A1(i)-A2(i)-A3(i)   !---与速度有关的量
      s2=A4(i)               !---与加速度有关的量
      !-----计算点压力
      SlamPre(i)=s1*Vel**2.0+s2*ace   
    end do
     
    
    !-----------此处确定计算点的压力最大值
    if( ExpCpModyCtr=='YES'   ) then
    
        temSimPreMax=maxval( SlamPre(1:temNumIntegration-1) )
        
        !-----计算接触点出的比例系数
        if( temSimPreMax<=0.0 ) then
            !-----20210913当前将系数强制令为0
            temExpModyRatio( 1:temNumIntegration )=1.0;
            
        else
            temCtPreRatio=temExpPreMax/temSimPreMax
            
            if( temCtPreRatio>=1.0 ) then
                !-----20210913当前将系数强制令为0
                temExpModyRatio( 1:temNumIntegration )=1.0;                
                
            else
                temModyB1=0.0;
                
                !-----确定对应压力计算点出的修正系数
                do i=1,temNumIntegration
            
                    s=temxi(i)/Ct    !----积分点比例系数
                    !----确定每个压力积分点对应的
                    if( s>temModyB1 ) then
                        t=(1.0/temCtPreRatio)**( 1.0/(1.0-temModyB1) )
            
                        temExpModyRatio(i)=t**(-( s-temModyB1 )  )
                        
                        if( temExpModyRatio(i)>1.0 ) temExpModyRatio(i)=1.0
                        
                    else
                
                        temExpModyRatio(i)=1.0;
                    end if
                end do
             end if
        end if
        
    end if
    
    !-------计算temxi=0 处的压力
    s=0.0;

    tsfx=0.0;
    tsdx=temDxi(1);

    ta1=1.0/temDecCtNc*Ct/sqrt( Ct**2.0-s**2.0 )
    ta2=0.5*Ct**2.0/(Ct**2.0-s**2.0 )/(1.0+tsdx**2.0)
    ta3=0.5*tsdx**2.0/(1.0+tsdx**2.0)
    ta4=sqrt( Ct**2.0-s**2.0 )+tsfx-Rht

    t=(ta1-ta2-ta3)*Vel**2.0+ta4*ace   !---计算点压力(temxi=0.0)
    SlamPre(0)=t

    !****************************
    !------确定临界点
    k=-1
    do i=temNumIntegration-1,0,-1
      if( SlamPre(i)>=0.0 ) then
        k=i
        exit
      end if
    end do

    if(k==-1 ) then   !----全是负压
      Force=0.0;
      SpreCtr=0;
      exit
    end if

    !-----砰击载荷是否为0的判别依据
    !-----到此一步，就说明有不为0的砰击载荷
    SpreCtr=1     
    !-----二分法确定压力为0的节点位置
    if(k==0) then
      s1=0.0;  s2=temxi(k+1);
    else
      s1=temxi(k);  s2=temxi(k+1);   !---初始上下限
    end if
    j=0
    do
      
      j=j+1
      s=(s1+s2)/2.0   !---temxi

      !-----插值fx
      call SlamPoint_Interpolation(Num0,Point0(:,1),Point0(:,2),s,tsfx  )
      !-----插值dx
      call SlamPoint_Interpolation(Num0,Point0(:,1),Dx(:),s,tsdx   )

      ta1=1.0/temDecCtNc*Ct/sqrt( Ct**2.0-s**2.0 )
      ta2=0.5*Ct**2.0/(Ct**2.0-s**2.0 )/(1.0+tsdx**2.0)
      ta3=0.5*tsdx**2.0/(1.0+tsdx**2.0)
      ta4=sqrt( Ct**2.0-s**2.0 )+tsfx-Rht

      t=(ta1-ta2-ta3)*Vel**2.0+ta4*ace   !---计算点压力

      if( abs(t)<=0.001 ) then
        NumPreNode=k+1;

        PreNode(NumPreNode,1)=s
        PreNode(NumPreNode,2)=tsfx

        SlamPre(NumPreNode)=0.0
        exit
      elseif( j>=50 ) then
        NumPreNode=k+1;

        PreNode(NumPreNode,1)=s
        PreNode(NumPreNode,2)=tsfx

        SlamPre(NumPreNode)=0.0
        exit
      end if

      if( t<0.0 ) then
        s2=s
      elseif( t>0.0 ) then
        s1=s
      end if
    end do

    !------------输出给定节点砰击压力
    !-----强制去除压力为负的节点
    do i=0,NumPreNode
      if(SlamPre(i)<0.0  ) SlamPre(i)=0.0
    end do

    
    !********************此处进行试验修正
    if( ExpCpModyCtr=='YES'   ) then  
        
        SlamPre(0)=SlamPre(0)*temExpModyRatio(1)
        
        do i=1,NumPreNode
            SlamPre(i)=SlamPre(i)*temExpModyRatio(i)
        end do
    end if
    
 
    exit

  end do

  end if

  !-----输出压力
  !-----计算的载荷为0的情况
    if( SpreCtr==0 ) then
      NumPreNode=1
 
      PreNode=0.0;
      PreNode(1,:)=Point0(2,:)

      SlamPre=0.0;
    end if


    !--------------20210421处理凹陷区域的节点压力
    temNumAbdent=SlamNumAbdent( NumLine )

    if( temNumAbdent>0 ) then
        
        do i=0,NumPreNode
            s=PreNode(i,2)   !---局部坐标系下，压力计算点的垂向坐标

            kk=0;
            do j=1,temNumAbdent
                s1=SlamAbdentBz(NumLine,j,1);    !---下界
                s2=SlamAbdentBz(NumLine,j,2);    !---上界

                if( ( s-s1 )*( s-s2 )<0.0 ) then
                    kk=1;
                    exit
                end if
            end do
            
            if( kk==1 ) then      !-----将属于凹陷区域的压力强制赋0
                SlamPre(i)=0.0;
            end if
        end do
    end if 



!!    !------积分
!!
!!    s=0.0;
!!
!!    do i=1,NumPreNode
!!      s1=PreNode(i-1,1)
!!      s2=PreNode(i,1)
!!
!!      t1=SlamPre(i-1)
!!      t2=SlamPre(i)
!!
!!
!!      s=s+(t1+t2)*(s2-s1)/2.0
!!    end do
!!
!!    Force2=s*2.0




  return
end subroutine Slam_Pressure
