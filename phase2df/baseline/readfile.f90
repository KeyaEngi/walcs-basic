subroutine readfile
    use constant
    use PanelGeometry
    use ShipHullVar
    use Slamming
    use verification
    use ArrayOperations
    
    use lift
    
    implicit none
    
    integer(4)::i,j,k
    
    real(8)::temp                  !临时变量
!与5相关的变量
    character(len=15)::tempC       !临时变量
    real(8)::tempD                 !临时变量
!与6相关的变量
    integer(4)::TemBline                                   !中部型线数目
    integer(4),allocatable,dimension(:)::TemNumPoint          !每条型线上的节点数（包括首尾）
    real(8),allocatable,dimension(:,:,:)::TemLinePoint     !每条型线上的型值点坐标（包括首尾，每条型线上节点不超过100个）
    
    
    
!与8相关的变量
    real(8),dimension(3)::tempNonx    !建模坐标原点
    real(8)::temNonRatio              !建模单位
    character(len=100)::tempNonbdf    !模型bdf名称
!与11有关的变量
    integer(4)::temSline0                                     !粗选砰击节点个数
    real(8),allocatable,dimension(:)::temStyp0x               !粗选砰击节点纵向坐标
    integer(4)::temltdiSNumP                                  !中纵剖线节点个数
    real(8),allocatable,dimension(:,:)::temltdiSP             !中纵剖线节点坐标
    integer(4)::temNumAngarea                                 !砰击倾角对应的划分区域个数
    real(8),allocatable,dimension(:)::temSAng_boundaryx       !砰击倾角对应的划分区域边界
    real(8),allocatable,dimension(:)::temSAng                 !各区域对应的砰击倾角
    integer(4)::temNumSWcof                                   !带宽修正系数对应的划分区域个数
    real(8),allocatable,dimension(:,:)::temSw_boundaryx       !带宽修正系数对应的划分区域边界
    real(8),allocatable,dimension(:)::temSW_coef              !各区域对应的带宽修正系数
    
    real(8)::temAngle                                         !考虑到船体纵倾影响后的砰击角
    integer(4)::temSNump0,temSNump1                           !粗等分弧长节点数30
    real(8),allocatable,dimension(:,:)::temSP0                !粗等分弧长节点坐标
    real(8),allocatable,dimension(:,:)::tem2DSp0,tem2DSp1             !局部坐标系中的粗等分弧长节点坐标
    real(8),dimension(3)::Local_cord,dp1,dp2                          !坐标系转换中间变量
    

    
    
!**********前面是变量定义，以下开始读取文件***********

ww(1)=1.0;ww(2)=1.0      !与高斯积分有关
cor2(1)=-1.0/3.0**0.5;   cor2(2)=1.0/3.0**0.5;

    
 
!1,读取路径信息
    open( 11,file='access.dat' )
    call iutmp(11)
    read(11,'(A)') InAccess
    read(11,'(A)') OutAccess
    close(11)
    
!2,读取工程信息文件
    open( 11,file=trim(adjustl(InAccess))//'\'//'config.hyp' )
    call iutmp(11)
    read(11,*) projname
    close(11)

!3,读取输入、输出控制参数
    open(11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.cfg')
    call iutmp(11)
    !输入
    read(11,'(A)') Non_Linear       !入射波非线性考虑与否
    read(11,'(A)') Elastomer        !刚体/弹性体
    read(11,'(A)') Wavectrl         !静水/规则波
    read(11,'(A)') Liftctrl         !升力计算方式MLM/BEM
    read(11,'(A)') Airlift          !气动升力考虑与否（水上飞机/滑行艇）
    !输出
    read(11,'(A)') PrtMotion        !是否输出运动响应
    read(11,'(A)') PrtSectload      !是否输出剖面载荷响应结果

    !20231116添加
    read(11,'(A)') Slamctrl
    
    close(11)
    
!4，读取环境信息
    open(11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.dhc')
    call iutmp(11)
    read(11,*) rou,depth    !水密度(ton/m^3)/水深
    read(11,*) U0           !航速/kn
    close(11)
    
    rou=rou*1000.0   !---kg/m^3
    U0=U0*0.5144    !航速单位转换为m/s
    
!5,读取垂线间长，首尾吃水
    open(11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.fsa')
    call iutmp(11)
    read(11,"(A)") tempC
    read(11,*) tempD
    read(11,*) k

    do i=1,k
        read(11,*) tempD
    end do

    !       垂线间长          艏吃水/艉吃水
    read(11,*) Lpp,tempD,tempD,Tf,Ta

    close(11)

 !6,读取几何型值文件
    open(11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.mid')
    open(12,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.bow')
    open(13,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.stn')
    call iutmp(11)
    call iutmp(12)
    call iutmp(13)
    !读取中部型值
    read(11,*) TemBline,temp     !---中部型线条数
    
    !----注：每条型线不超过100个点
    allocate( TemNumPoint(TemBline+2),TemLinePoint(TemBline+2,100,3 ) )
    TemNumPoint=0;  TemLinePoint=0.0;
    
    !----中部曲线型值点数
    read(11,*) (TemNumPoint(i+1),i=1,TemBline)
    !----先给定中部型线纵向位置(当前还是站号)
    read(11,*) (TemLinePoint(i+1,1,1),i=1,TemBline)
    
    do i=2,TemBline+1    !----将站号转化为实际坐标值
    j=TemNumPoint(i)   !---中部型线第i条曲线共有j各节点
    
    TemLinePoint(i,1,1)=TemLinePoint(i,1,1)*temp/20.0
    TemLinePoint(i,2:j,1)=TemLinePoint(i,1,1)
    
    end do
    
    do i=1,TemBline
      k=TemNumPoint(i+1)

      read(11,*) (TemLinePoint(i+1,j,2),j=1,k)
      read(11,*) (TemLinePoint(i+1,j,3),j=1,k)
    end do
    
    !读取首部型值
    read(12,*) TemNumPoint(TemBline+2)
    do j=1,TemNumPoint(TemBline+2)
        read(12,*) TemLinePoint(TemBline+2,j,1),TemLinePoint(TemBline+2,j,3)
    end do
    
    !读取尾部型值
    read(13,*) TemNumPoint(1)
    do j=1,TemNumPoint(1)
        read(13,*) TemLinePoint(1,j,1),TemLinePoint(1,j,3)
    end do
    
    TemBline=TemBline+2    !----总的型线条数(从艉到艏依次排列)
    
    close(11)
    close(12)
    close(13)
    
 
    
    !计算初始纵倾角
    TrimAng=(Tf-Ta)/Lpp
    d_mc=0.0;
    d_mc(5,1)=TrimAng
    
    !计算浮汝德数
    Fn=U0/sqrt(g0*Lpp)
    
    
    !6.5，读取升沉、纵摇阻尼耗散项，参照LT横摇阻尼设定20230519
    open(11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.rdm')
    call iutmp(11)
    read(11,*) BcoefH,BcoefP
    
    read(11,*) coefT
    
    close(11)
    
    
    
    
!7,读取计算控制参数文件
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.inp')
    call iutmp(11)
    read(11,*) Dtsim                   !时间步长
    read(11,*) WL_SmoothP              !水线光顺周期
    read(11,*) Nwh,Nl                  !---物面径向网格节点，周向网格节点
    
    read(11,*) SurgeCtr                !纵荡运动开放状态
    read(11,*) SwayCtr                 !横荡运动开放状态
    read(11,*) RollCtr                 !横摇运动开放状态
    read(11,*) YawCtr                  !首摇运动开放状态
    
    close(11)
    !确保水线光顺周期在合理的范围内
    if( WL_SmoothP<1.or.WL_SmoothP>40 )  WL_SmoothP=3
    
!8,若考虑入射波的非线性，则需要对bdf文件进行处理
!****此处与原LT输入格式有所不同****
    if( trim(adjustl(Non_Linear))=='NL' ) then
        open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.nonctrl' )
        call iutmp(11)
        
        read(11,*) nonlinearCtrl        !非线性计算方式控制参数(0在原有的型线上进行划分，1在bdf模型上进行划分)
        
        if( nonlinearCtrl==1 )then
            tempNonx=0.0;
            read(11,*) (tempNonx(i),i=1,3)   !读取建模原点坐标(建模原点在用户坐标系下的坐标)
        
            read(11,*) temNonRatio        !读取建模单位
        
            read(11,"(A)") tempNonbdf       !读取模型bdf名称
        
            
            !****开始读取bdf文件(法向全部调向内),参考LT，需检查****
            call Nonlinear_bdffile_read( tempNonbdf )
        
            bdfnode=bdfnode/temNonRatio     !转换为标准单位
        
            !注意：此时得到的节点尚处于patran坐标系下
            !先转换到用户坐标系下       
            do i=1,bdfnum_node
                bdfnode(i,1)=bdfnode(i,1)+tempNonx(1)
                bdfnode(i,2)=bdfnode(i,2)+tempNonx(2)
                bdfnode(i,3)=bdfnode(i,3)+tempNonx(3)
            end do
            
        end if
        
        close(11)
    
    end if
    
    
    
  

!9,读取质量文件
    open(12,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.dmf')
    call iutmp(12)
    read(12,*) Masssolve
    
    if(trim(adjustl(Masssolve))=='WHOLE')then
        read(12,*) MASS,Cog(1),Cog(2),Cog(3)
        read(12,*)
        read(12,*)
        read(12,*) TotalI11,TotalI22,TotalI33,TotalI13           
    else                !分段质量模型(默然为对称质量模型)
        read(12,*) MASS                           
        read(12,*) Cog(1),Cog(2),Cog(3)    
        read(12,*) NITEM,NBSECT           !质量段数，计算剖面数
        
        allocate(PointCor(1:3,1:NITEM),IX_R(1:NITEM),MM(1:NITEM),x12(1:NITEM,1:2))
        allocate(PointCor2(1:3,1:NITEM))
        PointCor=0.0;IX_R=0.0;MM=0.0;x12=0.0;PointCor2=0.0
        
        allocate(SecVec(1:3,1:NBSECT), SecZSC(1:NBSECT))   !计算载荷的剖面形心以及剪心z坐标
        allocate(SectCog(1:NBSECT,1:3),SectMatrix(1:NBSECT,1:6,1:6),SectMatrix2(1:NBSECT,1:6,1:6) )
        allocate(SectRECog(1:NBSECT,1:3),SectREMatrix(1:NBSECT,1:6,1:6) )
                
        allocate(SecVec2(1:3,1:NBSECT) )
        SecVec=0.0;SecZSC=0.0;SectCog=0.0;SectMatrix=0.0;SectMatrix2=0.0
        SectRECog=0.0;SectREMatrix=0.0;SecVec2=0.0        
        
        do I=1,NITEM
            read(12,*) temp,MM(I),x12(i,1),x12(i,2),PointCor(1,I),IX_R(I),PointCor(3,I) !IX_R横摇惯性半径，PointCor质量段中心坐标
        end do
        do I=1,NBSECT
            read(12,*) temp,SecVec(1,I),SecVec(2,I),SecVec(3,I),SecZSC(i)
        end do 
        
    endif
    close(12)
    
    SecVec2(:,:)=SecVec(:,:)     !用户坐标系下的剖面形心位置
    
    PointCor(2,:)=0.0           !默认分段关于xoz对称
    MM(:)=1000*MM(:)            !质量单位变为kg
    MASS=sum(MM(:))             !总质量
    Cog(1)=sum(MM(:)*PointCor(1,:))/MASS
    Cog(2)=0.0
    Cog(3)=sum(MM(:)*PointCor(3,:))/MASS     !---此时的重心位置还在用户坐标下
    
    Cog0(1:3)=Cog(1:3)   !---Cog0 用户坐标系下的重心坐标
    
    PointCor2(:,:)=PointCor(:,:)      !---PointCor2用户坐标系下，分段重心坐标
    
    !坐标转换(先转换是否有问题?需要验证)
    !---用户坐标系下，分段重心坐标、分段起始纵向坐标相对于整船重心的位置
    PointCor(1,:)=PointCor(1,:)-cog(1)
    PointCor(2,:)=PointCor(2,:)-cog(2)
    PointCor(3,:)=PointCor(3,:)-cog(3) 
    x12(:,1)=x12(:,1)-cog(1)
    x12(:,2)=x12(:,2)-cog(1)
    !---用户坐标系下，剖面形心,剪心位置相对于整船重心位置   
    SecVec(1,:)=SecVec(1,:)-Cog(1)
    SecVec(2,:)=SecVec(2,:)-Cog(2)
    SecVec(3,:)=SecVec(3,:)-Cog(3)
    SecZSC(:)=SecZSC(:)-Cog(3)
    
 
    

    
!10,读取全船网格以及湿表面网格四节点坐标
    !**注意，我们这里定义的Node与LT有所不同，LT里面是整个计算域的网格，这里仅仅是船体湿表面网格
    allocate( NodeH(1:Nwh,1:NL,1:3),Node(1:Nwh,1:NL,1:3) )
    allocate( Nodeb(1:Nwh,1:NL,1:3) )
    
    NodeH=0.0;  Node=0.0;  Nodeb=0.0;
    
    call ShipLine2Node( TemBline,TemNumPoint,TemLinePoint,Tf,Ta )
    
    !**检查船体网格划分**
    open(unit=201,file=trim(adjustl(OutAccess))//'\'//'ShipNode.txt')
        do i=1,Nwh
            do j=1,Nl
                write(201,"(i8,i8,3f16.5)") i,j,(Node(i,j,k),k=1,3)
            end do
        end do
    
    close(201)
    !******20230420水动力网格检查没问题*******
    
    
!11,水动升力数值准备20230407-20230410
    if( trim(adjustl(Liftctrl))=='MLM' )then
    !**读取原始水平型线    
    open(unit=11,file=trim(adjustl(InAccess))//'\'//'SlamP.SLlib')
    call iutmp(11)
    read(11,*) SlamLibNumZ    !型线条数(按高度划分)
    allocate( SlamLibNumPort(SlamLibNumZ),SlamLibNumStar(SlamLibNumZ) )
    !读取每条型线上的节点数（Port左舷，Star右舷，左舷节点+右舷节点=一层节点数）
    read(11,*) (SlamLibNumPort(i),i=1,SlamLibNumZ)
    read(11,*) (SlamLibNumStar(i),i=1,SlamLibNumZ)
    
    !分别寻找左右舷单条型线上的最多节点     
    j=SlamLibNumPort(1)   
    k=SlamLibNumStar(1)
    do i=2,SlamLibNumZ
        j=max(j,SlamLibNumPort(i) );
        k=max(k,SlamLibNumStar(i) );
    end do
    
    !定义储存左右舷节点坐标的变量
    allocate( SlamLibPortNode(SlamLibNumZ,j,3),SlamLibStarNode(SlamLibNumZ,k,3) )
    
    !依次读取每条型线上的节点信息
    do i=1,SlamLibNumZ       !左舷
    do j=1,SlamLibNumPort(i)
        read(11,*) (SlamLibPortNode(i,j,k),k=1,3)
    end do
    end do
    do i=1,SlamLibNumZ    !右舷
    do j=1,SlamLibNumStar(i)
        read(11,*) (SlamLibStarNode(i,j,k),k=1,3)
    end do
    end do
    
    close(11)
    
    !**读取粗选节点纵向位置；中纵剖线节点坐标；砰击倾角相关数据；带宽修正系数相关数据
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'TypePoint.SLtyp')
    call iutmp(11)
    
    !----读取粗选节点数据
    read(11,*) temSline0     !----节点数(粗选节点潜在位置，后续根据实际相交情况筛选)
    allocate( temStyp0x(temSline0)  )     !粗选节点的纵向位置
    !读取粗选节点的纵向位置
    do i=1,temSline0
        read(11,*) temStyp0x(i)
    end do
    
    !----读取中纵剖线数据
    read(11,*) temltdiSNumP       !中纵剖线节点数目
    allocate( temltdiSP(temltdiSNumP,3) )    !中纵剖线节点坐标
    !读取中纵剖线节点坐标
    do i=1,temltdiSNumP
        read(11,*) (temltdiSP(i,j),j=1,3)
    end do
    
    !----读取砰击倾角划分信息
    read(11,*) temNumAngarea    !划分砰击倾角个数
    allocate(  temSAng_boundaryx(temNumAngarea),temSAng(temNumAngarea )  )   !区域边界，区域砰击倾角
    temSAng_boundaryx=0.0; temSAng=0.0;
    !读取倾角区域边界
    do i=1,temNumAngarea-1    
        read(11,*) temSAng_boundaryx(i)
    end do
    !读取各区域砰击倾角大小
    do i=1,temNumAngarea   
        read(11,*) temSAng(i)
    end do
    
    !----切片带宽修正系数(按纵向区域给定不同的修正系数)：该纵向区域与前面划分的砰击区域不一样
    read(11,*) temNumSWcof  !---纵向区域段数
    allocate( temSw_boundaryx( temNumSWcof,2 ),temSW_coef( temNumSWcof ) )     !纵向区域边界（2代表前后边界）；修正系数值
    temSw_boundaryx=0.0;
    temSW_coef=0.0;
    !读取纵向区域边界、每个分段内的修正系数
    do i=1,temNumSWcof   
        read(11,*) temSw_boundaryx(i,1),temSw_boundaryx(i,2),temSW_coef(i)
    end do
    
    close(11)
    
    !LT这里有指定区域砰击载荷计算
    !指定点砰击载荷计算
    !船中纵剖线各点与波面的相对位置计算
    
    

    
    
    !20231115修改，添加船中纵剖线各点与波面的相对位置计算
    SlamRelaMctr='NO'
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'SlamRelaM.slf')
    call IUTMP(11)
    read(11,*) SlamRelaMctr
    if( trim(adjustl(SlamRelaMctr))=='YES' ) then
        read(11,*) SlamNumShipRMP   !----船体纵剖线上的点数
        allocate( SlamShipIniRMP(SlamNumShipRMP,3) )
        SlamShipIniRMP=0.0;
        !读取用户坐标系给定船体纵剖线节点坐标
        do i=1,SlamNumShipRMP
            read(11,*) (SlamShipIniRMP(i,j),j=1,3)
        end do
            
        read(11,*) SlamNumWaveRMP   !---波浪上的点数
        allocate( SlamWaveIniRMP(SlamNumWaveRMP,3) )
        SlamWaveIniRMP=0.0;
        !读取用户坐标系给定波线节点坐标
        do i=1,SlamNumWaveRMP
            read(11,*) (SlamWaveIniRMP(i,j),j=1,3)
        end do
    end if
    close(11)
    
    
    !****这里先不考虑这些模块，后续视需要添加
    
    
    !SlamLibType/SlamNumLine/SlamAngle/SlamWidthMcoef
    !**下面子程序确定典型节点坐标/确定典型节点数量/以及相应区间的砰击倾角/以及相应区间的切片宽度修正系数
    
    !temSline0:粗选潜在砰击节点数目  !temStyp0x:粗选潜在砰击节点纵坐标
    !temltdiSNumP:中纵剖线节点数目   temltdiSP:中纵剖线节点坐标
    !temNumAngarea:砰击倾角对应的划分区域个数
    !temSAng_boundaryx:砰击倾角对应的划分区域边界（这里不包括第一边界的前端，最后一边界的后端）
    !temSAng:不同区域对应的砰击倾角
    !temNumSWcof:带宽修正系数对应的划分区域个数  temSw_boundaryx:带宽修正系数对应的划分区域边界
    !temSW_coef:不同区域对应的带宽修正系数
    
    call Slam_IniTypePoint_get( temSline0,temStyp0x,temltdiSNumP,temltdiSP,&
        & temNumAngarea,temSAng_boundaryx,temSAng,&
        & temNumSWcof,temSw_boundaryx,temSW_coef )
    
    !**以下是与剖面提取以及剖面几何信息有关的全局变量
    SlamIntNumP=500             !优化曲线最多节点数（每条砰击曲线（半横剖线）样条插值点总数）
    
    allocate( SlamNumLP(SlamNumLine)  )          !每条曲线拥有的型值点数
    allocate( SlamIniType(SlamNumLine,3,3) )     !曲线典型节点坐标(用于截取曲线)
    
    allocate( SlamNode(SlamNumLine,SlamIntNumP,3) )     !平动坐标下的砰击节点(左舷)    
    allocate( SlamPoint(SlamNumLine,SlamIntNumP,2) )    !局部坐标系中，每条曲线型值点坐标(y,z),计算压力
    allocate( Slamdx(SlamNumLine,SlamIntNumP) )         !局部坐标系中，加密剖线斜率
    allocate( SlamNumAbdent(SlamNumLine )   )           !----单条砰击剖线中需要抛弃积分的区域个数
    allocate( SlamAbdentBz(SlamNumLine,100,2  )   )     !----每个放弃的区域上下限(局部坐标系下)
    
    allocate( SlamLineCase(SlamNumLine) )        !---砰击剖线生成情况
    allocate( SlamCase(SlamNumLine) )            !---砰击事件发生状态
    allocate( SlamIniPenetration(SlamNumLine) )  !---初始入侵距离
    allocate( SlamIniRise(SlamNumLine) )         !---初始入侵距离对应的压面抬升高度
    allocate( SlamRelaP(SlamNumLine) )           !---时域模拟中，起始时刻的典型节点与波面的相对关系
    
    allocate( MSecVec(SlamNumLine) )             !用户坐标系下关于重心的砰击剖线纵向位置，用于插值计算梁振动模态
    
    SlamNumLP=0;   SlamNode=0.0;   SlamPoint=0.0;
    Slamdx=0.0;    SlamIniType=0.0;
    SlamNumAbdent=0;  SlamAbdentBz=0.0;
    SlamLineCase=0.0;   SlamCase=0.0;   SlamIniPenetration=0.0;
    SlamIniRise=0.0;  SlamRelaP=0.0;
    
    MSecVec=0.0
    
    !**将各位置处的典型砰击节点坐标赋值给各剖线的第一个典型节点
    SlamIniType(:,1,1:3)=SlamLibType(:,1:3)    !典型节点，平动坐标系下（当前还是在用户坐标系下）
    
    !**用户坐标系下关于重心的砰击剖线纵向位置，用于插值计算梁振动模态**
    MSecVec(:)=SlamLibType(:,1)
    
    MSecVec(:)=MSecVec(:)-Cog(1)

    
    
    !**先将节点库中的节点坐标（输入的水平剖线节点坐标）转换到用户坐标系中关于重心上
    do i=1,SlamLibNumZ      !水平剖线层数
    do j=1,SlamLibNumPort(i)       !每层水平半剖线节点数（左舷）
        SlamLibPortNode(i,j,1)=SlamLibPortNode(i,j,1)-cog(1);
        SlamLibPortNode(i,j,2)=SlamLibPortNode(i,j,2)-cog(2);
        SlamLibPortNode(i,j,3)=SlamLibPortNode(i,j,3)-cog(3);
    end do
    end do
    do i=1,SlamLibNumZ    !----右舷
    do j=1,SlamLibNumStar(i)
        SlamLibStarNode(i,j,1)=SlamLibStarNode(i,j,1)-cog(1);
        SlamLibStarNode(i,j,2)=SlamLibStarNode(i,j,2)-cog(2);
        SlamLibStarNode(i,j,3)=SlamLibStarNode(i,j,3)-cog(3);
    end do
    end do

    !将数据库中的典型砰击节点坐标也转换到用户坐标系下关于重心的位置
    do i=1,SlamNumLine
        SlamLibType(i,1)=SlamLibType(i,1)-cog(1);
        SlamLibType(i,2)=SlamLibType(i,2)-cog(2);
        SlamLibType(i,3)=SlamLibType(i,3)-cog(3);
    end do
      
    !**开始计算砰击二维剖线  
    !下面一个do循环的主要思路，先提取粗节点，再将粗节点转到局部坐标系中，再对粗节点进行圆角的过度处理
    !然后进行粗节点的细化处理，随后再将局部局部坐标系中细化节点转移到用户坐标系中（相对于重心）
    
    do i=1,SlamNumLine     !对典型砰击节点进行循环
        !-----注意：纵倾弧度的正负与夹角正负是相反的即：SlamAngle(i)-(-d_mc(5,1) )
        !-----由于夹角与坐标系无关，因此直接将其放到用户坐标系中进行划分
        temAngle=SlamAngle(i)+d_mc(5,1)    !考虑到船体纵倾影响后的砰击角
        !开始提取粗等分弧长节点横剖线,得到用户坐标系下的空间二维剖线节点坐标
        temSNump0=30;   !---等弧长节点数
        allocate( temSP0(temSNump0,3 ) )
        temSP0=0.0;
        
        !SlamLineCase:砰击剖线生成情况
        !temAngle:考虑船体纵倾影响的砰击角
        !temSNump0：半横剖线粗等分弧长节点数（30）
         
        !temSP0:半横剖线粗等分弧长节点坐标（30）
        call Slam_CoarseCurve( i,SlamLineCase(i),temAngle,temSNump0,temSP0 )
        
        if( SlamLineCase(i)==0 ) then    !----没有提取出剖线
            deallocate( temSP0 )
            cycle
        end if 
        
        !开始细化型线
        !先将粗等分弧长节点坐标转换到以典型砰击节点为原点的局部坐标系中
        allocate( tem2DSp0(temSNump0,2 )  )     !---用于细化的局部坐标系节点
        tem2DSp0=0.0;
        
        Local_cord=0.0
        Local_cord(2)=temAngle-Pi/2.0
        
        do j=1,temSNump0     !对半横剖线粗等分弧长节点数（30）进行循环
            dp1(1)=temSP0(j,1)-SlamLibType(i,1)
            dp1(2)=temSP0(j,2)-SlamLibType(i,2)
            dp1(3)=temSP0(j,3)-SlamLibType(i,3)

            dp2(1:3)=VectorL2G( dp1(1:3),Local_cord(1:3) )

            tem2DSp0(j,1)=dp2(2)
            tem2DSp0(j,2)=dp2(3)    !得到局部坐标系下的，半横剖线粗等分弧长剖线节点（以每一个典型砰击节点作为原点）
        end do
        
        !再对提取的点组成的型线进行质量修正,圆角过渡
        allocate( tem2DSp1(temSNump0*10,2 ) )
        tem2DSp1=0.0;
        
        !temSNump0等弧长节点数30
        !tem2DSp0局部坐标系下的原始剖线节点；
          
        !temSNump1粗等分弧长节点进行质量检查并进行圆角过度后的节点数目（大于30）
        !tem2DSp1对粗等分弧长节点进行质量检查并进行圆角过度后的节点坐标
        
        call Slam_ChkCoarseCurve( temSNump0,tem2DSp0,temSNump1,tem2DSp1 )      !暂时未仔细阅读
        
        !最后进行剖线细化
        
        !temSNump1:粗等分弧长节点进行质量检查并进行圆角过度后的节点数目（大于30）
        !tem2DSp1:对粗等分弧长节点进行质量检查并进行圆角过度后的节点坐标
        !SlamIntNumP:每条砰击曲线（半横剖线）样条插值点总数（理论最多）：lt中给定500
          
        !SlamNumLP:每条半横剖线曲线型值点数（最多500）
        !SlamPoint:局部坐标系中，每条半横剖线细等分弧长型值点坐标(y,z),用于计算压力
        !Slamdx:局部坐标系中，加密剖线斜率
        !SlamNumAbdent:每条砰击剖线中需要抛弃积分的区域个数
        !SlamAbdentBz:每个放弃的区域上下限z坐标(局部坐标系下)
        
        !-----细化剖线（此子程序只看了简单情况的处理方式）
        call Slam_CurveOptimize(temSNump1,tem2DSp1(1:temSNump1,1:2),SlamIntNumP,SlamNumLP(i),&
                            & SlamPoint(i,:,:),Slamdx(i,:),SlamNumAbdent(i),SlamAbdentBz(i,:,:)  )
        
        
        !将局部坐标系中的细化节点转换到用户坐标系中（关于重心）
        do j=1,SlamNumLP(i)            
            dp1(1)=0.0;
            dp1(2)=SlamPoint(i,j,1);
            dp1(3)=SlamPoint(i,j,2);

            dp2(1:3)=VectorG2L( dp1(1:3),Local_cord(1:3) )

            SlamNode(i,j,1)=dp2(1)+SlamLibType(i,1)
            SlamNode(i,j,2)=dp2(2)+SlamLibType(i,2)
            SlamNode(i,j,3)=dp2(3)+SlamLibType(i,3)
        end do
        
        
        deallocate( tem2DSp0,temSP0,tem2DSp1 )
        
   
        
    end do
    
    !至此：已获得局部坐标系下细化半横剖线砰击节点坐标SlamPoint,获得局部坐标系中加密剖线节点处斜率
    !获得每条砰击剖线中需要抛弃积分的区域个数并获得局部坐标系中每个放弃的区域上下限z坐标  
    !此外，还获得用户坐标系下关于重心的细化半横剖线砰击节点坐标SlamNode
    
    !**此处加入砰击剖线实际长短限制,用于处理存在剖线交叉的情况
    call Slam_linelength_limit()          !这里起到作用了吗**
    
    !下面给定曲线典型节点，第二个节点和第三个节点坐标
    !do i=1,SlamNumLine     !对砰击剖线进行循环（每条剖线都有三个曲线典型节点）
    !    SlamIniType(i,2,1:3)=SlamNode(i,SlamNumLP(i),1:3)    !第二个典型节点坐标
    !    SlamIniType(i,3,:)=SlamIniType(i,2,:)
    !    SlamIniType(i,3,2)=-SlamIniType(i,3,2)         !第三个典型节点坐标
    !end do
    
    do i=1,SlamNumLine     !对砰击剖线进行循环
          SlamIniType(i,2,1:3)=0.0   !第二个典型节点坐标
          SlamIniType(i,3,:)=0.0
          SlamIniType(i,3,2)=0.0         !第三个典型节点坐标

    end do
    
    
    !20231116输出能够计算砰击载荷的剖线的典型节点位置
    open(unit=7001,file=trim(adjustl(OutAccess))//'\'//'SlamlinePX.txt')
         do i=1,SlamNumLine
             if( SlamLineCase(i)==0 ) cycle
             write(7001,"(3(f15.6,1x))") (SlamIniType(i,1,j),j=1,3)
         end do     
     close(7001)    
    
    
    
    open(unit=7003,file=trim(adjustl(OutAccess))//'\'//'SlamChkLineNum.txt')
    open(unit=7004,file=trim(adjustl(OutAccess))//'\'//'SlamChkLine.txt')
    do i=1,SlamNumLine
        if( SlamLineCase(i)==0 ) cycle
        write(7003,"(i8)") SlamNumLP(i)
    end do
    
    !20240321输出另一半的砰击剖线
    do i=1,SlamNumLine
        if( SlamLineCase(i)==0 ) cycle
        write(7003,"(i8)") SlamNumLP(i)
    end do
    
    do i=1,SlamNumLine
        if( SlamLineCase(i)==0 ) cycle
        do j=1,SlamNumLP(i) 
            write(7004,"(3(f15.6,1x))") (SlamNode(i,j,k),k=1,3)
        end do
    end do
    
    !20240321输出另一半的砰击剖线
    do i=1,SlamNumLine
        if( SlamLineCase(i)==0 ) cycle
        do j=1,SlamNumLP(i) 
            write(7004,"(3(f15.6,1x))") SlamNode(i,j,1),SlamNode(i,j,2)*-1.0,SlamNode(i,j,3)
        end do
    end do    
    
    
    
    
    close(7003)
    close(7004)
     
    

    !****至此：将能够用于计算砰击的剖面的型线全部处理完毕
    !****当前SlamPoint处于局部坐标系下;SlamNode处于用户坐标系下(关于重心的位置)，尚需转换到随船平动坐标系下
    !****由于存在无法提取剖线的情况，因此后续砰击计算相关参数需要有SlamCase来判断是否需要计算
    
 
    elseif( trim(adjustl(Liftctrl))=='BEM' )then
        
    
        
        
    else
        
      
        
        
    endif
    
    
    
!12,根据考虑气动升力与否进行相关变量读取
!程序走通之后再加

    if( trim(adjustl(Airlift))=='YES' )then
    !读取尾翼型值数据    
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.Sair')
    call iutmp(11)
    !读取翼元个数
    read(11,*) Swing_num
    allocate( Spoint_num(Swing_num) )
    Spoint_num=0.0;
    
    !读取各翼元上节点个数
    read(11,*) ( Spoint_num(i),i=1,Swing_num )
    !寻找翼元个数最大值
    j=Spoint_num(1)
    do i=2,Swing_num        
        j=max( j,Spoint_num(i) )
    end do
    
    !读取翼元节点坐标
    allocate( SUpoint(Swing_num,j,3),SDpoint(Swing_num,j,3) )
    SUpoint=0.0;SDpoint=0.0;
    
    allocate( Swing_point(Swing_num,j,3),Stemp_point(Swing_num,j,3),Swing_node(Swing_num,j,2) )     !中弧线上节点坐标（随动、临时、局部坐标系）
    Swing_point=0.0;Stemp_point=0.0;Swing_node=0.0
    
    do i=1,Swing_num          !上弧
        do j=1,Spoint_num(i)            
            read(11,*) ( SUpoint(i,j,k),k=1,3 )
        end do     
    end do
    
    do i=1,Swing_num          !下弧
        do j=1,Spoint_num(i)            
            read(11,*) ( SDpoint(i,j,k),k=1,3 )
        end do     
    end do
    
    close(11)
    
    !读取中部主机翼数据
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.Mair')
    call iutmp(11)
    
    !读取翼元个数
    read(11,*) Mwing_num
    allocate( Mpoint_num(Mwing_num) )
    Mpoint_num=0.0;
    
    !读取各翼元上节点个数    
    read(11,*) ( Mpoint_num(i),i=1,Mwing_num )
    
    !寻找翼元个数最大值
    j=Mpoint_num(1)
    do i=2,Mwing_num
        j=max( j,Mpoint_num(i) )        
    end do
    
    !读取翼元节点坐标
    allocate( MUpoint(Mwing_num,j,3),MDpoint(Mwing_num,j,3) )
    MUpoint=0.0;MDpoint=0.0;
    
    allocate( Mwing_point(Mwing_num,j,3),Mtemp_point(Mwing_num,j,3),Mwing_node(Mwing_num,j,2) )      !中弧线上节点坐标（随动、临时、局部坐标系）
    Mwing_point=0.0;Mtemp_point=0.0;Mwing_node=0.0
    
    do i=1,Mwing_num        !上弧
        do j=1,Mpoint_num(i)
            read(11,*) ( MUpoint(i,j,k),k=1,3 )    
        end do
    end do
    
    do i=1,Mwing_num       !下弧
        do j=1,Mpoint_num(i)
            read(11,*) ( MDpoint(i,j,k),k=1,3 )    
        end do
    end do
    
    close(11)
    
    
    
    
  
        
    end if
    
    

    
    
!13,根据弹性考虑与否，读取模态数据        
    if( trim(adjustl(Elastomer))=='YES' )then
        
    !读取edm文件；包含 弹性模态数目，对应的固有频率，广义质量，结构阻尼系数  
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.emd')
    call iutmp(11)
    
    read(11,*) NR    !----弹性体模态数目
    NR=MonDim+NR     !----总模态数目
    
    !---固有频率，广义质量，阻尼系数
    allocate( Ihome(1:NR),StruM(NR,NR),ceb(1:NR)  )
    allocate( StruB(NR,NR),StruC(NR,NR)  )
    Ihome=0.0; StruM=0.0; ceb=0.0;
    StruB=0.0; StruC=0.0;
    
    do i=7,NR
        read(11,*) Ihome(i),StruM(i,i),ceb(i)
    end do

    close(11)
    
    StruM=StruM*1000.0   !----弹性体模态结构质量(单位转换)
    
    !---读取位移振型文件 .deu
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.deu')
    call iutmp(11)
    read(11,*)
    read(11,*) StruSect   !----结构剖面数
    
    allocate(StruXN(1:StruSect), StruZN(1:StruSect),StruYN(1:StruSect),StruZSC(1:StruSect)   )
    StruXN=0.0;  StruZN=0.0;  StruYN=0.0;  StruZSC=0.0;
    
    read(11,*) (StruZN(i),i=1,StruSect)  !---结构中和轴高度
    read(11,*) (StruXN(i),i=1,StruSect)  !---结构剖面纵向位置
    
    StruZSC=StruZN   !---当前默认剖面剪心高度与中和轴高度一致
    
    !----开始读取位移整型
    allocate(drm(1:NR,1:StruSect,1:6) )
    drm=0.0;
    
    do i=7,NR
    do j=1,StruSect
        read(11,*) temp,(drm(i,j,k),k=1,6)
    end do
    end do
    
    close(11)
    
    !----将剖面位置信息全部转换到用户坐标系下关于重心的相对位置
    StruXN(:)=StruXN(:)-Cog(1)   !----结构剖面x坐标
    StruZN(:)=StruZN(:)-Cog(3)   !----结构剖面水平中和轴高度
    StruZSC(:)=StruZSC(:)-Cog(3) !----结构剖面剪心高度
    StruYN(:)=StruYN(:)-Cog(2)   !----结构纵向中和轴水平位置
    
    !---读取载荷振型文件 .feu
    allocate(frm(1:NR,1:StruSect,1:6) )
    frm=0.0;
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.feu')
    call iutmp(11)
    read(11,*)
    
    do i=7,NR
    do j=1,StruSect
        read(11,*) temp,(frm(i,j,k),k=1,6)   !----注意单位
    end do
    end do
    
    close(11)
    
    

    
    
 
    else
    
    

    end if
    
    
    !**计算节点弹性位移（这里需要计算吗？？是否入射势要考虑网格的弹性变形）
    allocate( Ur(Nwh,NL,Nr,3),dur(Nwh,NL,Nr,3,3)  ) 
    Ur=0.0;   dur=0.0;
    
    allocate(HERM(NR,NR) )
    HERM=0.0;
    
 
    
    call Shiphull
   
    
    Nnode=Nwh*Nl
    
!14,读取计算工况文件
    open(unit=11,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.lci')
    call iutmp(11)
    
    if(trim(adjustl(wavectrl))=='Calm' ) then   !静水
        IrreCtrl=0;
        
        read(11,*)    loadnum   !---工况数
        
        allocate(loadIDname(loadnum) )   !----工况编号
        allocate(amp_num(loadnum),ome_num(loadnum),beta_num(loadnum),start_num(loadnum),time_num(loadnum) )
        
        do i=1,loadnum
            read(11,*) loadIDname(i),time_num(i)
        end do

        close(11)
        
 !20231113修改，目的是后续添加不规则模块
    elseif( trim(adjustl(wavectrl))=='Regular' )then     !规则波
        
        IrreCtrl=0;
        
        read(11,*)    loadnum   !---工况数
        
        allocate(loadIDname(loadnum) )   !----工况编号
        allocate(amp_num(loadnum),ome_num(loadnum),beta_num(loadnum),start_num(loadnum),time_num(loadnum) )
        
        do i=1,loadnum                                        !---周期数
            read(11,*) loadIDname(i),amp_num(i),ome_num(i),beta_num(i),time_num(i)
            !----换浪向(0度为迎浪)
            if(beta_num(i)<=180.0 ) then
                beta_num(i)=180.0-beta_num(i)
            else
                beta_num(i)=beta_num(i)-180.0
            end if
        end do
        
        close(11)
        
    !20231113修改    
    elseif( trim(adjustl(wavectrl))=='Irregular' )then
        
        IrreCtrl=1;
        read(11,*)   loadnum   !---工况数
        
        allocate(loadIDname(loadnum) )   !----工况编号
        allocate( hs_num(loadnum),tz_num(loadnum),beta_num(loadnum),start_num(loadnum),time_num(loadnum)  )     !有义波高，平均跨零周期
        !有义波高，平均跨零周期,浪向，计算开始时间，计算总时间
        do i=1,loadnum
            read(11,*) loadIDname(i),hs_num(i),tz_num(i),beta_num(i),time_num(i)
            !----换浪向(0度为迎浪)
            if( beta_num(i)<=180.0 )then
                beta_num(i)=180.0-beta_num(i)                
            else
                beta_num(i)=beta_num(i)-180.0
            end if  
            
        end do
        
        close(11)

    
        
    end if
    
    
    
    
!15,运行界面显示
    if (trim(adjustl(Non_Linear))=='LT' ) then
        write(*,'(A,I6)')adjustl('程序LT 正在计算中...... ')
    else
        write(*,'(A,I6)')adjustl('程序NL 正在计算中...... ')
    end if
    
    write(*,'(1x,A,I6)')
    write(*,'(3x,A,I6)')adjustl('  计算节点数 : '),Nnode
    write(*,'(1x,A,I6)')
    write(*,'(3x,A,f6.3)')adjustl('  计算航速U0 : '),U0
    if(trim(adjustl(wavectrl))=='Calm' ) then
        write(*,'(10x,A)')adjustl('    静水 ')
 !20231113修改，目的是后续添加不规则模块        
    elseif( trim(adjustl(wavectrl))=='Regular' )then       
        write(*,'(10x,A)')adjustl('    规则波 ')
    elseif( trim(adjustl(wavectrl))=='Irregular' )then
        write(*,'(10x,A)')adjustl('    不规则波 ')  
    end if
    
    
    if( trim(adjustl(Liftctrl))=='MLM' )then
        write(*,"(/)")
        write(*,"(A)") 'MLM方法前期准备数据计算'
        
        !---注意：要求SlamNumIntegration比SlamNumIntC大
        SlamNumIntC=100          !---剖线半宽等分Ct的点数
        SlamNumIntegration=200   !---每次进行砰击力积分时，半个剖面上计算了砰击压力的点数
        
        allocate( SlamIntC(SlamNumLine,SlamNumIntC ) )     !---选定Ct插值点，不包括Ct=0
        allocate( SlamNC(SlamNumLine,SlamNumIntC ) )       !---Ct对应的侵入深度ht,不包括Ct=0
        allocate( SlamRiseC(SlamNumLine,SlamNumIntC ) )    !---Ct对应接触点位置处液面抬升高度,不包括Ct=0
        allocate( SlamDerNC(SlamNumLine,SlamNumIntC+1 ) )  !---D(ht)/D(Ct)，包括Ct=0
        
        do i=1,SlamNumLine
            if( SlamLineCase(i)==0 ) cycle
            write(*,"(A30,i8)") 'MLM前期数据准备,剖线编号',i
            
            k=SlamNumLP(i)
            
            call Slam_WagCondition(k,SlamPoint(i,1:k,1:2),Slamdx(i,1:k),SlamNumIntC,&
                              & SlamIntC(i,:),SlamNC(i,:),SlamDerNC(i,:),SlamRiseC(i,:) )
            
   
        end do
        
        !********20210830此处加上试验测得的砰击压力峰值系数分布
        NumExpSlamCp=14;
        allocate( ExpSlamCp(NumExpSlamCp,2) )
        ExpSlamCp=0.0;
        
        ExpSlamCp(1,1)=0.0;	  ExpSlamCp(1,2)=53.0;
        ExpSlamCp(2,1)=1.0;	  ExpSlamCp(2,2)=100.0;
        ExpSlamCp(3,1)=2.0;	  ExpSlamCp(3,2)=161.0;
        ExpSlamCp(4,1)=2.5;	  ExpSlamCp(4,2)=200.0;
        ExpSlamCp(5,1)=2.8;	  ExpSlamCp(5,2)=210.0;
        ExpSlamCp(6,1)=3.5;	  ExpSlamCp(6,2)=200.0;
        ExpSlamCp(7,1)=4.0;	  ExpSlamCp(7,2)=184.0;
        ExpSlamCp(8,1)=5.0;	  ExpSlamCp(8,2)=159.0;
        ExpSlamCp(9,1)=6.0;    ExpSlamCp(9,2)=132.0;
        ExpSlamCp(10,1)=7.5;  ExpSlamCp(10,2)=100.0;
        ExpSlamCp(11,1)=8.0;   ExpSlamCp(11,2)=92.0;
        ExpSlamCp(12,1)=10.0;  ExpSlamCp(12,2)=69.0;
        ExpSlamCp(13,1)=12.0;  ExpSlamCp(13,2)=53.0;
        ExpSlamCp(14,1)=14.0;  ExpSlamCp(14,2)=37.0;
        
        !-----将度数换为弧度
        do i=1,NumExpSlamCp
            ExpSlamCp(i,1)=ExpSlamCp(i,1)/180.0*Pi
        end do
    
    end if
    

 
    
return    
end subroutine readfile