module Slamming
    
    implicit none
    
    integer(4)::SlamLibNumZ                 !型线条数(按高度划分)
    integer(4),allocatable,dimension(:)::SlamLibNumPort,SlamLibNumStar      !---每层水平半剖线节点数
    real(8),allocatable,dimension(:,:,:)::SlamLibPortNode,SlamLibStarNode   !---每层节点坐标
    
    real(8),allocatable,dimension(:,:)::SlamLibType    !---用户坐标系下，关于重心的典型砰击节点坐标
    integer(4)::SlamNumLine                            !---用于计算砰击压力的曲线数目（典型砰击节点数目）
    real(8),allocatable,dimension(:)::SlamAngle        !---平动坐标系下，对应典型砰击节点位置处的剖面倾斜角
    real(8),allocatable,dimension(:)::SlamWidthMcoef   !---对应典型砰击节点位置处的切片宽度修正系数
    
    !剖线几何信息
    integer(4)::SlamIntNumP                            !优化曲线最多节点数（每条砰击曲线（半横剖线）样条插值点总数）
    integer(4),allocatable,dimension(:)::SlamNumLP     !---每条曲线拥有的型值点数
    real(8),allocatable,dimension(:,:,:)::SlamIniType  !----曲线典型节点坐标(用于截取曲线)
    real(8),allocatable,dimension(:,:,:)::SlamNode     !---平动坐标下的砰击节点(左舷)
    real(8),allocatable,dimension(:,:,:)::SlamPoint    !局部坐标系中，每条曲线型值点坐标(y,z),计算压力
    real(8),allocatable,dimension(:,:)::Slamdx         !局部坐标系中，加密剖线斜率
    integer(4),allocatable,dimension(:)::SlamNumAbdent      !----单条砰击剖线中需要抛弃积分的区域个数
    real(8),allocatable,dimension(:,:,:)::SlamAbdentBz      !---放弃区域的上下界(Z,在局部坐标系下)
    
    integer(4),allocatable,dimension(:)::SlamLineCase       !----砰击剖线生成情况
    integer(4),allocatable,dimension(:)::SlamCase           !---砰击事件发生状态
    real(8),allocatable,dimension(:)::SlamIniPenetration    !----初始入侵距离
    real(8),allocatable,dimension(:)::SlamIniRise           !---初始入侵距离对应的液面抬升高度
    integer(4),allocatable,dimension(:)::SlamRelaP          !---时域模拟中，起始时刻的典型节点与波面的相对关系
    
    real(8),allocatable,dimension(:,:,:)::SlamStripTyp        !---切片左右界限典型节点(定位)
    integer(4),allocatable,dimension(:,:)::SlamStripBexist    !---切片左右界限是否划分，0不存在，1存在
    real(8),allocatable,dimension(:,:)::SlamStripBAngle       !---切片左右界限倾角
    real(8),allocatable,dimension(:,:)::SlamWidth             !----切片宽度
    real(8),allocatable,dimension(:,:,:,:)::SlamStripBNode    !----切片左右边界节点坐标
    
    !------MLM前期数据
    integer(4)::SlamNumIntC    !-----剖面轮廓上插值Ct的个数
    integer(4)::SlamNumIntegration  !---计算砰击压力点个数
    
    real(8),allocatable,dimension(:,:)::SlamIntC  !---选定Ct
    real(8),allocatable,dimension(:,:)::SlamNC    !---Ct对应的Ht
    real(8),allocatable,dimension(:,:)::SlamRiseC  !---Ct对应液面抬升W(Ct)
    real(8),allocatable,dimension(:,:)::SlamDerNC   !---Nc对Ct的导数
    
    !========20210830低底升角压力修正数据库
    integer(4)::NumExpSlamCp                !----试验中的压力系数分布点数目
    real(8),allocatable,dimension(:,:)::ExpSlamCp     !---试验中的压力峰值系数
    
    
    !************************砰击网格节点对应的模态位移及导数
    !----砰击节点处的各模态位移及位移导数(供插值用)(左舷)
    real(8),allocatable,dimension(:,:,:,:)::SlamUr
    real(8),allocatable,dimension(:,:,:,:,:)::SlamDur
    !----砰击节点处的各模态位移及位移导数(供插值用)(右舷)
    real(8),allocatable,dimension(:,:,:,:)::SlamRUr
    real(8),allocatable,dimension(:,:,:,:,:)::SlamRDur
    
    !***************20231115输出纵剖线上节点相对位置
    character(len=15)::SlamRelaMctr      !---砰击计算中纵剖线上的点的相对波面位置输出控制
    integer(4)::SlamNumShipRMP,SlamNumWaveRMP   !---用于计算纵剖线上的点个数以及波面上的点个数
    real(8),allocatable,dimension(:,:)::SlamShipRMP,SlamWaveRMP          !---点坐标(随船平动)
    real(8),allocatable,dimension(:,:)::SlamShipIniRMP,SlamWaveIniRMP    !---点坐标(与随船频率坐标系平行，但不减去浮心位置(与静水面持平))
    
    

    
    
    
    
end module Slamming