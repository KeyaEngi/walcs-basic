module Constant
    
    
    implicit none
    
    character(len=300)::InAccess,OutAccess,projname    !路径以及工程名
    
    !计算输入控制参数
    character(len=15)::Non_Linear       !入射波非线性考虑与否
    character(len=15)::Elastomer        !刚体/弹性体
    character(len=15)::Wavectrl         !静水/规则波/（不规则波20231113）
    character(len=15)::Liftctrl         !升力计算方式MLM/BEM
    character(len=15)::Airlift          !气动升力考虑与否（水上飞机/滑行艇）
    !计算输出控制参数
    character(len=15)::PrtMotion        !是否输出运动响应
    character(len=15)::PrtSectload      !是否输出剖面载荷响应结果
    
    character(len=15):: Slamctrl        !---每条型线上的砰击载荷输出与否20231116 'YES'/'NO'
    
    integer(4)::MonDim=6                !运动自由度（刚体）
    integer(4)::NR                      !弹性体模态数目
    
    
    real(8):: rou,depth,U0              !水密度，水深，航速
    
    real(8):: Lpp,Tf,Ta                 !垂线间长，首吃水，尾吃水
    
    real(8):: TrimAng                   !纵倾角
    real(8)::d_mc(6,1)             !该变量用于坐标旋转变换
    
    real(8)::g0=9.81                    !重力加速度
    real(8)::Fn                         !傅汝德数
    
    real(8)::Dtsim                      !时间步长
    integer(4)::WL_SmoothP              !水线光顺周期
    
    character(len=15)::SurgeCtr         !纵荡运动开放状态
    character(len=15)::SwayCtr          !横荡运动开放状态
    character(len=15)::RollCtr          !横摇运动开放状态
    character(len=15)::YawCtr           !首摇运动开放状态
    
    real(8)::pi=4.0*atan(1.0)
    real(8)::DPi=8.0*atan(1.0)
    
    
    integer(4)::IrreCtrl              !静水、规则波、不规则波控制参数
    integer(4)::loadnum                 !----需要计算的工况总数
    integer(4),allocatable,dimension(:)::loadIDname
    real(8),allocatable,dimension(:)::amp_num,ome_num,beta_num,start_num,time_num       !波幅,频率,浪向,计算起始时间,模拟时长
    
    real(8),allocatable,dimension(:)::hs_num,tz_num       !不规则波有义波高/平均跨零周期
    
    real(8)::amp,ome,omee,head,wavek,wavet,wavete,wavel,wllpp     !波幅，频率，遭遇频率，浪向，波数，周期，遭遇周期，波长，波长船长比
    
    integer(4)::Nramp        !平滑函数作用的步数
    integer(4)::Nsimu        !总的迭代步数
    
    !**添加非线性计算方式控制参数**
    integer(4)::nonlinearCtrl          !----控制非线性计算方式（0在原有的型线上进行划分，1在bdf模型上进行划分）

    
    !升沉、纵摇阻尼耗散项20230519
    real(8)::BcoefH,BcoefP           !升沉、纵摇阻尼耗散系数
    
    real(8)::coefT           !弹性控制系数（砰击载荷计及梁剖面速度程度）
    

    
    
    
    
    
end module Constant