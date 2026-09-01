module ShipHullVar
    
    implicit none
    
    character(100)::Masssolve         !质量模型整体/分段
    real(8)::MASS                     !总质量
    real(8),dimension(3)::Cog         !随船平动坐标系下重心坐标
    real*8 :: Zg,Xg,Yg                !重心坐标(随船平动坐标系下)
    real(8),dimension(3)::Cog0        !用户坐标系下重心坐标
    real(8)::TotalI11, TotalI22, TotalI33, TotalI13       !整体质量模型对应的惯性半径/惯性矩
    !与分段质量相关变量
    integer(4)::NITEM, NBSECT         !分段模型质量分段数/剖面数
    !-----用户坐标系(关于重心)，分段质心位置；       分段纵向起止位置
    real(8),allocatable,dimension(:,:)::PointCor,PointCor2,x12     
    real(8),allocatable,dimension(:)::IX_R ,MM    !分段横摇惯性半径；分段质量
    !与剖面载荷相关变量
    real(8),allocatable,dimension(:,:)::SecVec     !用户坐标系(关于重心)剖面形心位置
    real(8),allocatable,dimension(:,:)::SecVec2    !---随船平动系，关于重心，分段质心位置
    real(8),allocatable,dimension(:)::SecZSC       !----用户坐标系(关于重心)剪心位置
    
    real(8),allocatable,dimension(:,:)::SectCog,SectRECog  !---分段质心
    real(8),allocatable,dimension(:,:,:)::SectMatrix,SectREMatrix  !---分段质量矩阵(关于重心的)
    real(8),allocatable,dimension(:,:,:)::SectMatrix2   !---关于剖面形心的分段质量矩阵 
    
    !弹性相关变量
    real(8),allocatable,dimension(:)::Ihome,ceb             !---固有频率,阻尼系数
    real(8),allocatable,dimension(:,:)::StruM,StruB,StruC   !---结构质量,阻尼,刚度矩阵
    integer(4)::StruSect                                    !---结构剖面数
    real(8),allocatable,dimension(:)::StruXN,StruZN,StruYN,StruZSC    !---结构剖面中和轴垂向高度，中和轴纵向位置，剖面扭转中心
    real(8),allocatable,dimension(:,:,:)::drm,frm			!梁模型的各个弹性模态在各横剖面的位移
    
    real(8),allocatable,dimension(:,:,:,:)::Ur        !----各模态节点位移
    real(8),allocatable,dimension(:,:,:,:,:)::dUr     !---各模态节点位移的梯度
    
    real(8)::Volx,Voly,Volz,vol                           !从不同方向积分求的体积
    real(8)::cobx,coby,cobz                           !体积分别对x，y，z的一阶静矩
    real(8),dimension(3)::Cob                         !浮心坐标(随船平动坐标系下)
    real(8)::WetArea                                  !---湿表面积
    
    
    real*8,allocatable,dimension(:,:)::Mij            !刚性质量矩阵
    
    real(8),dimension(6)::Initial_Mg                  !----关于重心的重力载荷
    real(8),dimension(6)::Initial_Fs                  !----初始静浮力
    
    real(8),allocatable,dimension(:)::Initial2_Mg
    
    
    
    real(8),allocatable,dimension(:,:,:)::PointCorE   !---分段质心对应剖面的固有振型
    real(8),allocatable,dimension(:)::PointCorYn,PointCorZn,PointCorZsc
    
    real(8),allocatable,dimension(:)::ForMot,ForceI
    
    real*8,dimension(1:6,1:6)::HRM        !刚体运动6自由度静水恢复力矩阵
    real(8),allocatable,dimension(:,:)::HERM     !----包含弹性体模态在内的静水恢复力矩阵
    
    !******计算剖面载荷相关变量20230425******
    real(8),allocatable,dimension(:,:,:)::Eloadr   !---载荷模态
    
    !******计算计算与梁振动运动相关的模态，用于考虑水动升力弹性20230426******
    real(8),allocatable,dimension(:,:,:)::MEloadr     !梁振动模态
    real(8),allocatable,dimension(:)::MSecVec         !用户坐标系下关于重心的砰击剖线纵向位置，用于插值计算梁振动模态
    
    real(8),allocatable,dimension(:,:)::MEsecLoad     !各砰击剖面梁振动
    
    
  
    
    
    
    
end module ShipHullVar