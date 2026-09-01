module lift
    
    implicit none
    
    integer(4)::Swing_num,Mwing_num       !机翼划分翼元剖面个数(尾翼、主机翼)
    integer(4),allocatable,dimension(:)::Spoint_num,Mpoint_num            !各翼元剖面上节点数
    
    real(8),allocatable,dimension(:,:,:)::SUpoint,SDpoint,MUpoint,MDpoint              !上、下弧节点坐标
        
    real(8),allocatable,dimension(:,:,:)::Swing_point,Mwing_point                !中弧线节点坐标(随船平动坐标系)
    
    real(8),allocatable,dimension(:,:,:)::Stemp_point,Mtemp_point               !中弧线坐标转换临时变量
    
    real(8),allocatable,dimension(:,:,:)::Swing_node,Mwing_node                  !中弧线节点坐标(局部坐标系)
    
    real(8),allocatable,dimension(:)::SarfaL0,MarfaL0                !各翼元剖面零升力攻角(尾翼、主机翼)
    
    real(8)::arfa               !攻角
        
    real(8),allocatable,dimension(:)::Sca,Mca                         !弦长(尾翼、主机翼)
    
    real(8),allocatable,dimension(:)::air_force                !气动升力
    
    
    
    
    
end module lift