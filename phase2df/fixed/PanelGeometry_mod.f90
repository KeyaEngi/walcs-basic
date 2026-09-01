module PanelGeometry
    
    implicit none
    
    !考虑入射波非线性相关变量
    !----壳单元id，节点ID，单元对应属性id
    integer(4)::bdfnum_ele,bdfnum_node
    integer(4),allocatable,dimension(:)::bdfele_kind
    !----单元对应节点id
    integer(4),allocatable,dimension(:,:)::bdfele_node
    !----节点坐标
    real(8),allocatable,dimension(:,:)::bdfnode
    
    real(8),allocatable,dimension(:,:,:)::NodeH,Node    !随船平动坐标系，用户坐标系平均湿表面网格
    !real(8)::Node_num      !平均湿表面网格数目
    real(8),allocatable,dimension(:,:,:)::Nodeb         !随船平动坐标系全船网格
    !real(8)::Nodeb_num     !全船网格数目
    integer(4)::Nwh,NL      !船体径向节点数，船体周向节点数
    
    !-----瞬时湿表面网格截取下限(节点阈值)
    real(8),allocatable,dimension(:,:)::InstBreakNode
    real(8)::breakKexi
    
    integer(4)::Nnode                    !湿表面网格节点
    
    

    
    
    
    
    
    
end module PanelGeometry