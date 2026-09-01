subroutine ShipLine2Node2
    use PanelGeometry
    
    implicit none
    
    integer(4)::i,j,k
    character(10)::tempC     !临时变量
    real(8)::temp            !临时变量
    
    open(13,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.gdf')
    call iutmp(13)
    read(13,*) tempC
    read(13,*) Nodeb_num,temp,temp,temp,temp
    allocate( Nodeb(Nodeb_num,4,3) )
    Nodeb=0.0
    do i=1,Nodeb_num
        read(13,*) ( (Nodeb(i,j,k),k=1,3),j=1,4 )   !用隐do循环读取，需要检查是否有问题
    end do
    
    close(13)
    
    open(13,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.fsa')
    call iutmp(13)
    read(13,*) tempC
    read(13,*) temp     !basic里面的纵倾角计算的准确吗，是否可以直接用
    read(13,*) Node_num,temp
    allocate( NodeH(Node_num,4,3),Node(Node_num,4,3) )
    NodeH=0.0;Node=0.0
    
    do i=1,Node_num
        read(13,*) ( ( NodeH(i,j,k),k=1,3 ),j=1,4 )
    end do
    
    Node(:,:,:)=NodeH(:,:,:)
    close(13)
    
    !****原LT里面这里已经对全船网格以及湿表面网格进行坐标前两步的坐标转换处理****
    !****本程序中放到后面一起处理****
    
    
    
end subroutine ShipLine2Node2