subroutine Instant_bdf_wetsurface( y,it,t,smtf,InstForce )
    use Constant
    use ShipHullVar
    use PanelGeometry
    use ArrayOperations
    
    !20231113修改
    use IrreWaveVar
    
    implicit none
    
    real(8),dimension(1:2*NR)::y        !-----瞬时摇荡运动
    real(8)::t
    integer(4)::it     !----步数
    real(8)::smtf
    
    real(8),dimension(NR)::InstForce
    
    
    !瞬时位置
    real(8),dimension(1:6)::tem_Dmc,tem_Dmc2
    real(8),dimension(3)::Inst_Cog           !-----瞬时重心位置
    real(8),dimension(bdfnum_node,3  )::bdfnode2      !---瞬时节点位置
    real(8),dimension(bdfnum_node)::Inst_feta    !-----瞬时水线位置
    
    !划分网格专用变量
    real(8),dimension(4,3)::temNode2
    real(8),dimension(4)::temInstfeta
    
    integer(4)::temNumNew_ele   !---重构的单元数量
    integer(4),dimension(2)::temNewKind  !---每个重构单元类型(3或4)
    real(8),dimension(2,4,3)::temNewnode  !---重构单元的节点坐标
    
    
    
    integer(4)::New_bdfnumele   !----重新划分后，单元数量
    integer(4),dimension(bdfnum_ele*2)::New_bdfele_kind   !---新网格类型
    real(8),dimension(bdfnum_ele*2,4,3)::New_bdfnode      !---新网格节点坐标
    
    real(8),dimension(bdfnum_ele*2,4,NR,3)::Newele_Ur   !---新网格中心点各模态下的位移
    
    
    
    !---重新划分后的网格 平动坐标系下，关于重心(初始浮态下)
    real(8),dimension(bdfnum_ele*2,4,3)::New_bdfnode2
    !---重新划分后的网格 平动坐标系下，关于重心(瞬时位置)
    real(8),dimension(bdfnum_ele*2,4,3)::New_bdfnode3
    
    
    !-----模态求解有关变量(临时变量)
    integer(4)::TemNumNewNode
    real(8),allocatable,dimension(:,:)::TemAllNode
    real(8),allocatable,dimension(:,:,:)::TemAllUr
    
    !---与静水压力和入射力有关的量
    real(8),dimension(bdfnum_ele*2,4)::NewInst_feta  !---新节点位置处对应波面起伏
    real(8),dimension(bdfnum_ele*2,4)::NewInst_PIS        !---瞬时节点入射力与静压力的压力
    real(8),dimension(bdfnum_ele*2,4)::NewInst_dtphiI
    
    real(8),dimension(3)::tmpx
    real(8)::tmpt
    
    
    !***输出瞬时湿表面网格***
    character(len=10)::mesh_name
    
    
    
    integer(4)::i,j,k,ii,jj,kk,iii,jjj,kkk
    
    !******以上为变量定义*******
    
    !**刚体运动模态的主坐标**
    tem_Dmc=0.0;
    
    tem_Dmc(1:6)=y(1:6);
    
    tem_Dmc2=tem_Dmc
    tem_Dmc2(1:6)=tem_Dmc2(1:6)+d_mc(1:6,1)
    
    !平动坐标系下重心瞬时位置
    Inst_Cog(1)=Cog(1)+tem_Dmc(1)
    Inst_Cog(2)=Cog(2)+tem_Dmc(2)
    Inst_Cog(3)=Cog(3)+tem_Dmc(3)
    
    !**确定船体节点瞬时位置**
    bdfnode2(1:bdfnum_node,1:3)=bdfnode(1:bdfnum_node,1:3)
    
    bdfnode2(:,1)=bdfnode2(:,1)-Cog(1)
    bdfnode2(:,2)=bdfnode2(:,2)-Cog(2)
    bdfnode2(:,3)=bdfnode2(:,3)-Cog(3)    !----随船平动坐标系下，关于重心
    
    do i=1,bdfnum_node
        bdfnode2(i,1:3)=VectorL2G( bdfnode2(i,1:3),tem_Dmc(4:6) )
    end do
    
    !船体瞬时节点位置(平动坐标系下)
    bdfnode2(:,1)=bdfnode2(:,1)+Cog(1)+tem_Dmc(1)
    bdfnode2(:,2)=bdfnode2(:,2)+Cog(2)+tem_Dmc(2)
    bdfnode2(:,3)=bdfnode2(:,3)+Cog(3)+tem_Dmc(3)
    
    !**根据瞬时节点位置，计算瞬时节点位置对应入射波波面高度**
    Inst_feta=0.0;
    do i=1,bdfnum_node
    !20231113修改
        tmpx(1:3)=bdfnode2(i,1:3)
        
        if( IrreCtrl==0 )then
        tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*t   
        Inst_feta(i)=amp*cos(tmpt )
        else
            Inst_feta(i)=0.0;
            do jj=1,IrreNum
                tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*t + Irrepha(jj)
                Inst_feta(i)=Inst_feta(i)+Irreamp(jj)*cos( tmpt )
            end do            
        end if
   
    end do
    Inst_feta=Inst_feta*smtf      !----瞬时波面起伏
    
    !**开始划分瞬时网格**
    New_bdfnumele=0;
    New_bdfele_kind=0;
    New_bdfnode=0.0;
    TemNumNewNode=0;
    do i=1,bdfnum_ele
        
        !**找到网格的4个节点以及节点瞬时位置对应入射波起伏**
        temNode2=0.0;
        temInstfeta=0.0;
        do j=1,bdfele_kind(i)
            temNode2(j,:)=bdfnode2(bdfele_node(i,j),: ); !瞬时位置

            temInstfeta(j)=Inst_feta(bdfele_node(i,j) ); !瞬时位置水面  
        end do
        
        !**开始截取网格**
        !新的网格数目，每个新网格对应类型，每个网格对应节点坐标
        temNumNew_ele=0;   temNewKind=0;   temNewnode=0.0;
        
        call get_bdf_InstantWater(bdfele_kind(i),temNode2,temInstfeta,temNumNew_ele,temNewKind,temNewnode )
        
        do j=1,temNumNew_ele
            
            New_bdfnumele=New_bdfnumele+1;
            New_bdfele_kind(New_bdfnumele )=temNewKind(j);               !----单元类型
            New_bdfnode(New_bdfnumele,1:4,1:3 )=temNewnode(j,1:4,1:3);   !----瞬时节点位置
            
            TemNumNewNode=TemNumNewNode+temNewKind(j)  !---记录节点数量

        end do

    end do
    
    
    !***输出瞬时湿表面网格，用于检查20230516****
    
    call qinttostr( it-1,mesh_name,10 )
    
    !open(unit=5000,file=trim(adjustl(OutAccess))//'\'//'intmesh'//'\'//trim(adjustl(mesh_name))//'.txt')
    !
    !do i=1,New_bdfnumele
    !    if( New_bdfele_kind(i)==4 ) then
    !        write(5000,"(12(f15.6,1x))") ((New_bdfnode(i,j,k),k=1,3),j=1,4)
    !    elseif( New_bdfele_kind(i)==3 ) then
    !        write(5000,"(12(f15.6,1x))") ((New_bdfnode(i,j,k),k=1,3),j=1,3),(New_bdfnode(i,3,k),k=1,3)
    !    end if   
    !end do
    !close(5000)
    
    
    
    
    !**生成用于计算运动压力积分的网格节点**
    !随船平动坐标系下关于重心(此时单元法向量可以沿用New_bdfele_normal)
    New_bdfnode2=0.0;
    do i=1,New_bdfnumele
        !---单元节点坐标转换
        do j=1,New_bdfele_kind(i)
            tmpx(1:3)=New_bdfnode(i,j,1:3)-Inst_Cog(1:3)
            tmpx(1:3)=VectorG2L(tmpx(1:3),tem_Dmc(4:6) )

            New_bdfnode2(i,j,1:3)=tmpx(1:3)
        end do
    end do
    
    !**生成与瞬时网格同倾斜，但是关于瞬时重心的网格节点**
    New_bdfnode3=0.0;
    do i=1,New_bdfnumele
        !---单元节点坐标转换
        do j=1,New_bdfele_kind(i)
            New_bdfnode3(i,j,1:3)=New_bdfnode(i,j,1:3)-Inst_Cog(1:3)
        end do
    end do
    
    !**计算节点坐标模态**
    allocate( TemAllNode( TemNumNewNode,3 ) )
    allocate( TemAllUr(TemNumNewNode,NR,3)  )
    
    !**确定新划分单元节点处的模态位移(取用刚体模态的,采用节点node3  )**
    !单元节点处刚体模态
    TemAllNode=0.0;
    ii=0;
    do i=1,New_bdfnumele
        !---单元节点坐标转换
        do j=1,New_bdfele_kind(i)
            ii=ii+1;
            TemAllNode(ii,1:3)=New_bdfnode3(i,j,1:3)
        end do
    end do
    
    TemAllUr=0.0;
    call bdfEmesh( TemNumNewNode,TemAllNode(1:TemNumNewNode,1:3 ),TemAllUr(1:TemNumNewNode,1:NR,1:3)  )
    
    !在此返还到关于单元的排列方式
    Newele_Ur=0.0;
    ii=0;
    do i=1,New_bdfnumele
        !---单元节点坐标转换
        do j=1,New_bdfele_kind(i)
            ii=ii+1;
            Newele_Ur(i,j,1:6,1:3)=TemAllUr(ii,1:6,1:3)   !----给定刚体模态
        end do
    end do
    
    !单元节点处弹性体模态
    TemAllNode=0.0;
    ii=0;    
    do i=1,New_bdfnumele
        !---单元节点坐标转换
        do j=1,New_bdfele_kind(i)
            ii=ii+1;
            TemAllNode(ii,1:3)=New_bdfnode2(i,j,1:3)
        end do
    end do
    
    TemAllUr=0.0;
    call bdfEmesh( TemNumNewNode,TemAllNode(1:TemNumNewNode,1:3 ),TemAllUr(1:TemNumNewNode,1:NR,1:3)  )
    ii=0;
    do i=1,New_bdfnumele
        !---单元节点坐标转换
        do j=1,New_bdfele_kind(i)
            ii=ii+1;
        do jj=7,NR
            Newele_Ur(i,j,jj,1:3)=TemAllUr(ii,jj,1:3)  
            !-----将弹性体模态对应节点位移进行旋转变换到瞬时位置(是否需要此步骤?)
            Newele_Ur(i,j,jj,1:3)=VectorL2G( Newele_Ur(i,j,jj,1:3),tem_Dmc(4:6) )
        end do
        end do
    end do
    
    
    !**开始进行非线性入射力计算**
    NewInst_PIS=0.0;
    NewInst_feta=0.0;
    NewInst_dtphiI=0.0;
    
    do i=1,New_bdfnumele
        
        do iii=1,New_bdfele_kind(i)
            
            !先计算每个节点对应的波面起伏(节点瞬时位置处)
            !20231113修改
            tmpx(1:3)=New_bdfnode(i,iii,1:3);
            
            if( IrreCtrl==0 )then
                tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*t   
                NewInst_feta(i,iii)=amp*cos(tmpt )
            else
                NewInst_feta(i,iii)=0.0;
                do jj=1,IrreNum
                    tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*t + Irrepha(jj)
                    NewInst_feta(i,iii)=NewInst_feta(i,iii)+Irreamp(jj)*cos( tmpt )  
                end do
                
            end if
            
            NewInst_feta(i,iii)=NewInst_feta(i,iii)*smtf
            
            !强制限制重新划分的节点不准超过波面起伏高度
            if(tmpx(3)>NewInst_feta(i,iii) ) then
                tmpx(3)=NewInst_feta(i,iii)
            end if
            !20231113修改
            if( tmpx(3)>=0.0 )then      !超过静水面的部分
                if( IrreCtrl==0 )then
                    tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)-omee*t
                    NewInst_dtphiI(i,iii)=amp*g0*cos( tmpt )
                else
                    NewInst_dtphiI(i,iii)=0.0;
                    do jj=1,IrreNum
                        tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)-Irreomee(jj)*t+Irrepha(jj)
                        NewInst_dtphiI(i,iii)=NewInst_dtphiI(i,iii)+Irreamp(jj)*g0*cos( tmpt ) 
                    end do   
                end if
                
                NewInst_dtphiI(i,iii)=NewInst_dtphiI(i,iii)*smtf
                NewInst_dtphiI(i,iii)=NewInst_dtphiI(i,iii)-g0*tmpx(3)
                NewInst_PIS(i,iii)=-NewInst_dtphiI(i,iii)  !---注意正负号
                
            else
                !20231113修改
                if( IrreCtrl==0 )then
                    tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)-omee*t
                    NewInst_dtphiI(i,iii)=amp*g0*exp(wavek*tmpx(3))*cos( tmpt ) 
                else
                    NewInst_dtphiI(i,iii)=0.0;
                    do jj=1,IrreNum
                        tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)-Irreomee(jj)*t+Irrepha(jj)
                        NewInst_dtphiI(i,iii)=NewInst_dtphiI(i,iii)+Irreamp(jj)*g0*exp(Irrek(jj)*tmpx(3))*cos( tmpt )    
                    end do    
                end if
                
                NewInst_dtphiI(i,iii)=NewInst_dtphiI(i,iii)*smtf
                NewInst_dtphiI(i,iii)=NewInst_dtphiI(i,iii)-g0*tmpx(3)

                NewInst_PIS(i,iii)=-NewInst_dtphiI(i,iii)  !---注意正负号

            end if

        end do
    
    end do
    
    InstForce=0.0;
    call Instant_bdf_SurfaceIntegral(New_bdfnumele,New_bdfele_kind(1:New_bdfnumele),NewInst_PIS(1:New_bdfnumele,1:4),&
                 & New_bdfnode3(1:New_bdfnumele,1:4,1:3 ),Newele_Ur(1:New_bdfnumele,1:4,1:NR,1:3),InstForce  )

    
    
    
    
end subroutine Instant_bdf_wetsurface
    
    