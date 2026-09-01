subroutine Instant_wetsurface( y,it,t,smtf,ForceIS )
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
    
    real(8),dimension(NR)::ForceIS     !----入射力矩静浮力的合力
    
    real(8),dimension(Nwh,Nl)::Inst_PIS   !---瞬时节点入射力与静压力的压力
    real(8),dimension(Nwh,Nl,3)::Node1H     !-----划分后的用于计算载荷的网格(用户坐标系下，关于重心)
    
    real(8),dimension(1:Nwh,1:Nl)::Inst_feta    !-----瞬时水线位置
    real(8),dimension(Nl,3)::breaknode0         !----节点阈值
    real(8),dimension(Nl)::breakInst_feta
    
    
    
    real(8),dimension(1:6)::tem_Dmc,tem_Dmc2
    real(8),dimension(Nwh,Nl,3)::Node0
    real(8),dimension(Nwh,Nl,3)::Node1      !-----划分后的瞬时湿表面网格
    real(8),dimension(Nwh,Nl,3)::Node2      !-----瞬时湿表面网格，关于重心，平动坐标系下
    real(8),dimension(Nwh,Nl,3)::Node3    
    integer(4)::submerge                  !整船淹湿状况
    
    real(8),dimension(Nwh,Nl,NR,3)::InstUr,InstUr2
    real(8),dimension(Nwh,Nl,NR,3,3)::InstDur

    real(8),dimension(3)::Inst_Cog         !-----瞬时重心位置
    
    real(8),dimension(Nwh,Nl)::Inst_dtphiI

    real(8),dimension(3)::vec1,vec2
    real(8),dimension(3)::tmpx
    real(8)::tmpt
    
    
    integer(4)::i,j,k,ii,jj,kk
    !*******以上是变量定义********
    
    !**刚体运动模态的主坐标**
    tem_Dmc=0.0
    tem_Dmc(1:6)=y(1:6)
    
    tem_Dmc2=tem_Dmc
    tem_Dmc2(1:6)=tem_Dmc2(1:6)+d_mc(1:6,1)    !考虑船体的初始倾角
    
    !**随船平动坐标系下重心瞬时位置**
    Inst_Cog(1)=Cog(1)+tem_Dmc(1)
    Inst_Cog(2)=Cog(2)+tem_Dmc(2)
    Inst_Cog(3)=Cog(3)+tem_Dmc(3)
    
    !**船体瞬时节点坐标**
    Node0=Nodeb          !随船平动坐标系下的全船网格
    
    Node0(:,:,1)=Node0(:,:,1)-Cog(1)
    Node0(:,:,2)=Node0(:,:,2)-Cog(2)
    Node0(:,:,3)=Node0(:,:,3)-Cog(3)    !----随船平动坐标系下，关于重心的全船网格
    
    do i=1,Nwh
    do j=1,Nl     
        Node0(i,j,1:3)=VectorL2G(Node0(i,j,1:3),tem_Dmc(4:6) )
    end do
    end do
    
    !船体瞬时节点位置
    Node0(:,:,1)=Node0(:,:,1)+Cog(1)+tem_Dmc(1)
    Node0(:,:,2)=Node0(:,:,2)+Cog(2)+tem_Dmc(2)
    Node0(:,:,3)=Node0(:,:,3)+Cog(3)+tem_Dmc(3)
    
    !**阈值节点瞬时位置**
    breaknode0=InstBreakNode       !第一层节点与第二层节点中间的节点
    
    breaknode0(:,1)=breaknode0(:,1)-Cog(1);
    breaknode0(:,2)=breaknode0(:,2)-Cog(2);
    breaknode0(:,3)=breaknode0(:,3)-Cog(3);
    
    do j=1,Nl
        breaknode0(j,1:3)=VectorL2G( breaknode0(j,1:3),tem_Dmc(4:6) )
    end do
    
    !阈值节点瞬时位置
    breaknode0(:,1)=breaknode0(:,1)+Cog(1)+tem_Dmc(1)
    breaknode0(:,2)=breaknode0(:,2)+Cog(2)+tem_Dmc(2)
    breaknode0(:,3)=breaknode0(:,3)+Cog(3)+tem_Dmc(3)
    
    !**根据瞬时节点位置，计算瞬时节点位置对应入射波波面高度**
    Inst_feta=0.0;
    do i=1,Nwh
    do j=1,Nl
            tmpx(1:3)=Node0(i,j,1:3)
            !20231113修改
            if( IrreCtrl==0 )then
                tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*t   
                Inst_feta(i,j)=amp*cos(tmpt )
            
            else
                Inst_feta(i,j)=0.0;
                do jj=1,IrreNum
                    tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*t + Irrepha(jj)
                    Inst_feta(i,j)=Inst_feta(i,j)+Irreamp(jj)*cos( tmpt )   
                end do                
            end if
            
    end do
    end do
    
    Inst_feta=Inst_feta*smtf       !瞬时船体节点，同x,y位置处对应的瞬时波面起伏
    
    breakInst_feta=0.0;
    do j=1,Nl
        tmpx(1:3)=breaknode0(j,1:3);
!20231113修改
        if( IrreCtrl==0 )then
            tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*t   
            breakInst_feta(j)=amp*cos(tmpt )
        else
            Inst_feta(i,j)=0.0;
            do jj=1,IrreNum
                tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*t + Irrepha(jj)
                breakInst_feta(j)=breakInst_feta(j)+Irreamp(jj)*cos( tmpt )     
            end do
            
        end if
        
    end do
    breakInst_feta=breakInst_feta*smtf   !瞬时阈值节点，同x,y位置处对应的瞬时波面起伏
    
    
    !**划分后的瞬时湿表面网格**
    
    Node1=0.0;     !瞬时湿表面网格(平动坐标系下)
    
    !Nwh:径向网格数；Nl：周向网格数；Node0：瞬时船体网格坐标；
    !Inst_feta:瞬时船体网格，同x，y位置处对应的波面高度
    !breaknode0:瞬时阀值节点坐标；breakInst_feta：瞬时阀值节点坐标，同x，y位置处对应的波面高度
    
    !Node1:划分后的瞬时湿表面网格; submerge:整船淹湿状况
    
    call get_InstantWater2( Nwh,Nl,Node0,Inst_feta,breaknode0,breakInst_feta,Node1,submerge )
    
    !**得到平动坐标系初始状态下关于瞬时重心的瞬时湿表面网格节点坐标(关于重心)**用于计算弹性体节点模态位移
    
    Node2=0.0;
    do i=1,Nwh
    do j=1,Nl
        vec1(1:3)=Node1(i,j,1:3)-Inst_Cog(1:3)
        vec2(1:3)=VectorG2L(vec1,tem_Dmc(4:6) )
        
        Node2(i,j,1:3)=vec2(1:3)
    end do
    end do
    
    !**平动坐标系下关于瞬时重心**
    Node3=0.0;
    do i=1,Nwh
    do j=1,Nl
        Node3(i,j,1:3)=Node1(i,j,1:3)-Inst_Cog(1:3)
    end do
    end do
    
    !**计算节点模态位移**
    !弹性体
    InstUr2=0.0;  InstDur=0.0;
    
    call Inst_Emesh(Nwh,Nl,Node2,InstUr2,InstDur)     !-----取瞬时节点位移，对应弹性模态的
    
    InstUr=0.0;
    
    !进行坐标转换
    do ii=7,NR
      do i=1,Nwh
      do j=1,Nl
          !-----将平动坐标系初始位置的弹性模态节点位移作为瞬时节点上弹性模态的位移
          !-----再将得到的位移整体旋转一个纵倾角，以和瞬时网格对应
          InstUr2(i,j,ii,1:3)=VectorL2G(InstUr2(i,j,ii,1:3),tem_Dmc(4:6) )

          InstUr(i,j,ii,1:3)=InstUr2(i,j,ii,1:3);
      end do
      end do
    end do
    
    !刚体模态对应节点位移
    do i=1,Nwh
    do j=1,Nl
        !----节点坐标(关于重心)         
        tmpx(1:3)=Node3(i,j,1:3)
        
        !------6个刚体运动模态下节点的位移 
        InstUr(i,j,1,1)=1.0;
        InstUr(i,j,2,2)=1.0;
        InstUr(i,j,3,3)=1.0;
        InstUr(i,j,4,2)=-tmpx(3);    InstUr(i,j,4,3)=tmpx(2);
        InstUr(i,j,5,1)=tmpx(3);     InstUr(i,j,5,3)=-tmpx(1);
        InstUr(i,j,6,1)=-tmpx(2);    InstUr(i,j,6,2)=tmpx(1);

    end do
    end do
    
    !**开始计算瞬时静水恢复力以及瞬时入射波力**
    !注意瞬时时刻广义力的求法，是否需要修改？？
    if( submerge==1 )then
        
        Inst_PIS=0.0;
        Inst_feta=0.0;
        Inst_dtphiI=0.0;
        
        do i=1,Nwh
            do j=1,NL
                !**先计算每个节点对应的波面起伏**
                !20231113修改
                tmpx(1:3)=Node1(i,j,1:3)
                if( IrreCtrl==0 )then
                    tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*t   
                    Inst_feta(i,j)=amp*cos(tmpt )
                else
                    Inst_feta(i,j)=0.0;
                    do jj=1,IrreNum
                        tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*t + Irrepha(jj)
                        Inst_feta(i,j)=Inst_feta(i,j)+Irreamp(jj)*cos( tmpt )  
                    end do
                    
                end if
                
                Inst_feta(i,j)=Inst_feta(i,j)*smtf
                
                !**计算偏phi偏t，强制限制重新划分的节点不准超过波面起伏高度**
                if( tmpx(3)>Inst_feta(i,j) )then
                    
                    Inst_PIS(i,j)=0.0;
                    
                else
                    
                    if( tmpx(3)>=0.0 )then     !超过静水面的部分
                        !20231113修改
                        if( IrreCtrl==0 )then
                            tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)-omee*t
                            Inst_dtphiI(i,j)=amp*g0*cos( tmpt )
                        
                        else
                            Inst_dtphiI(i,j)=0.0;
                            do jj=1,IrreNum
                                tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)-Irreomee(jj)*t+Irrepha(jj)
                                Inst_dtphiI(i,j)=Inst_dtphiI(i,j)+Irreamp(jj)*g0*cos( tmpt ) 
                            end do
  
                        end if
                        
                        Inst_dtphiI(i,j)=Inst_dtphiI(i,j)*smtf
                        Inst_dtphiI(i,j)=Inst_dtphiI(i,j)-g0*tmpx(3)
                        Inst_PIS(i,j)=-Inst_dtphiI(i,j)     !---注意正负号
     
                    else
                        !20231113修改
                        if( IrreCtrl==0 )then
                            tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)-omee*t
                            Inst_dtphiI(i,j)=amp*g0*exp(wavek*tmpx(3))*cos( tmpt )
                        
                        else
                            Inst_dtphiI(i,j)=0.0;
                            do jj=1,IrreNum
                                tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)-Irreomee(jj)*t+Irrepha(jj)
                                Inst_dtphiI(i,j)=Inst_dtphiI(i,j)+Irreamp(jj)*g0*exp(Irrek(jj)*tmpx(3))*cos( tmpt )    
                            end do
                            
                        end if
                        
                        Inst_dtphiI(i,j)=Inst_dtphiI(i,j)*smtf
                        Inst_dtphiI(i,j)=Inst_dtphiI(i,j)-g0*tmpx(3)
                        Inst_PIS(i,j)=-Inst_dtphiI(i,j)     !---注意正负号
       
                    end if
        
                end if

            end do
   
        end do
        
        ForceIS=0.0;
        call SurfaceIntegral1(Nwh,Nl,Node3,Inst_PIS,ForceIS,InstUr )
 
    else
        Inst_PIS=0.0;
        ForceIS=0.0;
    end if
    
    
   return
    
    
end subroutine Instant_wetsurface