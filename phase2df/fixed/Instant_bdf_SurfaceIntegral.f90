

!----对重新划分得到的bdf型单元进行压力积分
!----积分方式 单元平均压力积分(以单元中心点处压力为参考值)
subroutine Instant_bdf_SurfaceIntegral(NumN,eleKind,ele_PIS,ele_node,ele_Ur,Force)

    use ArrayOperations
    
    use Constant,only: rou,NR

    implicit none

    integer(4)::NumN   !---单元数
    integer(4),dimension(NumN)::eleKind   !----单元类型
    real(8),dimension(NumN,4)::ele_PIS  !---单元节点处合成压力
    real(8),dimension(NumN,4,3)::ele_node !---单元节点坐标
    real(8),dimension(NumN,4,NR,3)::ele_Ur  !---单元法节点在各模态下的位移
    real*8 :: Force(1:NR)

    real(8),dimension(NR)::temForce

    real(8),dimension(4,3)::temnode4
    real(8),dimension(4)::temPIS4
    real(8),dimension(4,NR,3)::temUr4
    real(8)::temarea
    real(8)::temP
    integer(4)::i,j,k
    
    !-----注意，所有bdf单元法向均已调整指向船内部

    Force=0.0;
    do i=1,NumN
         
         if( eleKind(i)==3  ) then
             temnode4=0.0; temPIS4=0.0; temUr4=0.0; temUr4=0.0;
             temnode4(1:3,1:3)=ele_node(i,1:3,1:3)   !----节点坐标
             temPIS4(1:3)=ele_PIS(i,1:3)   !----节点压力
             temUr4(1:3,1:NR,1:3)=ele_Ur(i,1:3,1:NR,1:3)   !---节点模态位移
             
             call ele3node_SurfaceIntegral( NR,temnode4(1:3,:),temPIS4(1:3),temUr4(1:3,:,:),temForce  )

             Force(:)=Force(:)+temForce(:);
         elseif( eleKind(i)==4 ) then
             temnode4=0.0; temPIS4=0.0; temUr4=0.0; temUr4=0.0;
             temnode4(1:4,1:3)=ele_node(i,1:4,1:3)   !----节点坐标
             temPIS4(1:4)=ele_PIS(i,1:4)   !----节点压力
             temUr4(1:4,1:NR,1:3)=ele_Ur(i,1:4,1:NR,1:3)   !---节点模态位移
             !-----单个4边形积分得到的力
             call ele4node_SurfaceIntegral( NR,temnode4,temPIS4,temUr4,temForce  )

             Force(:)=Force(:)+temForce(:);

         end if
    end do

    Force=rou*Force

    return
end subroutine Instant_bdf_SurfaceIntegral


!------三角形单元积分
subroutine ele3node_SurfaceIntegral( NR,node,Ps,Ur,temForce )
    
    use verification 
    implicit none
    integer(4)::NR
    real(8),dimension(3,3)::node
    real(8),dimension(3)::Ps
    real(8),dimension(3,NR,3)::Ur
    real(8),dimension(NR)::temForce

    real(8),dimension(3)::pasiq
    real(8),dimension(3,NR,3)::eur

    real(8),dimension(NR,3)::Ceur
    real(8)::Cp
    real(8),dimension(3)::normal
    real(8)::area
    real(8)::pl
    
    real(8),dimension(3)::ls
    integer(4)::i,j,k

    !-----单元节点坐标
    x1=node(1,1);y1=node(1,2);z1=node(1,3)
    x2=node(2,1);y2=node(2,2);z2=node(2,3)
    x3=node(3,1);y3=node(3,2);z3=node(3,3)

    pasiq(1:3)=Ps(1:3)   !---节点压力
    eur(1:3,1:NR,1:3)=Ur(1:3,1:NR,1:3)   !---节点位移

    !------计算法向量
    DxDxi=x2-x1;DyDxi=y2-y1;DzDxi=z2-z1
    DxDet=x3-x1;DyDet=y3-y1;DzDet=z3-z1
    hs1_x(1) = DyDxi*DzDet - DyDet*DzDxi
    hs1_y(1) = DzDxi*DxDet - DzDet*DxDxi
    hs1_z(1) = DxDxi*DyDet - DxDet*DyDxi            
    hs1(1)=(hs1_x(1)**2+hs1_y(1)**2+hs1_z(1)**2)**0.5

    if( abs( hs1(1) )<=1.0e-6 ) then
        normal(1:3)=0.0
    else
        normal(1)=hs1_x(1)/hs1(1)   !-----法向量单位化
        normal(2)=hs1_y(1)/hs1(1)
        normal(3)=hs1_z(1)/hs1(1)
    end if

    !------确定单元面积
    ls(1)=sqrt( (x2-x1)**2.0+(y2-y1)**2.0+(z2-z1)**2.0 )
    ls(2)=sqrt( (x3-x2)**2.0+(y3-y2)**2.0+(z3-z2)**2.0 )
    ls(3)=sqrt( (x1-x3)**2.0+(y1-y3)**2.0+(z1-z3)**2.0 )

    pl=( ls(1)+ls(2)+ls(3) )/2.0
    !-----单元面积
    area=pl*(pl-ls(1) )*(pl-ls(2) )*(pl-ls(3) )
    if( area<0.0 ) area=0.0
    area=sqrt(area )

    !-----确定中心点振型位移
    Ceur=0.0;
    Ceur(:,:)=( eur(1,:,:)+eur(2,:,:)+eur(3,:,:) )/3.0
    !-----确定中心点压力
    Cp=(pasiq(1)+pasiq(2)+pasiq(3) )/3.0

    temForce=0.0;
    Do k=1,NR
        temForce(k)=temForce(k)-( Ceur(k,1)*normal(1)+Ceur(k,2)*normal(2)+Ceur(k,3)*normal(3) )*area*Cp
    end do
    !------注意，还未乘以rou

    return
end subroutine ele3node_SurfaceIntegral


!-------四边形单元积分(含模态,平面4点等参元+2节点Gauss积分)
subroutine ele4node_SurfaceIntegral(NR,node,Ps,Ur,temForce)
    
    use verification
    use GauLegCoe
    implicit none

    integer(4)::NR
    real(8),dimension(4,3)::node
    real(8),dimension(4)::Ps
    real(8),dimension(4,NR,3)::Ur

    real(8),dimension(NR)::temForce

    real(8),dimension(4)::pasiq
    real(8),dimension(4,NR,3)::eur

    real(8)::Preq
    real(8),dimension(NR,3)::eurq
    
    real(8),dimension(3)::vanq
    integer(4)::i,j,k,jj,j1

    !-----单元节点坐标
    x1=node(1,1);y1=node(1,2);z1=node(1,3)
    x2=node(2,1);y2=node(2,2);z2=node(2,3)
    x3=node(3,1);y3=node(3,2);z3=node(3,3)
    x4=node(4,1);y4=node(4,2);z4=node(4,3)

    pasiq(1:4)=Ps(1:4)   !---节点压力
    eur(1:4,1:NR,1:3)=Ur(1:4,1:NR,1:3)   !---节点位移

    DxDxi=x2-x1;DyDxi=y2-y1;DzDxi=z2-z1
    DxDet=x4-x1;DyDet=y4-y1;DzDet=z4-z1
    hs1_x(1) = DyDxi*DzDet - DyDet*DzDxi
    hs1_y(1) = DzDxi*DxDet - DzDet*DxDxi
    hs1_z(1) = DxDxi*DyDet - DxDet*DyDxi            
    hs1(1)=(hs1_x(1)**2+hs1_y(1)**2+hs1_z(1)**2)**0.5
        
    DxDxi=x3-x2;DyDxi=y3-y2;DzDxi=z3-z2
    DxDet=x1-x2;DyDet=y1-y2;DzDet=z1-z2
    hs1_x(2) = DyDxi*DzDet - DyDet*DzDxi
    hs1_y(2) = DzDxi*DxDet - DzDet*DxDxi
    hs1_z(2) = DxDxi*DyDet - DxDet*DyDxi            
    hs1(2)=(hs1_x(2)**2+hs1_y(2)**2+hs1_z(2)**2)**0.5
          
    DxDxi=x4-x3;DyDxi=y4-y3;DzDxi=z4-z3
    DxDet=x2-x3;DyDet=y2-y3;DzDet=z2-z3
    hs1_x(3) = DyDxi*DzDet - DyDet*DzDxi
    hs1_y(3) = DzDxi*DxDet - DzDet*DxDxi
    hs1_z(3) = DxDxi*DyDet - DxDet*DyDxi            
    hs1(3)=(hs1_x(3)**2+hs1_y(3)**2+hs1_z(3)**2)**0.5
        
    DxDxi=x1-x4;DyDxi=y1-y4;DzDxi=z1-z4
    DxDet=x3-x4;DyDet=y3-y4;DzDet=z3-z4
    hs1_x(4) = DyDxi*DzDet - DyDet*DzDxi
    hs1_y(4) = DzDxi*DxDet - DzDet*DxDxi
    hs1_z(4) = DxDxi*DyDet - DxDet*DyDxi            
    hs1(4)=(hs1_x(4)**2+hs1_y(4)**2+hs1_z(4)**2)**0.5

       temForce=0.0;
	   Do jj=1,2
	   Do j1=1,2
		   ph1=(1.+cor2(j1))*(1.+cor2(jj))/4.   !----形函数
		   ph2=(1.-cor2(j1))*(1.+cor2(jj))/4.
		   ph3=(1.-cor2(j1))*(1.-cor2(jj))/4.
		   ph4=(1.+cor2(j1))*(1.-cor2(jj))/4.           

           xq=ph1*x1+ph2*x2+ph3*x3+ph4*x4
	       yq=ph1*y1+ph2*y2+ph3*y3+ph4*y4
           zq=ph1*z1+ph2*z2+ph3*z3+ph4*z4

           hs_x=(ph1*hs1_x(1)+ph2*hs1_x(2)+ph3*hs1_x(3)+ph4*hs1_x(4))/4.        
           hs_y=(ph1*hs1_y(1)+ph2*hs1_y(2)+ph3*hs1_y(3)+ph4*hs1_y(4))/4.         
           hs_z=(ph1*hs1_z(1)+ph2*hs1_z(2)+ph3*hs1_z(3)+ph4*hs1_z(4))/4. 
           hs=(ph1*hs1(1)+ph2*hs1(2)+ph3*hs1(3)+ph4*hs1(4))/4.

           !-----压力积分节点
           Preq=ph1*pasiq(1)+ph2*pasiq(2)+ph3*pasiq(3)+ph4*pasiq(4)
           !-----模态位移积分节点
           eurq(:,:)=ph1*eur(1,:,:)+ph2*eur(2,:,:)+ph3*eur(3,:,:)+ph4*eur(4,:,:)

           vanq(1)=hs_x;
           vanq(2)=hs_y;
           vanq(3)=hs_z;

           Do k=1,NR
               temForce(k)=temForce(k)-ww(j1)*ww(jj)*( eurq(k,1)*vanq(1)+eurq(k,2)*vanq(2)+eurq(k,3)*vanq(3) )*Preq
           end do

       end do
       end do
       !------注意，此时还未乘以rou

    return
end subroutine ele4node_SurfaceIntegral













