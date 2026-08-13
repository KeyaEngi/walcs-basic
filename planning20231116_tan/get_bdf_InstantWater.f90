

!------bdf网格文件的瞬时湿表面划分
!------思路 对于完全在水面以下的单元和完全出水的单元不处理
!------     对于与水线有交线的单元，先将4边形单元划分为两个三角形单元

subroutine get_bdf_InstantWater(Num0,node0,insfate0,NumNew_ele,NewKind,Newnode  )

    implicit none

    integer(4)::Num0   !-----单元拥有节点数
    real(8),dimension(4,3)::node0   !----节点坐标
    real(8),dimension(4)::insfate0  !----节点对应水面高度
    real(8)::area0
    real(8),dimension(3)::xav0
    real(8),dimension(3)::normal0     !-----初始法向

    integer(4)::NumNew_ele   !---重构的单元数量
    integer(4),dimension(2)::NewKind  !---每个重构单元类型(3或4)
    real(8),dimension(2)::NewArea
    real(8),dimension(2,3)::Newxav
    real(8),dimension(2,3)::NewNormal
    real(8),dimension(2,4,3)::Newnode  !---重构单元的节点坐标

    integer(4)::temeleexist   !----新单元存在状态(0不存在 1存在)
    integer(4)::temelekind    !----新单元拥有节点数目
    real(8),dimension(4,3)::temnewnode   !---新单元节点坐标
    real(8)::temnewelearea    !----新单元投影面积
    real(8),dimension(3)::temnewxav        !----新单元中心点坐标

    integer(4)::numsubmerge

    real(8),dimension(3,3)::temnode
    real(8),dimension(3)::teminsfeta

    integer(4)::temNum0
    integer(4)::i,j,k
    real(8)::s,t


!**********************************
    temNum0=Num0

    numsubmerge=0;
    do i=1,temNum0
        if(node0(i,3)<=insfate0(i)  ) then
            numsubmerge=numsubmerge+1;
        end if
    end do


    NumNew_ele=0;
    NewKind=0;
    Newnode=0.0;
    if( numsubmerge==0  ) then  !完全出水
        NumNew_ele=0;

    elseif( numsubmerge==temNum0 ) then  !完全入水
        NumNew_ele=1;
        !-----新单元拥有的节点数、单元面积、单元中心点坐标、单元法向量
        NewKind(1)=temNum0;
        Newnode(1,1:temNum0,1:3)=node0(1:temNum0,1:3) !--单元节点坐标
    else  !--单元与水线有交点，需要重新划分
        
        !---先将4边形单元划分为两个三角形单元
        if( temNum0==3 ) then
            temnode=0.0;  teminsfeta=0.0;
            temnode(1:3,1:3)=node0(1:3,1:3);
            teminsfeta(1:3)=insfate0(1:3)

            call bdfele_breaktri3(temnode,teminsfeta,temeleexist,temelekind,temnewnode,temnewelearea,&
                  &temnewxav )
            
            if( temeleexist==1 ) then
                NumNew_ele=1;
                NewKind(1)=temelekind;  !--新单元节点数
                Newnode(1,1:temelekind,1:3)=temnewnode(1:temelekind,1:3)
            elseif( temeleexist==0 ) then
                NumNew_ele=0;
            end if

        elseif( temNum0==4 ) then
            !---思路 先将四边形单元划分为两个三角形单元
            do i=1,2
                if(i==1) then
                    temnode=0.0;  teminsfeta=0.0;
                    temnode(1:3,1:3)=node0(1:3,1:3);
                    teminsfeta(1:3)=insfate0(1:3)                    
                elseif(i==2) then
                    temnode=0.0;  teminsfeta=0.0;
                    temnode(1,1:3)=node0(1,1:3);
                    temnode(2,1:3)=node0(3,1:3);
                    temnode(3,1:3)=node0(4,1:3);
                    teminsfeta(1)=insfate0(1);                       
                    teminsfeta(2)=insfate0(3);
                    teminsfeta(3)=insfate0(4);
                end if

                call bdfele_breaktri3(temnode,teminsfeta,temeleexist,temelekind,temnewnode,temnewelearea,&
                      &temnewxav )

                if( temeleexist==1 ) then
                    NumNew_ele=NumNew_ele+1;
                    NewKind(NumNew_ele)=temelekind;  !--新单元节点数
                    Newnode(NumNew_ele,1:temelekind,1:3)=temnewnode(1:temelekind,1:3)
                end if
                
            end do

        end if  

    end if

    return
end subroutine get_bdf_InstantWater

!-----分割三角形单元
subroutine bdfele_breaktri3(node0,zz,eleexist,elekind,newnode,newelearea,&
               &newxav )

    implicit none

    real(8),dimension(3,3)::node0
    real(8),dimension(3)::zz

    integer(4)::eleexist   !----新单元存在状态(0不存在 1存在)
    integer(4)::elekind    !----新单元拥有节点数目
    real(8),dimension(4,3)::newnode   !---新单元节点坐标
    real(8)::newelearea    !----新单元投影面积
    real(8),dimension(3)::newxav        !----新单元中心点坐标

    integer(4)::Numsubmerge
    integer(4),dimension(3)::subnode

    real(8),dimension(3,3)::node1
    real(8),dimension(3)::zz1
    real(8),dimension(3)::x1,x2,x3,x4,x5,x6
    real(8),dimension(5)::ls
    integer(4)::i,j,k
    real(8)::s,t,p1,p2,w1,w2,w
    real(8)::Kexi


    !-----先判断水下节点数目以及各节点与水线的相对位置关系
    !-----0 节点在水线以上 1 节点在水线以下
    Numsubmerge=0;
    subnode=0;

    do i=1,3
        if(node0(i,3)<=zz(i)  ) then
            Numsubmerge=Numsubmerge+1;
            subnode(i)=1;
        end if
    end do

    eleexist=0;
    if(Numsubmerge==0  ) then   !---完全出水
        eleexist=0;
        elekind=0;
        newnode=0.0;
        newelearea=0.0;
        newxav=0.0;
    elseif( Numsubmerge==1 ) then  !---只有一个点在水下,重构的是三角形单元
        eleexist=1;
        elekind=3;

        newnode=0.0;
        !-----标准情况：1号节点在水下，2和3号节点在水上
        if(subnode(1)==1  ) then  !---不需要挪动节点顺序
            node1(:,:)=node0(:,:);
            zz1(:)=zz(:);
        elseif( subnode(2)==1 ) then  !---挪动节点顺序
            node1(1,:)=node0(2,:);
            node1(2,:)=node0(3,:);
            node1(3,:)=node0(1,:);
            zz1(1)=zz(2);
            zz1(2)=zz(3);
            zz1(3)=zz(1);
        elseif( subnode(3)==1 ) then  !---挪动节点顺序
            node1(1,:)=node0(3,:);
            node1(2,:)=node0(1,:);
            node1(3,:)=node0(2,:);
            zz1(1)=zz(3);
            zz1(2)=zz(1);
            zz1(3)=zz(2);            
        end if

        newnode(1,:)=node1(1,:);
        !-----计算交点 1/2边与水线交点为新的2号点
        s=abs( zz1(1)-node1(1,3) )
        t=abs( node1(2,3)-zz1(2) )

        if( abs(s+t)<=1.0e-6 ) then
            Kexi=0.5;
        else
            Kexi=s/(s+t);
        end if

        !----新的2号节点
        newnode(2,:)=node1(1,:)*(1.0-Kexi)+node1(2,:)*Kexi

        !-----1/3边与水线交点为新的3号点
        t=abs( node1(3,3)-zz1(3) )
        
        if( abs(s+t)<=1.0e-6 ) then
            Kexi=0.5;
        else
            Kexi=s/(s+t);
        end if
        newnode(3,:)=node1(1,:)*(1.0-Kexi)+node1(3,:)*Kexi

        newnode(4,:)=newnode(3,:);
        !-----新单元节点坐标构造完毕

    elseif( Numsubmerge==2 ) then   !---有两个点在水线下，重构出来的是四边形单元
        eleexist=1;
        elekind=4;

        newnode=0.0;
        !----标准情况 1号点在水上，2和3号点在水下
        if( subnode(1)==0 ) then
            node1(:,:)=node0(:,:);
            zz1(:)=zz(:);            
        elseif( subnode(2)==0 ) then
            node1(1,:)=node0(2,:);
            node1(2,:)=node0(3,:);
            node1(3,:)=node0(1,:);
            zz1(1)=zz(2);
            zz1(2)=zz(3);
            zz1(3)=zz(1);            
        elseif( subnode(3)==0 ) then
            node1(1,:)=node0(3,:);
            node1(2,:)=node0(1,:);
            node1(3,:)=node0(2,:);
            zz1(1)=zz(3);
            zz1(2)=zz(1);
            zz1(3)=zz(2); 
        end if

        newnode(2,:)=node1(2,:);
        newnode(3,:)=node1(3,:);
        !-----计算交点 1/2边与水线交点为新的1号点
        s=abs( zz1(1)-node1(1,3) )
        t=abs( node1(2,3)-zz1(2) )
        
        if( abs(s+t)<=1.0e-6 ) then
            Kexi=0.5;
        else
            Kexi=s/(s+t);
        end if
        newnode(1,:)=node1(1,:)*(1.0-Kexi)+node1(2,:)*Kexi

        !-----1/3边与水线交点为新的3号点
        t=abs( node1(3,3)-zz1(3) )
        
        if( abs(s+t)<=1.0e-6 ) then
            Kexi=0.5;
        else
            Kexi=s/(s+t);
        end if
        newnode(4,:)=node1(1,:)*(1.0-Kexi)+node1(3,:)*Kexi
        !----新单元节点构造完毕

    elseif( Numsubmerge==3 ) then  !----三个点全部在水面以下，直接使用原始单元信息
        eleexist=1;
        elekind=3;

        newnode=0.0;        
        newnode(1:3,1:3)=node0(1:3,1:3)
        newnode(4,:)=newnode(3,:)
        !----新单元节点构造完毕

    end if

    !-----计算新单元的面积、中心点坐标


    if( eleexist==1 ) then

        x1(:)=newnode(1,:);
        x2(:)=newnode(2,:);
        x3(:)=newnode(3,:);
        x4(:)=newnode(4,:);
        ls(1)=sqrt((x1(1)-x2(1))**2.0+(x1(2)-x2(2))**2.0+(x1(3)-x2(3))**2.0 )
        ls(2)=sqrt((x2(1)-x3(1))**2.0+(x2(2)-x3(2))**2.0+(x2(3)-x3(3))**2.0 )
        ls(3)=sqrt((x3(1)-x4(1))**2.0+(x3(2)-x4(2))**2.0+(x3(3)-x4(3))**2.0 )
        ls(4)=sqrt((x4(1)-x1(1))**2.0+(x4(2)-x1(2))**2.0+(x4(3)-x1(3))**2.0 )
        ls(5)=sqrt((x1(1)-x3(1))**2.0+(x1(2)-x3(2))**2.0+(x1(3)-x3(3))**2.0 )

        p1=( ls(1)+ls(2)+ls(5) )/2.0;
        p2=( ls(3)+ls(4)+ls(5) )/2.0;
        w1=p1*(p1-ls(1))*(p1-ls(2))*(p1-ls(5))
        w2=p2*(p2-ls(3))*(p2-ls(4))*(p2-ls(5))
        if( w1<=0.0 ) w1=0.0;
        if( w2<=0.0 ) w2=0.0;
        w1=sqrt(w1);
        w2=sqrt(w2);

        !----面积
        w=w1+w2;
        !------先判断单元面积是否符合要求
        if( w<=1.0e-7 ) then
            eleexist=0;
            elekind=0;
            newnode=0.0;
            newelearea=0.0;
            newxav=0.0;
        else
            newelearea=w;
        end if

        !----此处较为随意(需要优化)
        if( elekind==3 ) then
            newxav(:)=(x1(:)+x2(:)+x3(:)  )/3.0
        elseif( elekind==4 ) then
            newxav(:)=(x1(:)+x2(:)+x3(:)+x4(:)  )/4.0
        end if

    end if

    return
end subroutine bdfele_breaktri3



