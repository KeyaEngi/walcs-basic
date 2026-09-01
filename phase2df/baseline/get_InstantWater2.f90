!输入：Nwh:径向网格数；Nl：周向网格数；Node0：瞬时船体网格坐标；
!Inst_feta:瞬时船体网格，同x，y位置处对应的波面高度；
!breaknode0:瞬时阀值节点坐标；breakInst_feta：瞬时阀值节点坐标，同x，y位置处对应的波面高度
!breakKexi=1.0\2.0
    
!输出：Node1:划分后的瞬时水线网格；  submerge:整船淹湿状况
    
    

subroutine get_InstantWater2( Nwh,Nl,Node0,Inst_feta,breaknode0,breakInst_feta,Node1,submerge )
    
    use PanelGeometry,only:breakKexi
    implicit none

    integer(4)::Nwh,Nl
    real(8),dimension(Nwh,Nl,3)::Node0
    real(8),dimension(Nwh,Nl)::Inst_feta
    real(8),dimension(Nl,3)::breaknode0    !----节点阈值
    real(8),dimension(Nl)::breakInst_feta

    real(8),dimension(Nwh,Nl,3)::Node1   !---划分后的瞬时水线网格
    integer(4)::submerge      !------整船淹湿状况


    real(8),dimension(Nwh,Nl,3)::InstP0   !----筛选的水线以下的节点
    integer(4),dimension(Nl)::NumInstP0

    real(8),dimension(Nwh,3)::temInsP0
    integer(4)::temNumInsP0
    real(8),dimension(Nwh)::length
    real(8)::delt

    real(8)::kexi

    integer(4)::i,j,k,ii,jj,kk,iii,jjj,kkk
    real(8)::Z0,Z01,Z02,s,t,z1,z2,x1,x2,y1,y2,x3,y3,z3,ss,tt
    real(8),dimension(3)::temx1,temx2



    !-------先判断节点与波面的相对位置
    NumInstP0=0;
    InstP0=0.0;

    do j=1,Nl      !对周向网格循环
        
        if( breaknode0(j,3)>=breakInst_feta(j) ) then  !----此时认为阈值节点已出水
            NumInstP0(j)=2;
            InstP0(1,j,1:3)=Node0(1,j,1:3);
            InstP0(2,j,1:3)=breaknode0(j,1:3);
        else

            !------从1到Nwh寻找交点
            kk=0;
            do i=1,Nwh
                Z0=Node0(i,j,3)   !----节点高度
                Z01=Inst_feta(i,j)  !---节点对应水面高度

                if( Z01>=Z0 ) then
                    kk=kk+1;
                else
                    exit
                end if
            end do

            if( kk==0 ) then
                NumInstP0(j)=2;
                InstP0(1,j,1:3)=Node0(1,j,1:3);
                InstP0(2,j,1:3)=breaknode0(j,1:3);    
            elseif( kk==1 ) then
                !----根据kk找交点
                z1=Node0( kk,j,3 );    s=Inst_feta( kk,j );
                z2=Node0( kk+1,j,3 );  t=Inst_feta( kk+1,j );

                ss=abs( s-z1 );
                tt=abs( z2-t );

                if(ss+tt<=1.0e-6) then
                    Kexi=0.5;
                else
                    Kexi=ss/(ss+tt);
                end if
                
                if( Kexi<=breakKexi ) then
                    NumInstP0(j)=2;
                    InstP0(1,j,1:3)=Node0(1,j,1:3);
                    InstP0(2,j,1:3)=breaknode0(j,1:3);                      
                else
                    NumInstP0(j)=2;
                    InstP0(1,j,1:3)=Node0(1,j,1:3);
                    InstP0(2,j,1:3)=Node0(1,j,1:3)*(1.0-Kexi)+Node0(2,j,1:3)*Kexi
                end if
                      
            elseif( kk==Nwh ) then   !----完全入水
                NumInstP0(j)=Nwh;

                InstP0(1:Nwh,j,1:3)=Node0(1:Nwh,j,1:3);
            else

                !----根据kk找交点
                z1=Node0( kk,j,3 );    s=Inst_feta( kk,j );
                z2=Node0( kk+1,j,3 );  t=Inst_feta( kk+1,j );

                ss=abs( s-z1 );
                tt=abs( z2-t );

                if(ss+tt<=1.0e-6) then
                    Kexi=0.5;
                else
                    Kexi=ss/(ss+tt);
                end if

                do i=1,kk
                    NumInstP0(j)=NumInstP0(j)+1;
                    InstP0(i,j,1:3)=Node0(i,j,1:3);
                end do
                NumInstP0(j)=NumInstP0(j)+1;
                InstP0(kk+1,j,1:3)=Node0(kk,j,1:3)*(1.0-Kexi)+Node0(kk+1,j,1:3)*Kexi
            end if
        end if

        !-------保险
        if( NumInstP0(j)>=2.and.NumInstP0(j)<=Nwh ) then

        else
            NumInstP0(j)=2;
            InstP0(1,j,1:3)=Node0(1,j,1:3);
            InstP0(2,j,1:3)=breaknode0(j,1:3);            
        end if

    end do 

    Node1=0.0;
    !-----开始等弧长划分节点
    do j=1,Nl       !对周向节点进行循环
        
        temNumInsP0=NumInstP0(j);
        temInsP0=0.0;
        temInsP0( 1:temNumInsP0,1:3 )=InstP0( 1:temNumInsP0,j,1:3 )

        !----累加弦长
        length=0.0;
        t=0.0;
        do i=2,temNumInsP0 
            temx1(1:3)=temInsP0(i-1,1:3);
            temx2(1:3)=temInsP0(i,1:3);

            s=sqrt( (temx2(1)-temx1(1))**2.0+(temx2(2)-temx1(2))**2.0+(temx2(3)-temx1(3))**2.0 )
            t=t+s;
            length(i)=t;
        end do

        !-----更新节点
        delt=length( temNumInsP0 )/real( Nwh-1 )  !---等弧长

        Node1(1,j,1:3)=temInsP0(1,1:3);   !----第一个节点
        do i=2,Nwh-1
            t=delt*real( i-1 )

            do k=1,temNumInsP0-1
                ss=length(k);
                tt=length(k+1);

                if( (t-ss)*(t-tt)<=0.0 ) then
                    if( abs(tt-ss)<=1.0e-6 ) then
                        Kexi=0.5;
                    else
                        Kexi=(t-ss)/(tt-ss)
                    end if

                    Node1(i,j,1:3)=temInsP0(k,1:3)*(1.0-Kexi)+temInsP0(k+1,1:3)*Kexi
                    exit
                end if
            end do
        end do

        Node1(Nwh,j,1:3)=temInsP0(temNumInsP0,1:3)

    end do
    submerge=1;

    return
end subroutine get_InstantWater2


