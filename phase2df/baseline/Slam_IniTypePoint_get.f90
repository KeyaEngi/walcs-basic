
!输入：砰击粗选节点数目；粗选节点纵向位置；中纵剖线节点数目；中纵剖线节点坐标；划分砰击倾角个数；
!      划分砰击倾角个数对应的边界；砰击倾角大小；切片带宽修正系数对应的区域段数；切片带宽修正系数对应的区域边界;
!      切片带宽修正系数值
!输出：SlamLibType（典型砰击节点坐标）；SlamNumLine（典型砰击节点数量）；SlamAngle（典型砰击节点对应区间砰击倾角）
!      SlamWidthMcoef（典型砰击节点位置对应带宽修正系数）

!功能：就是通过插值的方式获得中纵剖线上的砰击典型点坐标，并求出典型砰击节点对应的砰击角及带宽修正系数，因为粗选的砰击
    !节点纵坐标可能会超出船长的范围内，所以典型砰击节点数目相对于粗选砰击节点数目会有所减少

!-----模块目的：获得中纵剖线上的砰击典型点坐标

subroutine Slam_IniTypePoint_get(num0,pointx0,NumLp,LPoint,numAngarea,Angboundary,Ang_value,&
                               &  NumSWcof,Sw_boundaryx,SW_coef  )
    
    use Slamming,only:SlamLibType,SlamNumLine,SlamAngle,SlamWidthMcoef
    use Constant,only:Pi,OutAccess
    implicit none

    
    real(8)::temAngle   !---给定斜率
    integer(4)::num0
    real(8),dimension(num0)::pointx0
    integer(4)::NumLp
    real(8),dimension(NumLp,3)::LPoint
    integer(4)::numAngarea    !-----剖线倾角的分布区域份数
    real(8),dimension(numAngarea)::Angboundary     !----各区域的边界
    real(8),dimension(numAngarea)::Ang_value   !---区域内倾角大小(当前的还是deg)
    integer(4)::NumSWcof      !-----剖线切片宽度修正系数分布区域份数
    real(8),dimension( NumSWcof,2 )::Sw_boundaryx    !----各区域边界
    real(8),dimension( NumSWcof )::SW_coef    !---各区域切片宽度修正系数

    integer(4)::NumType
    real(8),dimension(num0,3)::temTypeP
    real(8),dimension(num0)::temTypeAngle

    integer(4)::temNumP
    real(8),dimension(NumLp,3 )::temP

    real(8)::lk,lb
    real(8)::tems
    real(8)::temx
    real(8),dimension(3)::x1,x2,x3,x4
    real(8)::s,t,s1,s2,t1,t2
    integer(4)::i,j,k
    real(8)::Kexi


    !砰击倾角大小先变换为弧度
    do i=1,numAngarea       !对砰击倾角对应的划分区域进行循环
        Ang_value(i)=Ang_value(i)/180.0*Pi
    end do


    NumType=0;
    temTypeP=0.0;
    temTypeAngle=0.0;
    do i=1,num0     !对粗选的砰击节点数进行循环
        temx=pointx0(i);    !粗选砰击节点纵向位置

        !---identify incline degree
        !----定义粗选砰击节点对应的砰击倾角temAngle
        if( numAngarea==1 ) then       !如果砰击倾角对应的划分区域个数只有一个的话
            temAngle=Ang_value(1)    
        else
            
            do j=1,numAngarea
                if( j==1 ) then
                    s2=Angboundary(1)    !第一砰击区域对应的右边界

                    if( temx<=s2  ) then    !所有粗选砰击节点位于第一砰击区域右边界以前的点，砰击角均为第一区域的砰击角
                        temAngle=Ang_value(1)
                        exit
                    end if
                elseif( j==numAngarea ) then
                    s1=Angboundary(numAngarea-1 )    !最后一个砰击区域对应的左边界

                    if( temx>s1  ) then   !所有粗选砰击节点位于最后一个砰击区域左边界以后的点，砰击角均为最后区域的砰击角
                        temAngle=Ang_value(numAngarea)
                        exit
                    end if
                else
                    s1=Angboundary(j-1);     !中间区域左边界
                    s2=Angboundary(j);       !中间区域右边界

                    if(  temx>s1.and.temx<=s2  ) then
                        temAngle=Ang_value(j);
                        exit
                    end if
                end if
            end do
        end if


        temP=0.0;
        temNumP=0;
        do j=1,NumLp-1              !对中纵剖线节点（数目-1）进行循环
            x1(:)=LPoint(j,:);      !中纵剖线节点坐标
            x2(:)=LPoint(j+1,:);

            if( abs(temAngle-Pi/2.0)<=1.0e-4 ) then  !---竖直剖面（区域砰击角为90度）
                s=( temx-x1(1) )*(temx-x2(1) );    
                if(s<=0.0) then          !粗选砰击节点纵向位置在两中纵剖线节点之间，或与任意一中纵剖线节点纵向位置重合
                    temNumP=temNumP+1;       !计数

                    if( x1(1)==x2(1) ) then    !如果三点纵向位置相同，则取垂向坐标较低的点
                        if( x1(3)<=x2(3) ) then     
                            temP(temNumP,:)=x1(:);
                        else
                            temP(temNumP,:)=x2(:);
                        end if
                    else                      !如果三点纵向坐标不完全相同，则线性插值
                        Kexi=(temx-x1(1))/(x2(1)-x1(1));
                        temP(temNumP,:)=x1(:)*(1.0-Kexi)+x2(:)*Kexi;
                    end if
                end if
            else
                !----计算直线斜率lk以及0点截距lb
                lk=tan( temAngle );
                lb=-lk*temx;
                !----确定与x1,x2同x的剖线上的节点坐标
                x3=0.0; x4=0.0;
                x3(1)=x1(1);
                x4(1)=x2(1);
                x3(3)=lk*x3(1)+lb;
                x4(3)=lk*x4(1)+lb;

                s=(x1(3)-x3(3))*(x2(3)-x4(3))
                if( s<=0.0 ) then
                    if( abs(x1(3)-x3(3))+abs(x2(3)-x4(3))==0.0 ) then
                        temNumP=temNumP+1;
                        temP(temNumP,:)=(x1(:)+x2(:))/2.0
                    else
                        Kexi=abs(x1(3)-x3(3))/( abs(x1(3)-x3(3))+abs(x2(3)-x4(3)) )
                        temNumP=temNumP+1;
                        temP(temNumP,:)=x1(:)*(1.0-Kexi)+x2(:)*Kexi;
                    end if
                end if
            end if
        end do

        !----选最下面的交点作为典型节点
        !也就是说一个粗选砰击节点纵坐标最多对应一个典型砰击节点
        if(temNumP==0 ) cycle

        x1(:)=temP(1,:);
        do j=2,temNumP
            if( temP(j,3)<x1(3) ) then
                x1(:)=temP(j,:)
            end if
        end do

        NumType=NumType+1;      !典型砰击节点数量
        temTypeP(NumType,:  )=x1(:);         !典型砰击节点坐标

        temTypeAngle( NumType )=temAngle;        !典型砰击节点对应区域对应的砰击角
    end do

    !-----------确定典型节点数量以及节点位置
    SlamNumLine=NumType;
    allocate( SlamLibType(SlamNumLine,3 ) )

    SlamLibType=0.0;
    SlamLibType(1:SlamNumLine,1:3)=temTypeP(1:SlamNumLine,1:3)


    allocate( SlamAngle(SlamNumLine ) )
    SlamAngle=0.0;
    SlamAngle(1:SlamNumLine )=temTypeAngle(1:SlamNumLine )


    !*************************************
    !-------根据选定的砰击典型节点位置计算对应切片宽度修正系数
    allocate( SlamWidthMcoef(SlamNumLine ) )
    SlamWidthMcoef=1.0

    do i=1,SlamNumLine    !典型砰击节点数量
        
        temx=SlamLibType( i,1 )  !---典型砰击节点纵向位置x

        do j=1,NumSWcof
            s=Sw_boundaryx(j,1);
            t=Sw_boundaryx(j,2);

            if( (temx-s)*(temx-t)<=0.0 ) then
                SlamWidthMcoef(i)=SW_coef(j);      !对应切片宽度修正系数
                exit
            end if
        end do
    end do


    open(unit=4001,file=trim(adjustl(OutAccess))//'\'//'SlamtypePoint.txt')
    write(4001,*) SlamNumLine

    do i=1,SlamNumLine
        write(4001,"(3(f15.6,1x))") (SlamLibType(i,j),j=1,3)
    end do
    close(4001)


    return
end subroutine Slam_IniTypePoint_get








