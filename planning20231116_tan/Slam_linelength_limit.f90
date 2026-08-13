!输入：SlamNumLine（典型砰击节点数目）；SlamNumLP（每条砰击曲线拥有的型值点数，细化后）；SlamNode（砰击节点左舷坐标，用户坐标系关于重心）
!SlamAngle（考虑到船体纵倾影响的砰击角）；SlamLibType（典型砰击节点坐标）；SlamLineCase（砰击剖线生成情况，对应典型砰击节点位置是否生成砰击剖线）    
    

!获取用户坐标系中各砰击剖线起始节点坐标，终止节点坐标temSeffectP   
    
    
!----砰击剖线长短限制模块
!----思路 以船中为分界线
!----前部区域从艏到中依次判断长短
!----后部区域从艉到中依次判断长短

subroutine Slam_linelength_limit()
    
    use Slamming,only:SlamNumLine,SlamNumLP,SlamNode,SlamAngle,SlamLibType,SlamLineCase
    use Constant,only:Pi
    implicit none

    !-----剖线右舷点
    real(8),dimension( SlamNumLine,2,3 )::temSeffectP     !各砰击剖线起始节点坐标，终止节点坐标

    real(8)::midSlamline
    integer(4)::NumLP

    real(8)::temAng1,temAng2

    real(8),dimension( SlamNumLine,3 )::teminteractP
    !----交点存在状态 0不存在 1存在
    integer(4),dimension(SlamNumLine)::teminteractcase

    integer(4)::limitcase

    real(8)::temx
    real(8),dimension(3)::x1,x2,x3,x4
    real(8)::lk1,lk2,lb1,lb2
    integer(4)::i,j,k
    real(8)::s,t,s1,s2
    real(8)::Kexi


    midSlamline=(SlamLibType(1,1)+SlamLibType(SlamNumLine,1) )/2.0      !大致估计船中x位置

    temSeffectP=0.0;

!***************************************
    !-----先确定船体前部剖线实际长短限制
    !-----艏部第一条剖线
    if(SlamLineCase(SlamNumLine)==0) then      !对应位置未生成砰击剖线
        temSeffectP(SlamNumLine,1,1:3)=SlamLibType(SlamNumLine,1:3);     !船首第一条剖线对应位置
        temSeffectP(SlamNumLine,2,1:3)=temSeffectP(SlamNumLine,1,1:3);
    else        !对应位置成功生成砰击剖线
        temSeffectP(SlamNumLine,1,1:3)=SlamLibType(SlamNumLine,1:3);    !船首第一条剖线起始节点坐标
        NumLP=SlamNumLP(SlamNumLine);
        temSeffectP(SlamNumLine,2,1:3)=SlamNode(SlamNumLine,NumLP,1:3);     !船首第一条剖线终止节点坐标
    end if


    do i=SlamNumLine-1,1,-1      !从首部到尾部进行循环  
        temx=SlamLibType(i,1);     !典型砰击节点x坐标

        if( temx<midSlamline ) exit   !中部位置限制

        if( SlamLineCase(i)==0 ) then       !未成功生成剖线
            temSeffectP(i,1,1:3)=SlamLibType(i,1:3);
            temSeffectP(i,2,1:3)=temSeffectP(i,1,1:3)
            cycle
        end if

        temAng1=SlamAngle(i);    !考虑船体纵倾影响的砰击角

        !----先给定第一个有效限制点
        temSeffectP(i,1,1:3)=SlamLibType(i,1:3);


        !--先确定第i条剖线与第i+1到第SlamNumLine条剖线的交点
        teminteractP=0.0;
        teminteractcase=0;     !记录第i条剖线与第j条剖线是否有交点，且交点在第j条剖线有效区间内（0无1有）
        do j=i+1,SlamNumLine
            temAng2=SlamAngle(j)

            x1=0.0;  !两剖线交点
            x2=temSeffectP(i,1,1:3);
            x3=temSeffectP(j,1,1:3);
            if( abs(temAng1-temAng2)<=1.0e-4  ) then    !两条剖线平行，无交点
                teminteractcase(j)=0;
            else
                if( abs(temAng1-Pi/2.0)<=1.0e-4 ) then
                    x1(1)=temSeffectP(i,1,1);

                    lk2=tan( temAng2 );
                    lb2=x3(3)-lk2*x3(1);

                    x1(3)=lk2*x1(1)+lb2
                elseif( abs(temAng2-Pi/2.0)<=1.0e-4 ) then
                    x1(1)=temSeffectP(j,1,1);

                    lk1=tan( temAng1 );
                    lb1=x2(3)-lk1*x2(1);

                    x1(3)=lk1*x1(1)+lb1;
                else
                    !----两条斜线
                    lk1=tan( temAng1 );
                    lb1=x2(3)-lk1*x2(1);                    
                    lk2=tan( temAng2 );
                    lb2=x3(3)-lk2*x3(1);

                    x1(1)=(lb2-lb1)/(lk1-lk2)
                    x1(3)=lk1*x1(1)+lb1
                end if

                !----确定交点是否在j条剖线的有效区间内
                s=temSeffectP(j,1,3);
                t=temSeffectP(j,2,3);
                if( (x1(3)-s)*(x1(3)-t)<=0.0 ) then
                    teminteractcase(j)=1;
                    teminteractP(j,1:3)=x1(1:3);
                else
                    teminteractcase(j)=0;
                end if
            end if

        end do

        !-----确定第i条剖线的上限制点
        limitcase=0;
        x1=0.0;
        do j=i+1,SlamNumLine
            
            if( teminteractcase(j)==1.and.limitcase==0 ) then
                limitcase=1;
                x1(1:3)=teminteractP(j,1:3);
            elseif( teminteractcase(j)==1 ) then
                if( x1(3)>teminteractP(j,3) ) then
                    x1(1:3)=teminteractP(j,1:3);    !由此找出最低交点
                end if
            end if
        end do
        !-----第i条剖线的最低相交点确定

        !-----开始确定限制点
        x2=0.0;
        if( limitcase==1 ) then
            NumLP=SlamNumLP(i);
            x2(1:3)=SlamNode(i,NumLP,1:3);  

            if( x2(3)<=x1(3) ) then
                temSeffectP(i,2,1:3)=x2(1:3)
            else
                do j=NumLP-1,1,-1
                    s1=SlamNode(i,j,3);
                    s2=SlamNode(i,j+1,3);
                    if( (x1(3)-s1)*(x1(3)-s2)<=0.0 ) then
                        SlamNumLP(i)=j+1;
                        temSeffectP(i,2,1:3)=x1(1:3)
                        exit
                    end if
                end do

            end if
        else
            NumLP=SlamNumLP(i);
            temSeffectP(i,2,1:3)=SlamNode(i,NumLP,1:3);
        end if



    end do


!******************************************
!-----确定船体后部剖线的长短限制
    !-----艉部第一条剖线
    if(SlamLineCase(1)==0) then
        temSeffectP(1,1,1:3)=SlamLibType(1,1:3);
        temSeffectP(1,2,1:3)=temSeffectP(1,1,1:3);
    else
        temSeffectP(1,1,1:3)=SlamLibType(1,1:3);
        NumLP=SlamNumLP(1);
        temSeffectP(1,2,1:3)=SlamNode(1,NumLP,1:3);
    end if


    do i=2,SlamNumLine
        temx=SlamLibType(i,1);

        if( temx>=midSlamline ) exit

        if( SlamLineCase(i)==0 ) then
            temSeffectP(i,1,1:3)=SlamLibType(i,1:3);
            temSeffectP(i,2,1:3)=temSeffectP(i,1,1:3)
            cycle
        end if

        temAng1=SlamAngle(i);

        !----先给定第一个有效限制点
        temSeffectP(i,1,1:3)=SlamLibType(i,1:3);


        !--先确定第i条剖线与第i+1到第SlamNumLine条剖线的交点
        teminteractP=0.0;
        teminteractcase=0;
        do j=i-1,1,-1
            temAng2=SlamAngle(j)

            x1=0.0;  !---交点
            x2=temSeffectP(i,1,1:3);
            x3=temSeffectP(j,1,1:3);
            if( abs(temAng1-temAng2)<=1.0e-4  ) then
                teminteractcase(j)=0;
            else
                if( abs(temAng1-Pi/2.0)<=1.0e-4 ) then
                    x1(1)=temSeffectP(i,1,1);

                    lk2=tan( temAng2 );
                    lb2=x3(3)-lk2*x3(1);

                    x1(3)=lk2*x1(1)+lb2
                elseif( abs(temAng2-Pi/2.0)<=1.0e-4 ) then
                    x1(1)=temSeffectP(j,1,1);

                    lk1=tan( temAng1 );
                    lb1=x2(3)-lk1*x2(1);

                    x1(3)=lk1*x1(1)+lb1;
                else
                    !----两条斜线
                    lk1=tan( temAng1 );
                    lb1=x2(3)-lk1*x2(1);                    
                    lk2=tan( temAng2 );
                    lb2=x3(3)-lk2*x3(1);

                    x1(1)=(lb2-lb1)/(lk1-lk2)
                    x1(3)=lk1*x1(1)+lb1
                end if

                !----确定交点是否在j条剖线的有效区间内
                s=temSeffectP(j,1,3);
                t=temSeffectP(j,2,3);
                if( (x1(3)-s)*(x1(3)-t)<=0.0 ) then
                    teminteractcase(j)=1;
                    teminteractP(j,1:3)=x1(1:3);
                else
                    teminteractcase(j)=0;
                end if
            end if

        end do

        !-----确定第i条剖线的上限制点
        limitcase=0;
        x1=0.0;
        do j=i-1,1,-1
            
            if( teminteractcase(j)==1.and.limitcase==0 ) then
                limitcase=1;
                x1(1:3)=teminteractP(j,1:3);
            elseif( teminteractcase(j)==1 ) then
                if( x1(3)>teminteractP(j,3) ) then
                    x1(1:3)=teminteractP(j,1:3);
                end if
            end if
        end do
        !-----第i条剖线的最低相交点确定

        !-----开始确定限制点
        x2=0.0;
        if( limitcase==1 ) then
            NumLP=SlamNumLP(i);
            x2(1:3)=SlamNode(i,NumLP,1:3);  

            if( x2(3)<=x1(3) ) then
                temSeffectP(i,2,1:3)=x2(1:3)
            else
                do j=NumLP-1,1,-1
                    s1=SlamNode(i,j,3);
                    s2=SlamNode(i,j+1,3);
                    if( (x1(3)-s1)*(x1(3)-s2)<=0.0 ) then
                        SlamNumLP(i)=j+1;
                        temSeffectP(i,2,1:3)=x1(1:3)
                        exit
                    end if
                end do

            end if
        else
            NumLP=SlamNumLP(i);
            temSeffectP(i,2,1:3)=SlamNode(i,NumLP,1:3);
        end if



    end do





    return
end subroutine Slam_linelength_limit







