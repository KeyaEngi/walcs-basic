!输入：Num0（半横剖线粗等分弧长节点个数30）；point0（局部坐标系下，半横剖线粗等分弧长节点yz坐标）
    
!输出：Num1（半横剖线粗等分弧长节点进行质量检查并进行圆角过度后的节点数目（大于30））；
       !Point1（半横剖线粗等分弧长节点进行质量检查并进行圆角过度后的节点坐标）

!功能：    
    
    
!----子程序目的，检查提取的粗节点是否需要修改

subroutine Slam_ChkCoarseCurve(Num0,point0,Num1,Point1)

    implicit none

    integer(4)::Num0
    real(8),dimension(Num0,2)::point0

    integer(4)::Num1
    real(8),dimension(Num0*10,2)::Point1
    !-----需要进行圆弧修正的节点选取
    integer(4)::NumCir
    integer(4),dimension(Num0)::ModiPID
    integer(4),dimension(Num0)::ModiPCase

    real(8),dimension(Num0,5,2)::ModiPoint

    real(8)::AngleB

    real(8),dimension(2)::temCp1,temCp2
    real(8)::lk1,lk2,lb1,lb2
    real(8),dimension(2)::temCirCen
    real(8)::temCirR
    real(8)::cb1,cb2

    real(8)::temAngle1,temAngle2
    real(8),dimension(2)::x1,x2,x3
    real(8),dimension(2)::vel1,vel2
    integer(4)::i,j,k,ii,jj,kk,i1,i2
    real(8)::s,t,s1,s2,t1,t2
    real(8)::Kexi
    real(8)::Pi


!!    open(unit=2001,file='ChkCoarse.txt')



    Pi=atan(1.0)*4.0

    AngleB=-sqrt(3.0)/2.0   !----以135度？？作为判断是否需要圆角过度的临界条件（不是150度吗）
    !------思路
    !------根据夹角差来计算
    NumCir=0;      !用来记录需要进行圆角过度的节点数量
    ModiPID=0;     !用来记录需要进行圆角过度的节点编号
    do i=2,Num0-1
        x1(:)=point0(i-1,:);
        x2(:)=point0(i,:);
        x3(:)=point0(i+1,:);

        !-----作向量相乘来判断(单位化)
        vel1(:)=x1(:)-x2(:);
        s=sqrt( vel1(1)**2.0+vel1(2)**2.0 )
        vel1=vel1/s;

        vel2(:)=x3(:)-x2(:);
        s=sqrt( vel2(1)**2.0+vel2(2)**2.0 )
        vel2=vel2/s;

        !------向量点乘
        t=vel1(1)*vel2(1)+vel1(2)*vel2(2);

        if( t<=AngleB ) cycle      !角度小于150度需要过度

        !------确定需要进行圆角过度的节点编号
        NumCir=NumCir+1;
        ModiPID(NumCir )=i

!!        write(2001,"(i8)") i
    end do



    !----开始过度圆角
    ModiPCase=0;
    ModiPoint=0.0;

    do ii=1,NumCir
        i=ModiPID(ii);

        x1(:)=point0(i-1,:);    !----节点
        x2(:)=point0(i,:);
        x3(:)=point0(i+1,:);

        !------先计算两个线段的长度
        s1=sqrt((x1(1)-x2(1))**2.0+(x1(2)-x2(2))**2.0 )
        s2=sqrt((x3(1)-x2(1))**2.0+(x3(2)-x2(2))**2.0 )

        t=min(s1,s2  )
        t=t/5.0*2.0          !-----需要用圆角替换的直线段长度

        !----1,2点之间的圆弧交点
        Kexi=t/s1;
        temCp1(:)=x2(:)*(1.0-Kexi)+x1(:)*Kexi;
        !----2,3之间的圆弧交点
        Kexi=t/s2;
        temCp2(:)=x2(:)*(1.0-Kexi)+x3(:)*Kexi;

!!        write( 2001,"(2(f15.6,1x))" ) (x1(j),j=1,2)
!!        write( 2001,"(2(f15.6,1x))" ) (temCp1(j),j=1,2)
!!        write( 2001,"(2(f15.6,1x))" ) (x2(j),j=1,2)
!!        write( 2001,"(2(f15.6,1x))" ) (temCp2(j),j=1,2)
!!        write( 2001,"(2(f15.6,1x))" ) (x3(j),j=1,2)
!!        stop
        !-----确定圆心
        !-----先求斜率,垂线斜率
        kk=0;
        temCirCen=0.0;
        if( abs(x1(1)-x2(1))<=1.0e-8 ) then   !----其中一条垂线是水平线
            lk1=0.0
            temCirCen(2)=temCp1(2);

            if( abs(x3(1)-x2(1))<=1.0e-8 ) then   !----第二条线的垂线也是水平线

            elseif( abs(x3(2)-x2(2))<=1.0e-8 ) then   !----第二条线的垂线是竖直线
                temCirCen(1)=temCp2(1)
                !----圆心确定后，定半径
                temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                kk=1;
            else
                lk2=tan( atan( (x3(2)-x2(2))/(x3(1)-x2(1)) )+Pi/2.0 )  !---第二条线的垂线斜率
                !----计算对应纵轴交点lb2
                lb2=temCp2(2)-lk2*temCp2(1)

                !---计算圆心
                temCirCen(1)=(temCirCen(2)-lb2 )/lk2
                temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                kk=1;
            end if

        elseif( abs(x1(2)-x2(2))<=1.0e-8 ) then   !----其中一条垂线是竖直线
            temCirCen(1)=temCp1(1)
            
            if( abs(x3(1)-x2(1))<=1.0e-8 ) then   !----第二条线的垂线是水平线
                temCirCen(2)=temCp2(2)
                temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                kk=1;
            elseif( abs(x3(2)-x2(2))<=1.0e-8 ) then   !----第二条线的垂线是竖直线

            else
                lk2=tan( atan( (x3(2)-x2(2))/(x3(1)-x2(1)) )+Pi/2.0 )  !---第二条线的垂线斜率
                !----计算对应纵轴交点lb2
                lb2=temCp2(2)-lk2*temCp2(1)

                !---计算圆心
                temCirCen(2)=lk2*temCirCen(1)+lb2
                temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                kk=1;
            end if

        else
            lk1=tan( atan( (x1(2)-x2(2))/(x1(1)-x2(1)) )+Pi/2.0 )
            !-----计算第一条垂线对应纵轴交点
            lb1=temCp1(2)-temCp1(1)*lk1;

            if( abs(x3(1)-x2(1))<=1.0e-8 ) then   !----第二条线的垂线是水平线
                temCirCen(2)=temCp2(2);


                temCirCen(1)=( temCirCen(2)-lb1 )/lk1;   !---确定圆心
                temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                kk=1;
            elseif( abs(x3(2)-x2(2))<=1.0e-8 ) then   !----第二条线的垂线是竖直线
                temCirCen(1)=temCp2(1);

                temCirCen(2)=lk1*temCirCen(1)+lb1;
                temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                kk=1;
            else
                
                lk2=tan( atan( (x3(2)-x2(2))/(x3(1)-x2(1)) )+Pi/2.0 )  !---第二条线的垂线斜率
                !----计算对应纵轴交点lb2
                lb2=temCp2(2)-lk2*temCp2(1)

                if( abs(lk2-lk1)<=1.0e-8 ) then
                    
                else
                    !----计算交点
                    temCirCen(1)=(lb1-lb2)/(lk2-lk1)
                    temCirCen(2)=lk2*temCirCen(1)+lb2
                    temCirR=sqrt( (temCirCen(1)-temCp2(1) )**2.0+(temCirCen(2)-temCp2(2) )**2.0 )
                    kk=1;
                end if
            end if

        end if

        !------至此没问题

        if( kk==0 ) then
          ModiPCase(ii)=0;
          cycle
        end if


        ModiPCase(ii)=1;


        !----确定圆弧边界对应的夹角(从x轴正向开始计算)
        vel1(:)=temCp1(:)-temCirCen(:);


        if( abs(vel1(1))<=1.0e-8 ) then   !----数值轴上
            if( vel1(2)>=0.0 ) then
              cb1=Pi/2.0
            else
              cb1=Pi/2.0*3.0
            end if
        else
            if( abs(vel1(2))<=1.0e-8  ) then   !---水平轴上
                if( vel1(1)<0.0 ) then
                  cb1=Pi
                else
                  cb1=0.0
                end if
            else
                lk1=atan( abs(vel1(2)/vel1(1)) )
                if( vel1(1)>0.0.and.vel1(2)>0.0 ) then
                  cb1=lk1;
                elseif( vel1(1)<0.0.and.vel1(2)>0.0 ) then
                  cb1=Pi-lk1;
                elseif( vel1(1)<0.0.and.vel1(2)<0.0 ) then
                  cb1=Pi+lk1;
                elseif( vel1(1)>0.0.and.vel1(2)<0.0 ) then
                  cb1=2.0*Pi-lk1;
                end if
            end if
        end if


        vel2(:)=temCp2(:)-temCirCen(:);
        if( abs(vel2(1))<=1.0e-8 ) then   !----数值轴上
            if( vel2(2)>=0.0 ) then
              cb2=Pi/2.0
            else
              cb2=Pi/2.0*3.0
            end if
        else
            if( abs(vel2(2))<=1.0e-8  ) then   !---水平轴上
                if( vel2(1)<0.0 ) then
                  cb2=Pi
                else
                  cb2=0.0
                end if
            else
                lk2=atan( abs(vel2(2)/vel2(1)) )
                if( vel2(1)>0.0.and.vel2(2)>0.0 ) then
                  cb2=lk2;
                elseif( vel2(1)<0.0.and.vel2(2)>0.0 ) then
                  cb2=Pi-lk2;
                elseif( vel2(1)<0.0.and.vel2(2)<0.0 ) then
                  cb2=Pi+lk2;
                elseif( vel2(1)>0.0.and.vel2(2)<0.0 ) then
                  cb2=2.0*Pi-lk2;
                end if
            end if
        end if


        !-----
        if( cb2-cb1>Pi ) then
            cb2=cb2-2.0*Pi
        elseif( cb2-cb1<-Pi ) then
            cb1=cb1-2.0*Pi
        end if

        !----开始插值节点
        


        do j=1,3
          Kexi=real(j)/4.0

          s=cb1*(1.0-Kexi)+cb2*Kexi

          ModiPoint(ii,j+1,1)=temCirR*cos(s)+temCirCen(1)
          ModiPoint(ii,j+1,2)=temCirR*sin(s)+temCirCen(2)

        end do

        ModiPoint(ii,1,:)=temCp1(:);
        ModiPoint(ii,5,:)=temCp2(:);

    end do

    !----开始替换节点
    Point1=0.0;
    if( NumCir==0 ) then
        
        Point1(1:Num0,:)=Point0(1:Num0,:)
        Num1=Num0;

    else

        kk=0;
        do i=1,NumCir
            ii=ModiPID(i)

            if( i==1 ) then
                i1=1;
                i2=ii-1;
            else
                i1=ModiPID(i-1)+1;
                i2=ii-1;
            end if
            !----前面的节点
            do jj=i1,i2
              kk=kk+1;
              Point1(kk,:)=point0(jj,:)
            end do

            !----需要替换的点
            if( ModiPCase(i)==0 ) cycle
            do jj=1,5
                kk=kk+1;
                Point1(kk,:)=ModiPoint(i,jj,:)
            end do
        end do

        !----后续的点加上
        i1=ModiPID(NumCir)+1;
        i2=Num0;

        do jj=i1,i2
            kk=kk+1;
            Point1(kk,:)=point0(jj,:)
        end do

        Num1=kk;
    end if



    return

end subroutine Slam_ChkCoarseCurve