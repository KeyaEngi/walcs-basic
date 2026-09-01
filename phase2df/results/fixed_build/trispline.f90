
!-----采用普通三次样条方法插值振型函数

subroutine trispline( n,point1,point2,ds1,dsn,NumN,InterX,Ns,Nds,Ndds)

    implicit none

    integer(4)::n    !----节点数
    real(8),dimension(n,1)::point1,point2   !----节点位置
    real(8)::ds1,dsn   !---首尾节点边界条件,一阶导数
    integer(4)::NumN    !----待插节点数
    real(8),dimension(NumN)::InterX   !---待插节点坐标x

    real(8),dimension(NumN)::Ns,Nds,Ndds   !---待插节点处函数值，一阶导，二阶导

    real(8),dimension(n,2)::point
    real(8),dimension(n-1,4)::func_xy

    real(8),dimension(n)::length    !----节点累加弦长
    real(8),dimension(n,2)::XT,YT   !----坐标的参数方程插值节点
    real(8),dimension(n-1,4)::func_xt,func_yt  !---分段系数


    real(8)::a1,a2,a3,a4,b1,b2,b3,b4

    real(8)::x1,y1
    real(8)::dx,dy,ddx,ddy   !---x,y对弦长的导数

    integer(4)::i,j,k,kk
    real(8)::s,t,s1,s2,t1,t2


    point(:,1)=point1(:,1);    !---插值节点x
    point(:,2)=point2(:,1);    !---插值节点y

    !-----节点x,y三次样条插值
    call spline3(n,point,ds1,dsn,func_xy )

    !----开始根据横坐标插值
    do i=1,NumN
        x1=InterX(i)

        !----确定区间
        do j=1,n-1
          s1=point(j,1);
          s2=point(j+1,1);

          if( x1<point(1,1) ) then
            k=1
            exit
          elseif( x1>point(n,1) ) then
            k=n-1
            exit
          elseif( (x1-s1)*(x1-s2)<=0.0 ) then
            k=j
            exit
          end if

        end do

        !-------插值
        a1=func_xy(k,1);   a2=func_xy(k,2);
        a3=func_xy(k,3);   a4=func_xy(k,4);

        Ns(i)=a1*x1**3.0+a2*x1**2.0+a3*x1+a4;

        Nds(i)=3.0*a1*x1**2.0+2.0*a2*x1+a3;

        Ndds(i)=6.0*a1*x1+2.0*a2

    end do

    return
end subroutine trispline








