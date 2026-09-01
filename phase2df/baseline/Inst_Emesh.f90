



subroutine Inst_Emesh(Nwh,Nl,Node0,InstUr,InstDur)

    use ShipHullVar,only:StruSect,StruXN,drm,StruZN,StruYN,StruZSC
    use Constant,only:NR
    
    implicit none

    integer(4)::Nwh,Nl
    real(8),dimension(Nwh,Nl,3)::Node0
    real(8),dimension(Nwh,Nl,NR,3)::InstUr
    real(8),dimension(Nwh,Nl,NR,3,3)::InstDur

    real(8),dimension(Nwh*Nl,1:3)::xav
    real(8),dimension(Nwh*Nl,1:NR,1:6)::Ns,Nds,Ndds   !---插值节点位移，位移对x一阶导数，位移对x二阶导数
    integer(4)::numN

    real(8),dimension(Nwh,Nl)::NSectZN,NSectYN,NSectZsc

    real(8),dimension(6)::drm1,drmN    !---6中变形的,节点边界条件

    real(8)::x1,x2,x3,y1,y2,y3,yy1,yy2,yy3,yyy1,yyy2,yyy3,z1,z2,z3
    real(8)::temZn,temYn,temZsc
    real(8)::temU,temV,temW,tema,temb,temc   !----节点对应横剖面的位移转角
    real(8)::temdU,temdV,temdW,temda,temdb,temdc
    real(8)::kexi
    integer(4)::i,j,k,ii,jj,kk
    
!*************************************
    numN=0;   xav=0.0;
    do i=1,Nwh
    do j=1,Nl
        numN=numN+1;
        xav(numN,1:3)=Node0(i,j,1:3);  !---构造待插值点
    end do
    end do


    Ns=0.0; Nds=0.0; Ndds=0.0;
    drm1=0.0; drmN=0.0;
    !------
    do i=7,NR
        !-----先给边界条件
        drm1(1)=(drm(i,2,1)-drm(i,1,1))/(StruXN(2)-StruXN(1));            
        drmN(1)=(drm(i,StruSect,1)-drm(i,StruSect-1,1))/(StruXN(StruSect)-StruXN(StruSect-1));   

        drm1(2)=(drm(i,2,2)-drm(i,1,2))/(StruXN(2)-StruXN(1)); 
        drmN(2)=(drm(i,StruSect,2)-drm(i,StruSect-1,2))/(StruXN(StruSect)-StruXN(StruSect-1));

        drm1(3)=(drm(i,2,3)-drm(i,1,3))/(StruXN(2)-StruXN(1));   
        drmN(3)=(drm(i,StruSect,3)-drm(i,StruSect-1,3))/(StruXN(StruSect)-StruXN(StruSect-1));


        drm1(4)=(drm(i,2,4)-drm(i,1,4))/(StruXN(2)-StruXN(1))
        drmN(4)=(drm(i,StruSect,4)-drm(i,StruSect-1,4))/(StruXN(StruSect)-StruXN(StruSect-1))
        drm1(5)=(drm(i,2,5)-drm(i,1,5))/(StruXN(2)-StruXN(1))
        drmN(5)=(drm(i,StruSect,5)-drm(i,StruSect-1,5))/(StruXN(StruSect)-StruXN(StruSect-1))
        drm1(6)=(drm(i,2,6)-drm(i,1,6))/(StruXN(2)-StruXN(1))
        drmN(6)=(drm(i,StruSect,6)-drm(i,StruSect-1,6))/(StruXN(StruSect)-StruXN(StruSect-1))

        !----开始对节点进行插值
        !----三次样条插值得到待插节点处的 数值，对x一阶导数，对x二阶导数
        do j=1,6
            call espl1(StruXN(:),drm(i,:,j),StruSect,drm1(j),drmN(j),xav(:,1),numN,Ns(:,i,j),Nds(:,i,j),Ndds(:,i,j) )
        end do

    end do

    !-----计算节点对应横剖面的 中和轴高度，剪心高度(考虑到结构的突变性，不能用样条插值)
    !-----分段线性插值
    NSectZN=0.0; NSectYN=0.0; NSectZsc=0.0;

    do i=1,Nwh
    do j=1,Nl
        x3=Node0(i,j,1);

        if(x3<StruXN(1) ) then
            x1=StruXN(1);     x2=StruXN(2)  !---结构剖面x坐标
            y1=StruZN(1);     y2=StruZN(2)  !---结构剖面水平中和轴高度Zn
            yy1=StruYN(1);    yy2=StruYN(2) !---结构剖面纵向中和轴Yn
            yyy1=StruZSC(1);  yyy2=StruZSC(2)  !---结构剖面扭转中心Zsc
        elseif(x3>StruXN(StruSect) ) then
            x1=StruXN(StruSect-1);     x2=StruXN(StruSect)  !---结构剖面x坐标
            y1=StruZN(StruSect-1);     y2=StruZN(StruSect)  !---结构剖面水平中和轴高度Zn
            yy1=StruYN(StruSect-1);    yy2=StruYN(StruSect) !---结构剖面纵向中和轴Yn
            yyy1=StruZSC(StruSect-1);  yyy2=StruZSC(StruSect)  !---结构剖面扭转中心Zsc            
        else
            do k=1,StruSect-1
                x1=StruXN(k);  x2=StruXN(k+1);
                if( (x1-x3)*(x2-x3)<=0.0 ) exit
            end do
            y1=StruZN(k);     y2=StruZN(k+1)  !---结构剖面水平中和轴高度Zn
            yy1=StruYN(k);    yy2=StruYN(k+1) !---结构剖面纵向中和轴Yn
            yyy1=StruZSC(k);  yyy2=StruZSC(k+1)  !---结构剖面扭转中心Zsc
        end if

        kexi=(x3-x1)/(x2-x1)
        NSectZN(i,j)=y1*(1.0-kexi)+y2*kexi;
        NSectYN(i,j)=yy1*(1.0-kexi)+yy2*kexi;
        NSectZsc(i,j)=yyy1*(1.0-kexi)+yyy2*kexi;
    end do
    end do

!**************************************
!---------开始计算湿表面网格节点弹性体模态位移，及梯度
!-----主义坐标系   默认此时计算得到的值在随船平动坐标系下
    
    do ii=7,Nr
        
        k=0;
        do i=1,Nwh
        do j=1,Nl
            !----节点坐标
            x1=Node0(i,j,1);    y1=Node0(i,j,2);    z1=Node0(i,j,3)
            !----节点对应横剖面 Zn,Yn,Zsc
            temZn=NSectZN(i,j); temYn=NSectYN(i,j); temZsc=NSectZsc(i,j)

            k=k+1;
            !----节点对应横剖面的位移转角
            temU=Ns(k,ii,1); temV=Ns(k,ii,2); temW=Ns(k,ii,3);
            tema=Ns(k,ii,4); temb=Ns(k,ii,5); temc=Ns(k,ii,6);
            !----节点对应横剖面的位移及转角的x的偏导
            temdU=Nds(k,ii,1); temdV=Nds(k,ii,2); temdW=Nds(k,ii,3);
            temda=Nds(k,ii,4); temdb=Nds(k,ii,5); temdc=Nds(k,ii,6);

            !-----弹性模态节点位移
            InstUr(i,j,ii,1)=temU+(z1-temZn)*temb-(y1-temYn)*temc
            InstUr(i,j,ii,2)=temV-(z1-temZsc)*tema
            InstUr(i,j,ii,3)=temW+(y1-temYn)*tema

            !-----弹性模态节点各方向位移的梯度
            InstDur(i,j,ii,1,1)=temdU+(z1-temZn)*temdb-(y1-temYn)*temdc
            InstDur(i,j,ii,1,2)=-temc
            InstDur(i,j,ii,1,3)=temb

            InstDur(i,j,ii,2,1)=temdV-(z1-temZsc)*temda
            InstDur(i,j,ii,2,2)=0.0
            InstDur(i,j,ii,2,3)=-tema

            InstDur(i,j,ii,3,1)=temdW+(y1-temYn)*temda
            InstDur(i,j,ii,3,2)=tema
            InstDur(i,j,ii,3,3)=0.0

        end do
        end do

    end do

!******************开始处理刚体运动模态

    do i=1,Nwh
    do j=1,Nl
        !----节点坐标(关于重心)
        x1=Node0(i,j,1);    y1=Node0(i,j,2);    z1=Node0(i,j,3)            
        
        !------6个刚体运动模态下节点的位移 
        InstUr(i,j,1,1)=1.0;
        InstUr(i,j,2,2)=1.0;
        InstUr(i,j,3,3)=1.0;
        InstUr(i,j,4,2)=-z1;    InstUr(i,j,4,3)=y1;
        InstUr(i,j,5,1)=z1;     InstUr(i,j,5,3)=-x1;
        InstUr(i,j,6,1)=-y1;    InstUr(i,j,6,2)=x1;

        !----节点位移梯度
        InstDur(i,j,4,2,3)=-1.0;   InstDur(i,j,4,3,2)=1.0;
        InstDur(i,j,5,1,3)=1.0;    InstDur(i,j,5,3,1)=-1.0;
        InstDur(i,j,6,1,2)=-1.0;   InstDur(i,j,6,2,1)=1.0;

    end do
    end do



    return
end subroutine Inst_Emesh