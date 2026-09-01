

!----子程序目的：给bdf文件中单元中心点计算弹性位移
subroutine bdfEmesh(temNum_ele,xav0,bdfele_Ur )

    use ShipHullVar,only:Cog,StruSect,StruXN,drm,StruZN,StruYN,StruZSC
!!    use PanelGeometry,only:bdfnum_ele,bdfele_apply,bdfele_xav,bdfele_Ur
    use Constant,only:NR

    implicit none

    integer(4)::temNum_ele
    real(8),dimension( temNum_ele,3)::xav0
    real(8),dimension(temNum_ele,NR,3 )::bdfele_Ur  !---新划分单元中心点各模态下的位移

    real(8),dimension( temNum_ele,3)::temxav
    real(8),dimension( temNum_ele,1:NR,1:6)::Ns,Nds,Ndds   !---插值节点位移，位移对x一阶导数，位移对x二阶导数
    integer(4)::numN
    real(8),dimension(temNum_ele)::NSectZN,NSectYN,NSectZsc

    real(8),dimension(6)::drm1,drmN    !---6中变形的,节点边界条件

    real(8)::x1,x2,x3,y1,y2,y3,yy1,yy2,yy3,yyy1,yyy2,yyy3,z1,z2,z3
    real(8)::temZn,temYn,temZsc
    real(8)::temU,temV,temW,tema,temb,temc   !----节点对应横剖面的位移转角
    real(8)::temdU,temdV,temdW,temda,temdb,temdc
    real(8)::kexi
    integer(4)::i,j,k,ii,jj,kk


    !注：待插值节点temxav原本就是随船平动坐标系下关于重心的相对位置

    numN=temNum_ele;
    temxav(1:NumN,1:3)=xav0(1:temNum_ele,1:3)


    !-----先插值振型
    Ns=0.0; Nds=0.0; Ndds=0.0;
    drm1=0.0; drmN=0.0;
    !------
    do i=7,NR
        !-----先给边界条件 (边界处的一阶导数)
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

!------更换样条插值方法
        do j=1,6
            call trispline( StruSect,StruXN(:),drm(i,:,j),drm1(j),drmN(j),numN,temxav(:,1),Ns(:,i,j),Nds(:,i,j),Ndds(:,i,j) )
        end do

    end do


!****************************
    !-----计算节点对应横剖面的 中和轴高度，剪心高度(考虑到结构的突变性，不能用样条插值)
    !-----分段线性插值
    NSectZN=0.0; NSectYN=0.0; NSectZsc=0.0;

    do i=1,numN
        x3=temxav(i,1);   !----节点对应x坐标

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
        NSectZN(i)=y1*(1.0-kexi)+y2*kexi;
        NSectYN(i)=yy1*(1.0-kexi)+yy2*kexi;
        NSectZsc(i)=yyy1*(1.0-kexi)+yyy2*kexi;
    end do


    !**********************
    !----开始计算节点处弹性模态位移
    !**********************
!-----主坐标系   默认此时计算得到的值在随船平动坐标系下
    
    do ii=7,Nr
        
        k=0;
        do i=1,numN

            !----节点坐标
            x1=temxav(i,1);    y1=temxav(i,2);    z1=temxav(i,3);
            !----节点对应横剖面 Zn,Yn,Zsc
            temZn=NSectZN(i); temYn=NSectYN(i); temZsc=NSectZsc(i)

            k=k+1;
            !----节点对应横剖面的位移转角
            temU=Ns(k,ii,1); temV=Ns(k,ii,2); temW=Ns(k,ii,3);
            tema=Ns(k,ii,4); temb=Ns(k,ii,5); temc=Ns(k,ii,6);
!!            !----节点对应横剖面的位移及转角的x的偏导
!!            temdU=Nds(k,ii,1); temdV=Nds(k,ii,2); temdW=Nds(k,ii,3);
!!            temda=Nds(k,ii,4); temdb=Nds(k,ii,5); temdc=Nds(k,ii,6);

            !-----弹性模态节点位移
            bdfele_Ur(i,ii,1)=temU+(z1-temZn)*temb-(y1-temYn)*temc
            bdfele_Ur(i,ii,2)=temV-(z1-temZsc)*tema
            bdfele_Ur(i,ii,3)=temW+(y1-temYn)*tema

        end do
    end do

    do i=1,numN
        !----节点坐标(关于重心)
        !----节点坐标
        x1=temxav(i,1);    y1=temxav(i,2);    z1=temxav(i,3);           
        
        !------6个刚体运动模态下节点的位移 
        bdfele_Ur(i,1,1)=1.0;
        bdfele_Ur(i,2,2)=1.0;
        bdfele_Ur(i,3,3)=1.0;
        bdfele_Ur(i,4,2)=-z1;    bdfele_Ur(i,4,3)=y1;
        bdfele_Ur(i,5,1)=z1;     bdfele_Ur(i,5,3)=-x1;
        bdfele_Ur(i,6,1)=-y1;    bdfele_Ur(i,6,2)=x1;

    end do


!!    !----检查节点位移
!!    open(unit=10002,file='CTRelacUr.txt')
!!    open(unit=10001,file='elacUr.txt')
!!        write(10002,"(i8,i8)") numN
!!
!!!        do ii=1,NR
!!!        write(10001,"( A,i8)") '模态',ii
!!        ii=7;
!!        do i=1,numN
!!            if(bdfele_apply(i)==0) cycle
!!
!!            if(abs(temxav(i,2) )>=1.0  ) cycle
!!
!!            write(10001,"(3(f15.6,1x),3(f16.5,1x) )") (temxav(i,k),k=1,3),(bdfele_Ur(i,ii,k),k=1,3)
!!
!!        end do
!!!        write(10001,"(/)")
!!!        end do
!!    close(10001)
!!    close(10002)





    return
end subroutine bdfEmesh





