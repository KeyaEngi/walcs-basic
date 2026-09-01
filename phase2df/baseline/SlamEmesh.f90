!输入：Num0剖线上的细化节点数；temPoint0随船平动坐标系下关于重心的各剖线砰击节点坐标
!Cog:重心坐标；NR：计算维度（刚体+弹性）；StruSect（结构剖面数，迁移矩阵法）  
!drm:梁模型的各个弹性模态在各横剖面的位移;StruZN:结构剖面水平中和轴高度
!StruYN:结构剖面中和轴纵向位置；StruZSC：结构剖面扭转中心
!输出：Sur，SRur：左右舷剖线节点模态位移； SDur，SRDur：左右舷剖线节点模态位移梯度   
    
    
!--------砰击剖线计算节点关于弹性体的转化
!----注：此模块只针对单独的一条剖线进行处理
!----    传递到此模块中的节点坐标已经是关于重心的位置了

subroutine SlamEmesh(Num0,temPoint0,Sur,SDur,SRur,SRDur)
  
  use ShipHullVar,only:Cog,StruSect,StruXN,drm,StruZN,StruYN,StruZSC
  use Constant,only:NR

  implicit none

  integer(4)::Num0
  real(8),dimension(Num0,3)::temPoint0   !---节点坐标，关于重心
  !----左舷
  real(8),dimension(Num0,NR,3)::Sur      !---各模态下，节点位移
  real(8),dimension(Num0,NR,3,3)::SDur   !---各模态下，节点位移梯度
  !----右舷
  real(8),dimension(Num0,NR,3)::SRur      !---各模态下，节点位移
  real(8),dimension(Num0,NR,3,3)::SRDur   !---各模态下，节点位移梯度

  real(8),dimension(Num0)::Spx    !---节点x坐标
  real(8),dimension(Num0,1:NR,1:6)::Ns,Nds,Ndds   !---插值得到的位移，对x一阶导，对x二阶导
  real(8),dimension(6)::drm1,drmN    !---6种变形的,节点边界条件

  real(8),dimension(Num0)::NSectZN,NSectYN,NSectZsc   !---插值得到节点对应中和轴位置(Z,Y)，剪心位置

  real(8)::x1,x2,x3,y1,y2,y3,yy1,yy2,yy3,yyy1,yyy2,yyy3,z1,z2,z3
  real(8)::temZn,temYn,temZsc
  real(8)::temU,temV,temW,tema,temb,temc   !----节点对应横剖面的位移转角
  real(8)::temdU,temdV,temdW,temda,temdb,temdc
  real(8)::kexi
  integer(4)::i,j,k,ii,jj,kk
  real(8)::s,t

  !*************************

  Ns=0.0; Nds=0.0; Ndds=0.0;
  drm1=0.0; drmN=0.0;
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
        call trispline( StruSect,StruXN(:),drm(i,:,j),drm1(j),drmN(j),Num0,temPoint0(:,1),Ns(:,i,j),Nds(:,i,j),Ndds(:,i,j) )
    end do

  end do

    !-----计算节点对应横剖面的 中和轴高度，剪心高度(考虑到结构的突变性，不能用样条插值)
    !-----分段线性插值
    NSectZN=0.0; NSectYN=0.0;  NSectZsc=0.0;

    do i=1,Num0
        x3=temPoint0(i,1)

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
        
        kexi=(x3-x1)/(x2-x1)   !---比例系数

        NSectZN(i)=y1*(1.0-kexi)+y2*kexi;
        NSectYN(i)=yy1*(1.0-kexi)+yy2*kexi;
        NSectZsc(i)=yyy1*(1.0-kexi)+yyy2*kexi;
    end do

    !**************************************
    !-----------开始计算节点对应模态位移
    Sur=0.0;   SDur=0.0;
    SRur=0.0;  SRDur=0.0;
    do ii=7,NR
        !----左舷
        do i=1,Num0
            !-----节点坐标
            x1=temPoint0(i,1);  y1=temPoint0(i,2);  z1=temPoint0(i,3); 
            !-----节点对应横剖面特性 Zn,Yn,Zsc
            temZn=NSectZN(i); temYn=NSectYN(i); temZsc=NSectZsc(i);
            !----节点对应横剖面的位移转角
            temU=Ns(i,ii,1); temV=Ns(i,ii,2); temW=Ns(i,ii,3);
            tema=Ns(i,ii,4); temb=Ns(i,ii,5); temc=Ns(i,ii,6);
            !----节点对应横剖面的位移及转角的x的偏导
            temdU=Nds(i,ii,1); temdV=Nds(i,ii,2); temdW=Nds(i,ii,3);
            temda=Nds(i,ii,4); temdb=Nds(i,ii,5); temdc=Nds(i,ii,6);

            !-----弹性体模态节点位移
            SUr(i,ii,1)=temU+(z1-temZn)*temb-(y1-temYn)*temc
            SUr(i,ii,2)=temV-(z1-temZsc)*tema
            SUr(i,ii,3)=temW+(y1-temYn)*tema

            !-----弹性模态节点各方向位移的梯度
            Sdur(i,ii,1,1)=temdU+(z1-temZn)*temdb-(y1-temYn)*temdc
            Sdur(i,ii,1,2)=-temc
            Sdur(i,ii,1,3)=temb

            Sdur(i,ii,2,1)=temdV-(z1-temZsc)*temda
            Sdur(i,ii,2,2)=0.0
            Sdur(i,ii,2,3)=-tema

            Sdur(i,ii,3,1)=temdW+(y1-temYn)*temda
            Sdur(i,ii,3,2)=tema
            Sdur(i,ii,3,3)=0.0
        end do

        !----右舷
        do i=1,Num0
            !-----节点坐标
            x1=temPoint0(i,1);  y1=-temPoint0(i,2);  z1=temPoint0(i,3); 
            !-----节点对应横剖面特性 Zn,Yn,Zsc
            temZn=NSectZN(i); temYn=NSectYN(i); temZsc=NSectZsc(i);
            !----节点对应横剖面的位移转角
            temU=Ns(i,ii,1); temV=Ns(i,ii,2); temW=Ns(i,ii,3);
            tema=Ns(i,ii,4); temb=Ns(i,ii,5); temc=Ns(i,ii,6);
            !----节点对应横剖面的位移及转角的x的偏导
            temdU=Nds(i,ii,1); temdV=Nds(i,ii,2); temdW=Nds(i,ii,3);
            temda=Nds(i,ii,4); temdb=Nds(i,ii,5); temdc=Nds(i,ii,6);

            !-----弹性体模态节点位移
            SRUr(i,ii,1)=temU+(z1-temZn)*temb-(y1-temYn)*temc
            SRUr(i,ii,2)=temV-(z1-temZsc)*tema
            SRUr(i,ii,3)=temW+(y1-temYn)*tema

            !-----弹性模态节点各方向位移的梯度
            SRdur(i,ii,1,1)=temdU+(z1-temZn)*temdb-(y1-temYn)*temdc
            SRdur(i,ii,1,2)=-temc
            SRdur(i,ii,1,3)=temb

            SRdur(i,ii,2,1)=temdV-(z1-temZsc)*temda
            SRdur(i,ii,2,2)=0.0
            SRdur(i,ii,2,3)=-tema

            SRdur(i,ii,3,1)=temdW+(y1-temYn)*temda
            SRdur(i,ii,3,2)=tema
            SRdur(i,ii,3,3)=0.0
        end do

    end do

    !-------刚体运动模态位移
    do i=1,Num0
        !-----节点坐标
        x1=temPoint0(i,1);  y1=temPoint0(i,2);  z1=temPoint0(i,3);       

        !------6个刚体运动模态下节点的位移 
        SUr(i,1,1)=1.0;
        SUr(i,2,2)=1.0;
        SUr(i,3,3)=1.0;
        SUr(i,4,2)=-z1;    SUr(i,4,3)=y1;
        SUr(i,5,1)=z1;     SUr(i,5,3)=-x1;
        SUr(i,6,1)=-y1;    SUr(i,6,2)=x1;

        !----节点位移梯度
        Sdur(i,4,2,3)=-1.0;   Sdur(i,4,3,2)=1.0;
        Sdur(i,5,1,3)=1.0;    Sdur(i,5,3,1)=-1.0;
        Sdur(i,6,1,2)=-1.0;   Sdur(i,6,2,1)=1.0;

    end do

    !-----右舷
    !-------刚体运动模态位移
    do i=1,Num0
        !-----节点坐标
        x1=temPoint0(i,1);  y1=-temPoint0(i,2);  z1=temPoint0(i,3);       

        !------6个刚体运动模态下节点的位移 
        SRUr(i,1,1)=1.0;
        SRUr(i,2,2)=1.0;
        SRUr(i,3,3)=1.0;
        SRUr(i,4,2)=-z1;    SRUr(i,4,3)=y1;
        SRUr(i,5,1)=z1;     SRUr(i,5,3)=-x1;
        SRUr(i,6,1)=-y1;    SRUr(i,6,2)=x1;

        !----节点位移梯度
        SRdur(i,4,2,3)=-1.0;   SRdur(i,4,3,2)=1.0;
        SRdur(i,5,1,3)=1.0;    SRdur(i,5,3,1)=-1.0;
        SRdur(i,6,1,2)=-1.0;   SRdur(i,6,2,1)=1.0;

    end do


  return
end subroutine SlamEmesh