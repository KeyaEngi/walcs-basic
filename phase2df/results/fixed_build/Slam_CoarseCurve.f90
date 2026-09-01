!输入：SecID典型砰击节点序号；TypeAngle考虑船体纵倾影响的砰击角；TemNumSectP2等分弧长节点数
!输出：temSlamCase：砰击剖线生成情况；TemSectP2：
    
!功能：根据砰击典型节点坐标，以及典型砰击节点对应的砰击角，在原始水平剖线上进行提取获取该位置处的半横剖线节点坐标，进而等分弧长获取
    !该位置处粗等分弧长半横剖线节点坐标
    
    
!**********************************
!------根据典型节点处的夹角提取砰击剖线粗节点


subroutine Slam_CoarseCurve( SecID,temSlamCase,TypeAngle,TemNumSectP2,TemSectP2 )
  
  use Slamming,only:SlamLibNumZ,SlamLibNumPort,SlamLibPortNode,SlamLibType
  !                水平剖线层数；每层半剖线节点数（左舷）；每层半剖线节点坐标（左舷）；用户坐标系下，关于重心的典型砰击节点坐标
  !这四个全局变量目前全部已知
  use Constant,only:Pi
  implicit none

  integer(4)::SecID     !----该截取的典型节点ID（或者砰击二维剖线ID）
  integer(4)::temSlamCase   !-----砰击状态文件
  real(8)::TypeAngle    !----剖面倾角(用户坐标系下)
  integer(4)::TemNumSectP   !----提取的砰击剖面节点个数
  real(8),dimension( SlamLibNumZ+2,3 )::TemSectP   !----提取的剖面节点坐标
  
  integer(4)::TemNumSectP2   !----优化后提取的砰击剖面节点个数
  real(8),dimension( TemNumSectP2,3 )::TemSectP2

  integer(4)::TemNumP
  real(8),allocatable,dimension( :,: )::TemP
  real(8)::TemZ,TemZ1
  real(8)::TemX
  real(8),allocatable,dimension(:)::TemPY

  real(8),allocatable,dimension(:)::length



  real(8),dimension(3)::x1,x2
  integer(4)::i,j,k,ii,jj,kk
  real(8)::Kexi
  real(8)::s,t,s1,s2


      TemNumSectP=0;      !提取的砰击剖面节点个数（带有纵倾角度的横剖面，注意角度90认为是垂直的横剖面）
      TemSectP=0.0;       !提取的砰击剖面节点坐标

      !------起始节点坐标
      TemNumSectP=1;
      TemSectP(TemNumSectP,1:3)=SlamLibType(SecID,1:3 )   !---起始节点坐标（认为对应位置的典型砰击节点为起始节点）
      TemZ=TemSectP(TemNumSectP,3)        !----起始点高度


      !----从下往上依次截剖线寻找节点
      do i=1,SlamLibNumZ     !对水平剖线层数进行循环
          TemZ1=SlamLibPortNode(i,1,3)   !----第i层的高度

          if( TemZ1>TemZ ) then   !---依次向上寻找节点
              
              !-----先确定Temz1处对应的纵向坐标temX
              if( abs( TypeAngle-Pi/2.0  )<=1.0e-4 ) then   !---竖直
                  TemX=SlamLibType(SecID,1  )     !Temz1处对应的纵向坐标     
              else   !----有较大的倾角
                  s=tan(TypeAngle  );   !--斜率
                  t=SlamLibType(SecID,3  )-s*SlamLibType(SecID,1  );   !---与纵轴交点

                  TemX=( TemZ1-t )/s
              end if
              
              !-----先构建临时点组
              TemNumP=SlamLibNumPort(i);        !每层水平型线原始节点数
              allocate(TemP(TemNumP,3 )  )      !每层水平型线原始节点坐标
              TemP(1:TemNumP,1:3)=SlamLibPortNode(i,1:TemNumP,1:3)

              allocate( TemPY(TemNumP ) )
              TemPY=0.0;
              !------开始寻找交点
              kk=0;
              TemPY=0.0;    !----待插y值
              do j=1,TemNumP-1
                  x1(:)=TemP(j,:);
                  x2(:)=Temp(j+1,:);

                  if( ( x1(1)-TemX )*(x2(1)-TemX )<=1.0e-10 ) then     !这里是不是有问题呀
                      kk=kk+1;
                      if( abs(x2(1)-x1(1))<=1.0e-6 ) then
                          TemPY(kk)=x1(2);
                      else
                          Kexi=(TemX-x1(1) )/(x2(1)-x1(1))
                          TemPY(kk)=x1(2)*(1.0-Kexi)+x2(2)*Kexi
                      end if
                  end if
              end do

              if(kk==0) then      !事实上也有可能不能插值找到节点y值

              elseif(kk==1) then
                TemNumSectP=TemNumSectP+1;
                TemSectP(TemNumSectP,1)=TemX;
                TemSectP(TemNumSectP,2)=TemPY(1);
                TemSectP(TemNumSectP,3)=TemZ1;
              else
                s=TemPY(1)
                do j=2,kk
                  s=max(s,TemPY(j))     !去较大的
                end do
                TemNumSectP=TemNumSectP+1;
                TemSectP(TemNumSectP,1)=TemX;
                TemSectP(TemNumSectP,2)=s;
                TemSectP(TemNumSectP,3)=TemZ1;
              end if

              deallocate( TemP,TemPY )
          end if
      end do

      !************判断是否截取到点（根据是否取到点判断砰击剖线的生成情况）
      temSlamCase=0;
      if( TemNumSectP==1 ) then
          temSlamCase=0;
      elseif( TemNumSectP>=2 ) then
          temSlamCase=1;
      end if

      !-----------------将截取的节点进行等分弧长处理

      allocate( length(TemNumSectP ) )      !砰击剖线（横半剖线）长度
      length=0.0;

      s=0.0;
      do i=2,TemNumSectP
          x1(:)=TemSectP(i-1,:);
          x2(:)=TemSectP(i,:);

          t=sqrt( (x1(1)-x2(1))**2.0+(x1(2)-x2(2))**2.0+(x1(3)-x2(3))**2.0 )
          s=s+t;

          length(i)=s;
      end do

      !-----开始等分节点
      TemSectP2=0.0;
      TemSectP2(1,:)=TemSectP(1,:);      !等分弧长后的第1个节点坐标（半横剖线最底部）
      TemSectP2(TemNumSectP2,:)=TemSectP(TemNumSectP,:);       !等分弧长后的第30个节点坐标（半横剖线最顶部）

      s=length(TemNumSectP )/real(TemNumSectP2-1 )
      !插值获得等分弧长后的第2-29个节点坐标
      do i=2,TemNumSectP2-1
          t=s*real(i-1);

          do j=1,TemNumSectP-1
              s1=length(j);
              s2=length(j+1);

              if( (t-s1)*(t-s2)<=0.0 ) then
                  if(abs(s2-s1)<=1.0e-8) then
                      TemSectP2(i,:)=TemSectP(j,:);
                      exit
                  else
                      x1(:)=TemSectP(j,:);
                      x2(:)=TemSectP(j+1,:);
                      Kexi=(t-s1)/(s2-s1);
                      TemSectP2(i,:)=x1(:)*(1.0-Kexi)+x2(:)*Kexi
                      exit
                  end if
              end if
          end do
      end do







      !**********至此 砰击剖线粗节点构造完毕

  return
end subroutine Slam_CoarseCurve



