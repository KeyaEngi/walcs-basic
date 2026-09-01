!输入：smtf（光滑函数）；NumT（一个完整的耐波性时间步中，需要计算的砰击次数）；time0,time1（砰击计算区间首尾时间）    
!Dis0,Dis1（区间起始两个时刻，船体摇荡运动的位移）    
!输出：t0SlamForce（区间起始时刻对应的整体砰击载荷）；MainForce（一个迭代步长里的砰击平均载荷）；
    

    
subroutine SlamCase3( smtf,NumT,time0,time1,Dis0,Dis1,MEsecLoad0,MEsecLoad1,t0SlamForce,MainForce )
    use ShipHullVar,only:Cog
    use Slamming
    use ArrayOperations
    use Constant,only:pi,U0,rou,NR


    implicit none
    real(8)::smtf
    integer(4)::NumT
    real(8)::time0     !---t0时刻
    real(8)::time1     !---t1时刻
    real(8),dimension(6)::Dis0   !---t0时刻船体摇荡运动位移
    !20230427
    real(8),dimension(6,SlamNumLine)::MDis0     !考虑梁振动的各砰击剖面t0时刻摇荡运动
    real(8),dimension(6)::Dis1   !---t1时刻船体摇荡运动位移
    real(8),dimension(NR)::t0SlamForce       !---t0时刻船体整体砰击载荷
    real(8),dimension(NR)::MainForce         !----一个耐波性时间步长里，平均的砰击载荷
    
    
    !与梁振动运动有关的变量
    real(8),dimension( 1:6,1:SlamNumLine )::MEsecLoad0,MEsecLoad1
    
    
    real(8),dimension(6)::DtDis  !---t0到t1时间段船体摇荡运动速度
    !20230427
    real(8),dimension(6,SlamNumLine)::MDtDis      !考虑梁振动各砰击剖面速度
    real(8)::SlamDeltT        !---时间步长
    real(8)::time             !---当前时刻
    
    real(8),dimension(6)::Dis    !---time时刻摇荡运动位移
    real(8),dimension(6)::Vel    !---time时刻摇荡运动速度
    
    !20230427
    real(8),dimension(6,SlamNumLine)::MDis     !---time时刻摇荡运动位移(考虑梁振动)  
    real(8),dimension(6,SlamNumLine)::MVel     !---time时刻摇荡运动速度(考虑梁振动)
    
    
    !******进行砰击事件开始与结束的速度判据阈值
    real(8)::temSlamHV
    real(8)::temSlamCV
    
    integer(4),dimension(SlamNumLine,2)::temSlamRelaP   !---上一时刻与本时刻关键点与波面高度关系
    integer(4),dimension(SlamNumLine)::temSlamCase      !---砰击事件发生与否准则  0未发生/1发生
    real(8),dimension(SlamNumLine)::temIniPenetration   !---砰击初始入侵距离
    real(8),dimension(SlamNumLine)::temIniRise          !---砰击初始入侵距离对应的液面抬升高度
    real(8),dimension(SlamNumLine,3)::temP0VerticalWP   !---与第一典型砰击节点同x，y位置处的波面节点坐标
    
    real(8),dimension(3)::temp0,temp1
    real(8)::temp0feta                   !---与temp0对应的入射波面升高
    real(8)::temfeta                   !---与temp0对应的入射波面升高
    real(8)::temAngle,temAngle1          !剖线倾角
    real(8),dimension(6)::temDis,temVel   !---船体运动临时存储数组
    !20230427
    real(8),dimension(6,SlamNumLine)::MtemDis,MtemVel
    
    real(8),dimension(3)::InteractionP    !---过剖线平面与中纵剖面的交线，与波面交点
    real(8),dimension(3)::InteractionRVel  !----交点处的船波相对运动速度
    real(8),dimension(3)::temSvel,temWvel    !---交点处船与波面的临时绝对速度
    real(8),dimension(3)::temNor            !剖线单位方向向量
    
    !******确定计算砰击压力所需的相对运动距离，速度，加速度
    real(8)::timeDt
    real(8),dimension(7)::temRelaL
    real(8),dimension(5)::temRelaL2
    real(8),dimension(3)::temRelaV
    
    real(8)::InpRelaL    !---作为压力计算输入的侵入距离，速度，加速度
    real(8)::InpRelaV
    real(8)::InpRelaA
    real(8)::InpIniRelaL
    real(8)::InpIniRise
    real(8)::DeltRait           !---计算相对运动
     
    
    integer(4)::temNumIntC
    real(8),dimension(SlamNumIntC )::temNC
    real(8),dimension(SlamNumIntC )::temRiseC
    
    !*****************接收的压力及节点信息
    integer(4)::temSlamNumN
    real(8),dimension(0:SlamNumIntegration,1:2)::temSPrePoint
    real(8),dimension(0:SlamNumIntegration)::temSpre
    
    
    
    !*****************接收的剖线砰击载荷
    real(8),dimension(SlamNumLine,NR )::temSlamForce
    real(8),dimension(NumT,NR)::TimeSlamForce       !----区间中砰击时刻对应的整体砰击力
    
    real(8),dimension(SlamNumLine,3)::SlineForce    !----单个砰击剖线上的砰击压力积分(在平动坐标系下的三个方向的力)
    
    !*****************20231115瞬时位置
    real(8),allocatable,dimension(:,:)::temSlamShipIniRMP,temSlamWaveIniRMP
    
    !20240116波面起伏
    real(8),allocatable,dimension(:)::Waveupdown
    
    
    real(8),dimension(3)::tempx
    
    integer(4)::i,j,k,it,kk,iii,ii,jj
    
    real(8)::s,t,s1,t1,s2,t2,s3,t3
    
    real(8)::Kexi
    real(8)::riteT,riteS
    
    !***以上是定义变量***
    temNumIntC=SlamNumIntC          !----半剖面插值Ct节点数量100
    
    !**给定砰击发生、结束判定**
    temSlamHV=0.03     !----砰击发生速度判据
    temSlamCV=0.03     !----砰击结束速度判据
    
    DeltRait=1.0       !相对运动？？
    
    !**计算砰击起始到结束时间段内的船体摇荡运动速度、时间步长**
    !20230427
    DtDis(:)=(Dis1(:)-Dis0(:) )/(time1-time0 )    !---速度
    do i=1,SlamNumLine
        do j=1,6
           MDtDis(j,i)=DtDis(j)+( MEsecLoad1(j,i)-MEsecLoad0(j,i) )/(time1-time0 )         
        end do
    end do

    SlamDeltT=(time1-time0 )/real( NumT )         !---时间步长
    
    
    temSlamCase(:)=SlamCase(:)       !砰击发生与否
    temIniPenetration(:)=SlamIniPenetration(:)    !初始浸入距离
    temIniRise(:)=SlamIniRise(:)     !初始浸入距离对应的压面抬升高度
    !20230427
    do i=1,SlamNumLine
        do j=1,6
            MDis0(j,i)=Dis0(j)+MEsecLoad0(j,i)        
        end do
    end do
    
    !**先计算time0时刻的波面与第一典型砰击节点相对位置关系**
    do i=1,SlamNumLine
        temp0=SlamIniType(i,1,1:3)      !剖面的第一典型节点
        !----计算剖面第一典型节点的瞬时(time)位置(平动坐标系下)，这里计算的是temp0时刻
        temp1(1:3)=temp0(1:3)-Cog(1:3);
        !20230427
        temp0(1:3)=VectorL2G(temp1(1:3),MDis0(4:6,i) )   !---因旋转运动产生的线位移
        temp0(1:3)=temp0(1:3)+Cog(1:3)+MDis0(1:3,i)
        
        !----计算与瞬时典型节点同x,y位置处的入射波面高度        
        !smtf:光滑函数；temp0：第一典型砰击节点的瞬时位置；time0：砰击起始时刻        
        !temp0feta:与temp0对应的入射波面高度，与temp0对应的波浪的瞬时绝对速度
        
        call Slam_WaveElevation( smtf,temp0,time0,temp0feta,temWvel )
        
        if( temp0feta>=temp0(3) )then     !---先判断time0时刻与波面的相对位置
            temSlamRelaP(i,1)=1           !上一时刻与本时刻关键点与波面高度关系（1代表波面高于典型砰击节点）  
        else
            temSlamRelaP(i,1)=0            
        end if
     
        
    end do
    
    
    !**开始计算**
    TimeSlamForce=0.0;       !区间中砰击时刻对应的整体砰击力
    
    do it=1,NumT              !对每一个耐波性时间步内的砰击步进数进行循环
        
        SlineForce=0.0        !每时刻计算前先初始化；单个砰击剖线上的砰击压力积分(在平动坐标系下的三个方向的力)
        
        time=time0+SlamDeltT*real(it)    !砰击实际时刻
        
        if(it==NumT) time=time1
        
        !Dis=Dis0+DtDis*(time-time0 )   !time时刻，船体摇荡运动的位移
        !Vel=DtDis                      !time时刻，船体摇荡运动的速度
        !20230427
        do i=1,SlamNumLine
            do j=1,6
                MDis(j,i)=MDis0(j,i)+MDtDis(j,i)*(time-time0 )
                MVel(j,i)=MDtDis(j,i)                
            end do
        end do
        
        
        !**先判断time（瞬时）时刻典型节点与波面高度关系**
        temP0VerticalWP=0.0;       !与第一典型砰击节点同x，y位置处的波面节点坐标
        do i=1,SlamNumLine            
            
            !先统一计算典型节点与波面高度的关系
            temp0=SlamIniType(i,1,1:3)   !剖面的第一典型节点
            !计算剖面第一典型节点的瞬时(time)位置(平动坐标系下)
            temp1(1:3)=temp0(1:3)-Cog(1:3);
            !20230427
            temp0(1:3)=VectorL2G(temp1(1:3),MDis(4:6,i) )
            temp0(1:3)=temp0(1:3)+Cog(1:3)+MDis(1:3,i)
          
            !计算与瞬时典型节点同x,y位置处的入射波面高度、波面的绝对速度
            call Slam_WaveElevation(smtf,temp0,time,temp0feta,temWvel  )
          
            if( temp0feta>=temp0(3) ) then  !---先判断time时刻与波面的相对位置
            temSlamRelaP(i,2)=1           !1代表波面高于典型砰击点
            else
            temSlamRelaP(i,2)=0
            end if
            
            !-----还需要保留与temp0同x,y的波面交点位置
            temP0VerticalWP(i,1:2)=temp0(1:2)
            temP0VerticalWP(i,3)=temp0feta       !与第一典型砰击节点同x，y位置处的波面节点坐标
          
        end do
        
        !**开始判别各剖线是否发生砰击并计算各剖线上的砰击力**
        temSlamForce=0.0;
        
        do i=1,SlamNumLine    !对剖线进行循环(判断剖线砰击发生与否，并计算剖线砰击力)
        
        !先判断剖线砰击发生与否    
        if( temSlamCase(i)==0 )then
            !**砰击事件发生第一判据：由水面以上入侵到水面以下**
            if( temSlamRelaP(i,1)==0.and.temSlamRelaP(i,2)==1 )then
                
                !与波面的交点默认为典型位置瞬时节点位置
                InteractionP(1:3)=temP0VerticalWP(i,1:3)    !与第一典型砰击节点同x，y位置处的波面节点坐标
                
                !**先计算合成速度**
                temSvel=0.0;  temWvel=0.0;      !交点处船体瞬时绝对速度、波面瞬时绝对速度
                
                !**求交点处的波面运动绝对速度**                    
                !smtf:光滑函数；InteractionP：与第一典型砰击节点同x，y位置处的波面节点坐标
                !time:砰击实际时刻
                
                !s:临时变量（无意义）;temWvel:交点处波面运动绝对速度                    
                call Slam_WaveElevation(smtf,InteractionP,time,s,temWvel  )
                
                !**求交点处船体瞬时绝对速度**
                !20230427
                temp1(1:3)=InteractionP(1:3)-Cog(1:3)-MDis(1:3,i)  !瞬时位置(关于重心)
                temSvel(1:3)=R_RCross(MVel(4:6,i),temp1(1:3) )     !因旋转运动产生的线速度
                temSvel(1:3)=temSvel(1:3)+MVel(1:3,i)              !关键点线速度
                
                !**这里如果考虑弹性的话是不是要考虑梁振动速度（先不考虑20230419）**
                !**开始进行速度合成**
                InteractionRVel=0.0;
                !用固结于船的点速度 减去 波面运动速度
                InteractionRVel(1:3)=temSvel(1:3)-temWvel(1:3)
                !考虑定常航速
                InteractionRVel(1)=InteractionRVel(1)+U0
                
                !**计算考虑纵倾影响的剖线瞬时倾角**
                temAngle1=SlamAngle(i)
                !20230427
                temAngle=temAngle1-MDis(5,i)    !剖线实时倾角
                
                !**构建倾斜线单位向量，注意正负号**
                if( abs(temAngle-Pi/2.0)<=1.0e-4 )then    !---竖直线
                    temNor=0.0;
                    temNor(3)=-1.0;
                else
                    temp0(3)=InteractionP(3)-1.0
                    temp0(2)=InteractionP(2)
                    temp0(1)=InteractionP(1)-1.0/tan(temAngle)

                    temNor(1:3)=temp0(1:3)-InteractionP(1:3)
                    !---标准化
                    s=sqrt( (temNor(1))**2.0+(temNor(2))**2.0+(temNor(3))**2.0 )
                    temNor(1:3)=temNor(1:3)/s
                                        
                end if
                
                !**计算合成速度在给定方向上的投影速度**
                t=dot_product( InteractionRVel(1:3),temNor(1:3) )
                
                !**根据速度判断砰击发生与否、并计算砰击初始浸入距离、砰击初始入侵距离对应的液面抬升高度**
                if( t>=temSlamHV )then
                    temSlamCase(i)=1;
                    temIniPenetration(i)=0.0;     !砰击初始浸入距离
                    temIniRise(i)=0.0;            !砰击初始入侵距离对应的液面抬升高度
                    
                    SlamCase(i)=1;
                    SlamIniPenetration(i)=0.0;
                    SlamIniRise(i)=0.0;
                        
                end if
                    
     
            
            elseif( temSlamRelaP(i,2)==1 )then       !上一砰击时刻第一典型砰击节点在水面以下，下一砰击时刻在水面以下
                
                temAngle1=SlamAngle(i)    !先预估一个夹角
                
                !**先计算time时刻的典型节点实时位置**
                temp0=SlamIniType(i,1,1:3)   !剖面的第一典型节点
                !计算剖面第一典型节点的瞬时(time)位置(平动坐标系下)
                temp0(1:3)=temp0(1:3)-Cog(1:3);
                !20230427
                temp0(1:3)=VectorL2G(temp0(1:3),MDis(4:6,i) )
                temp0(1:3)=temp0(1:3)+Cog(1:3)+MDis(1:3,i)
                
                temAngle=temAngle1-MDis(5,i)    !----剖线实时倾角(为什么要两步隔开)
                
                !**计算与瞬时典型节点同x,y位置处的入射波面高度**
                call Slam_WaveElevation(smtf,temp0,time,temp0feta,temWvel  )
                
                !**找过剖线平面与中纵剖面的交线，与波面交点**
                if( abs(temAngle-Pi/2.0)<=1.0e-4 )then    !---竖直线
                    InteractionP(1:2)=temp0(1:2);
                    InteractionP(3)=temp0feta;
                                    
                else    !---倾斜线
                    !----先计算斜线斜率和与纵轴的交点
                    s=tan(temAngle )
                    t=temp0(3)-s*temp0(1)
                    
                    !----根据斜线方程，逐步迭代寻找交点
                    k=0;  kk=0;
                    
                    s1=temp0(1);  t1=temp0(3);   !----起始点
                    
                    do
                        k=k+1
                        t2=temp0(3)+0.5*real(k)  !---先从temp0(3)往上累加，找出区间点
                        s2=(t2-t)/s
                        
                        tempx(1)=s2;  tempx(2)=temp0(2);  tempx(3)=t2;
                        call Slam_WaveElevation(smtf,tempx,time,temfeta,temWvel  )
                        
                        if( temfeta<=t2 )then      !---二分法
                            iii=0
                            do
                                iii=iii+1;
                                if(iii>35) then
                                    InteractionP(1:3)=tempx(1:3)
                                    kk=1
                                    exit                                       
                                end if
                                                               
                                s3=(s1+s2)/2.0;
                                t3=(t1+t2)/2.0;
                                
                                tempx(1)=s3;  tempx(2)=temp0(2);  tempx(3)=t3;
                                call Slam_WaveElevation(smtf,tempx,time,temfeta,temWvel  )
                                
                                if( abs(temfeta-t3 )<=1.0e-6 ) then
                                    InteractionP(1:3)=tempx(1:3)
                                    kk=1
                                    exit
                                end if

                                if( temfeta>t3 ) then
                                    s1=s3;
                                    t1=t3;
                                elseif(temfeta<t3 ) then
                                    s2=s3;
                                    t2=t3;
                                end if
                                          
                            end do
                        end if
                        if(kk==1) exit
                        
                        s1=s2;
                        t1=t2;
                   
                    end do
                end if
                
                !**计算合成速度**
                
                temSvel=0.0;  temWvel=0.0;
                
                !**求交点处的波面运动绝对速度**
                call Slam_WaveElevation(smtf,InteractionP,time,s,temWvel  )
                
                !**求交点处船体瞬时绝对速度**
                !20230427
                temp1(1:3)=InteractionP(1:3)-Cog(1:3)-MDis(1:3,i)  !---瞬时位置(关于重心)
                temSvel(1:3)=R_RCross(MVel(4:6,i),temp1(1:3) )  !---因旋转运动产生的线速度
                temSvel(1:3)=temSvel(1:3)+MVel(1:3,i)   !---关键点线速度
                
                !**开始进行速度合成**
                InteractionRVel=0.0;
                !----用固结于船的点速度 减去 波面运动速度
                InteractionRVel(1:3)=temSvel(1:3)-temWvel(1:3)
                !----考虑定常航速
                InteractionRVel(1)=InteractionRVel(1)+U0
                
                !**构建倾斜线单位向量，注意正负号**
                if( abs(temAngle-Pi/2.0)<=1.0e-4 ) then   !---竖直线
                    temNor=0.0;
                    temNor(3)=-1.0;
                else
                    temNor(1:3)=temp0(1:3)-InteractionP(1:3)
                    !---标准化
                    s=sqrt( (temNor(1))**2.0+(temNor(2))**2.0+(temNor(3))**2.0 )
                    temNor(1:3)=temNor(1:3)/s
                end if
                
                
                !**计算合成速度在给定方向上的投影速度**
                t=dot_product(InteractionRVel(1:3),temNor(1:3)  )
                
                !**判断砰击发生与否、并计算砰击初始浸入距离、砰击初始入侵距离对应的液面抬升高度**
                
                if( t>=temSlamHV.and.iii<=35 )then
                    temSlamCase(i)=1;
                    temIniPenetration(i)=sqrt( (temp0(1)-InteractionP(1) )**2.0+(temp0(2)-InteractionP(2) )**2.0+(temp0(3)-InteractionP(3) )**2.0 )
                    
                    SlamCase(i)=1;
                    SlamIniPenetration(i)=temIniPenetration(i);     !砰击初始入侵距离
                    
                    !**根据初始入侵距离插值对应的液面抬升高度**
                    temNC(:)=SlamNC(i,:)         !----侵入距离插值节点
                    temRiseC(:)=SlamRiseC(i,:)   !----液面抬升插值节点
                    
                    s3=temIniPenetration(i);
                    do ii=1,temNumIntC
                        
                        if(ii==1) then
                            s1=0.0;    !---入侵距离
                            t1=0.0;    !---液面抬升
                        else
                            s1=temNC(ii-1);
                            t1=temRiseC(ii-1);
                        end if
                        
                        s2=temNC(ii);
                        t2=temRiseC(ii);
                        
                        if( s3>=temNC(temNumIntC) ) then
                            temIniRise(i)=temRiseC(temNumIntC)
                            exit

                        elseif(s3<=0.0  ) then
                            temIniRise(i)=0.0
                            exit
                        end if
                        
                        if( (s3-s1)*(s3-s2)<=0.0 ) then
                            Kexi=(s3-s1)/(s2-s1)
                            temIniRise(i)=t1*(1.0-Kexi)+t2*Kexi
                            exit
                        end if
            
                    end do
                    
                    !-----插值完毕
                    SlamIniRise(i)=temIniRise(i)      !初始入侵距离对应的压面抬升高度
                          
                end if
                

            end if
                
            
        elseif( temSlamCase(i)==1 )then
            
            !-----当temSlamCase(i)==1时，认为从时段的起始时刻已经在发生砰击了
            !-----此时 应当判断本时刻是否砰击结束
                
            !----注意，既然已经发生砰击，就证明SlamAngle(i)已经给定了
            !----      因此，此时temAngle1不是乱给的
            !----      当前统一默认为Pi/2.0
            temAngle1=SlamAngle(i)
            
            !**先计算time时刻的典型节点实时位置**
            temp0=SlamIniType(i,1,1:3)   !---剖面的第一典型节点
            !----计算剖面第一典型节点的瞬时(time)位置(平动坐标系下)
            temp0(1:3)=temp0(1:3)-Cog(1:3);
            !20230427
            temp0(1:3)=VectorL2G(temp0(1:3),MDis(4:6,i) )
            temp0(1:3)=temp0(1:3)+Cog(1:3)+MDis(1:3,i)
            
            temAngle=temAngle1-MDis(5,i)    !----剖线实时倾角
            
            !**计算与瞬时典型节点同x,y位置处的入射波面高度**
            call Slam_WaveElevation(smtf,temp0,time,temp0feta,temWvel  )
            
            !**要确保水面一定在船关键点上面**
            
            if( temp0feta>temp0(3) )then
                
                !**先找交点**
                if( abs(temAngle-Pi/2.0)<=1.0e-4 ) then     !---竖直线
                    InteractionP(1:2)=temp0(1:2);
                    InteractionP(3)=temp0feta;
                    
                else       !---倾斜线
                    !----先计算斜线斜率和与纵轴的交点
                    s=tan(temAngle )
                    t=temp0(3)-s*temp0(1)

                    !----根据斜线方程，逐步迭代寻找交点
                    k=0;  kk=0;
                    s1=temp0(1);  t1=temp0(3);   !----起始点
                    do
                        k=k+1
                        t2=temp0(3)+0.5*real(k)  !---先从temp0(3)往上累加，找出区间点
                        s2=(t2-t)/s

                        tempx(1)=s2;  tempx(2)=temp0(2);  tempx(3)=t2;
                        call Slam_WaveElevation(smtf,tempx,time,temfeta,temWvel  )
                        
                        if( temfeta<=t2 )then     !---二分法
                            iii=0;
                            do
                                iii=iii+1;
                                if(iii>35) then
                                    InteractionP(1:3)=tempx(1:3)
                                    kk=1
                                    exit                                       
                                end if
                                
                                s3=(s1+s2)/2.0;
                                t3=(t1+t2)/2.0;
                                
                                tempx(1)=s3;  tempx(2)=temp0(2);  tempx(3)=t3;
                                call Slam_WaveElevation(smtf,tempx,time,temfeta,temWvel  )
                                
                                if( abs(temfeta-t3 )<=1.0e-6 ) then
                                    InteractionP(1:3)=tempx(1:3)
                                    kk=1
                                    exit
                                end if
                                
                                if( temfeta>t3 ) then
                                    s1=s3;
                                    t1=t3;
                                elseif(temfeta<t3 ) then
                                    s2=s3;
                                    t2=t3;
                                end if
                                
                            end do
                        end if
                    
                        if(kk==1) exit

                        s1=s2;
                        t1=t2;
                    
                    end do                    
                end if         !找交点
                
                !**计算交点处的相对速度**
                
                temSvel=0.0;  temWvel=0.0;
                
                !**交点处的波面运动绝对速度**
                call Slam_WaveElevation(smtf,InteractionP,time,s,temWvel  )
                
                !**交点处船体瞬时绝对速度**
                !20230427
                temp1(1:3)=InteractionP(1:3)-Cog(1:3)-MDis(1:3,i)  !---瞬时位置(关于重心)
                temSvel(1:3)=R_RCross(MVel(4:6,i),temp1(1:3) )  !---因旋转运动产生的线速度
                temSvel(1:3)=temSvel(1:3)+MVel(1:3,i)   !---关键点线速度
                
                !**开始进行速度合成**
                InteractionRVel=0.0;
                !----用固结于船的点速度 减去 波面运动速度
                InteractionRVel(1:3)=temSvel(1:3)-temWvel(1:3)
                !----考虑定常航速
                InteractionRVel(1)=InteractionRVel(1)+U0
                
                !**构建倾斜线单位向量，注意正负号**
                if( abs(temAngle-Pi/2.0)<=1.0e-4 ) then   !---竖直线
                    temNor=0.0;
                    temNor(3)=-1.0;
                else
                    temNor(1:3)=temp0(1:3)-InteractionP(1:3)
                    !---标准化
                    s=sqrt( (temNor(1))**2.0+(temNor(2))**2.0+(temNor(3))**2.0 )
                    temNor(1:3)=temNor(1:3)/s
                end if
                
                !**计算合成速度在给定方向上的投影速度**
                t=dot_product(InteractionRVel(1:3),temNor(1:3)  )
                
                !**根据速度判断砰击发生与否、并计算砰击初始浸入距离、砰击初始入侵距离对应的液面抬升高度**
                if(t<temSlamCV.or.iii>35) then
                    temSlamCase(i)=0;
                    temIniPenetration(i)=0.0;
                    temIniRise(i)=0.0;

                    SlamCase(i)=0;
                    SlamIniPenetration(i)=0.0;  
                    SlamIniRise(i)=0.0;
                end if  
                
                
            else
                temSlamCase(i)=0;
                temIniPenetration(i)=0.0;
                temIniRise(i)=0.0;

                SlamCase(i)=0;
                SlamIniPenetration(i)=0.0;  
                SlamIniRise(i)=0.0;
                
            end if

        end if              !先判断剖线砰击发生与否
        
        !**至此第i条剖线发生砰击的情况一判断完毕**temSlamCase=0or1
        !**根据本时刻曲线砰击发生情况,开始计算砰击压力**
        !**加两个限制条件,砰击剖线生成成功且砰击发生**
        
        !**开始计算剖线上的砰击力**
        if( temSlamCase(i)==1.and.SlamLineCase(i)==1 )then
            
            temRelaL=0.0;   
            temRelaV=0.0;
            
            do j=1,7
                
                if(j==1) riteT=-3.0
                if(j==2) riteT=-2.0
                if(j==3) riteT=-1.0
                if(j==4) riteT=0.0
                if(j==5) riteT=1.0
                if(j==6) riteT=2.0
                if(j==7) riteT=3.0
                
                !构造time前后七个时刻
                timeDt=time+riteT*SlamDeltT*DeltRait     !临时时间存储变量
                !----临时船体重心位移，速度存储
                !temDis=Dis+DtDis*riteT*SlamDeltT*DeltRait     !---timeDt时刻位移
                !temVel=DtDis                                  !---timeDt时刻速度
                !20230427
                do ii=1,SlamNumLine
                    do jj=1,6
                        MtemDis(jj,ii)=MDis(jj,ii)+MDtDis(jj,ii)*riteT*SlamDeltT*DeltRait
                        MtemVel(jj,ii)=MDtDis(jj,ii)
                    end do
                    
                end do
                
                !20230427
                temAngle=SlamAngle(i)-MtemDis(5,i)    !计算t时刻剖线实时倾角
                
                !**计算剖线第一典型点实时位置**
                temp0=SlamIniType(i,1,1:3)   !---剖面的第一典型节点
                !----计算剖面第一典型节点的瞬时位置(平动坐标系下)
                temp0(1:3)=temp0(1:3)-Cog(1:3);
                temp0(1:3)=VectorL2G(temp0(1:3),MtemDis(4:6,i) )
                temp0(1:3)=temp0(1:3)+Cog(1:3)+MtemDis(1:3,i)
                
                !**先判断典型节点与波面高度关系
                call Slam_WaveElevation(smtf,temp0,timeDt,temfeta,temWvel  )
                if( temp0(3)<=temfeta  ) then
                    riteS=1.0
                else
                    riteS=-1.0
                end if
                
                !**寻找倾斜线与波面的交点**
                if( abs(temAngle-Pi/2.0)<=1.0e-4 )then            !---竖直线          
                    InteractionP(1:2)=temp0(1:2);
                    InteractionP(3)=temfeta;                    
                else
                    !----先计算斜线斜率和与纵轴的交点
                    s=tan(temAngle )        !斜率
                    t=temp0(3)-s*temp0(1)    !z轴截距

                    !----根据斜线方程，逐步迭代寻找交点
                    k=0;  kk=0;

                    s1=temp0(1);  t1=temp0(3);   !----起始点
                    do
                        k=k+1
                        t2=temp0(3)+0.5*real(k)*riteS   !---先从temp0(3)往上累加，找出区间点
                        s2=(t2-t)/s

                        tempx(1)=s2;  tempx(2)=temp0(2);  tempx(3)=t2;
                        !计算tempx（x，y）对应位置处的波面高度及波面的绝对速度
                        call Slam_WaveElevation(smtf,tempx,timeDt,temfeta,temWvel  )
                        
                        if( riteS==1.0 .and. temfeta<=t2 )then   !二分法
                            do
                                s3=(s1+s2)/2.0;
                                t3=(t1+t2)/2.0;

                                tempx(1)=s3;  tempx(2)=temp0(2);  tempx(3)=t3;
                                call Slam_WaveElevation(smtf,tempx,timeDt,temfeta,temWvel  )
                                
                                if( abs(temfeta-t3 )<=1.0e-6 ) then
                                    InteractionP(1:3)=tempx(1:3)
                                    kk=1
                                    exit
                                end if
                                
                                if( temfeta>t3 ) then
                                    s1=s3;
                                    t1=t3;
                                elseif(temfeta<t3 ) then
                                    s2=s3;
                                    t2=t3;
                                end if
                            
                            end do
                        elseif( riteS==-1.0 .and. temfeta>=t2 )then
                            do
                                s3=(s1+s2)/2.0;
                                t3=(t1+t2)/2.0;

                                tempx(1)=s3;  tempx(2)=temp0(2);  tempx(3)=t3;
                                call Slam_WaveElevation(smtf,tempx,timeDt,temfeta,temWvel  )
                                
                                if( abs(temfeta-t3 )<=1.0e-6 ) then
                                    InteractionP(1:3)=tempx(1:3)
                                    kk=1
                                    exit
                                end if
                                
                                if( temfeta>t3 ) then
                                    s2=s3;
                                    t2=t3;
                                elseif(temfeta<t3 ) then
                                    s1=s3;
                                    t1=t3;
                                end if

                            end do
                        end if
                        
                        if(kk==1) exit
                    
                        s1=s2;
                        t1=t2;
                                       
                    end do
                end if
                
                !***********至此，与波面的交点InteractionP寻找完毕
                !-----计算相对运动距离
                s=sqrt( (InteractionP(1)-temp0(1) )**2.0+(InteractionP(2)-temp0(2) )**2.0+(InteractionP(3)-temp0(3) )**2.0 )
                                
                !---注意  正的表示在波面以下
                !         负的表示在波面以上
                temRelaL(j)=s*riteS          !-----计算不同时段的入侵距离   
                
            end do
            
            !**尝试,光顺**
            temRelaL2=0.0;
            temRelaL2(1:5)=temRelaL(2:6);
            
            !**瞬时的相对入侵深度、相对运动速度、相对运动加速度**
            InpRelaL=temRelaL2(3);
            InpRelaV=(-temRelaL2(5)+8.0*temRelaL2(4)-8.0*temRelaL2(2)+temRelaL2(1) )/(12.0*SlamDeltT*DeltRait )
            InpRelaA=( temRelaL2(4)-2.0*temRelaL2(3)+temRelaL2(2) )/(SlamDeltT*DeltRait )**2.0
            
            !**加限制条件**
            if( InpRelaL<=0.01 ) then
                InpRelaV=0.0;
                InpRelaA=0.0;
            end if
            if( InpRelaV<temSlamHV ) then
                InpRelaV=0.0;
                InpRelaA=0.0;
            end if
            
            InpIniRelaL=temIniPenetration(i)   !---本次砰击发生时的初始入侵距离
            InpIniRise=temIniRise(i)           !---初始入侵距离对应的压面抬升
            
            if( InpRelaL<InpIniRelaL ) then
                InpRelaL=InpIniRelaL
            end if
            
            
            !**计算选定曲线的砰击压力**
            temSlamNumN=0;
            temSPrePoint=0.0;
            temSpre=0.0;
            
            !i:对砰击剖线进行循环的变量；InpIniRelaL：初始入侵距离；InpIniRise：初始入侵距离对应的压面抬升
            !InpRelaL：相对入侵深度；InpRelaV：相对运动速度；InpRelaA：相对运动加速度
              
            !temSlamNumN：计算压力节点数量；temSPrePoint：计算压力节点坐标；temSpre：计算节点压力
            
            
            call Slam_Pressure( i,InpIniRelaL,InpIniRise,InpRelaL,InpRelaV,InpRelaA,temSlamNumN,temSPrePoint,temSpre)
            
            
            !**砰击压力在曲线上的积分(空间积分),并且在此处，已经乘以了密度和带宽**
            !Dis:船体瞬时摇荡位移；i对砰击剖线循环的变量；temSlamNumN：计算压力节点数
            !temSPrePoint：计算压力节点坐标；temSpre：计算点压力
              
            !temSlamForce：节点压力积分得到的剖线9个方向压力；SlineForce单个砰击剖线上的砰击压力积分(在平动坐标系下的三个方向的力)
            
            temSlamForce(i,:)=0.0;
            
            call Slam_PreSpacialIntegral(MDis(:,i),i,temSlamNumN,temSPrePoint(0:temSlamNumN,1:2),temSpre(0:temSlamNumN),temSlamForce(i,:),SlineForce(i,:)  )
            
            
        else
            temSlamForce(i,:)=0.0;    !当前时刻每条剖线上的砰击力
              
            SlineForce(i,:)=0.0;      !单个砰击剖线上的砰击压力积分(在平动坐标系下的三个方向的力)
            
       
        end if
        
        TimeSlamForce(it,1:NR)=TimeSlamForce(it,1:NR)+temSlamForce(i,1:NR)     !time时刻砰击整体载荷叠加
   
            
        end do                !对剖线进行循环(判断剖线砰击发生与否，并计算剖线砰击力)
        
        !**更新位置关系**
        temSlamRelaP(:,1)=temSlamRelaP(:,2)
        
        
        !20231116修改，输出单根砰击剖线的砰击载荷
        write(65,"(f12.3,1x,\)",advance='no') time
        write(66,"(f12.3,1x,\)",advance='no') time          
          
        do i=1,SlamNumLine
            if( SlamLineCase(i)==0 ) then
                write(65,"(e15.4,1x,\)",advance='no') 0.0
                write(66,"(e15.4,1x,\)",advance='no') 0.0                  
            else
                write(65,"(e15.4,1x,\)",advance='no') SlineForce(i,1)/1000.0
                write(66,"(e15.4,1x,\)",advance='no') SlineForce(i,3)/1000.0
                  
            end if
        end do
        write(65,"(/)")
        write(66,"(/)")
        
        

        
        
    end do                    !对每一个耐波性时间步内的砰击步进数进行循环
    
    
    !**计算平均砰击载荷(随船平动下，用于计算运动响应)**
    MainForce=0.0;
    do j=1,NR
        do it=1,NumT
            if(it==1) then
                t1=t0SlamForce(j)
                t2=TimeSlamForce(it,j)
            else
                t1=TimeSlamForce(it-1,j)
                t2=TimeSlamForce(it,j)
            end if
                
            MainForce(j)=MainForce(j)+(t1+t2)/2.0*SlamDeltT    !---砰击力的冲量
  
        end do
    end do
        
    !**计算一段时间内的平均砰击力**
    MainForce=MainForce/(time1-time0 )
        
    !**更新下一时间步需要用到的t0时刻的整体砰击载荷**
    t0SlamForce(1:NR)=TimeSlamForce(NumT,1:NR)
        
    !do i=1,SlamNumLine
    !    if( SlamCase(i)==1 )then
    !        write(*,"(A20,1x,i8)") 'Slamming line',i
    ! 
    !    end if
    !
    !end do
    
    !20231115修改，添加船中纵剖线各点与波面的相对位置计算
    if( trim(adjustl(SlamRelaMctr))=='YES' )then
        allocate( temSlamShipIniRMP(SlamNumShipRMP,3),temSlamWaveIniRMP(SlamNumWaveRMP,3) )
        temSlamShipIniRMP=0.0;
        temSlamWaveIniRMP=0.0;
        
        !20240116波面起伏
        allocate( Waveupdown(2) )
        Waveupdown=0.0;
        
        !-----先计算实时位置
        Dis(1:6)=Dis1(1:6)    !----船体摇荡运动位移
        time=time1
        
        do i=1,SlamNumShipRMP      !---确定实船上的点坐标
            !----先确定平动坐标系下的实时位置
            temp0(1:3)=SlamShipRMP(i,1:3)      !---平动坐标系下的节点坐标(船体)
            !----计算节点瞬时位移
            temp1(1:3)=temp0(1:3)-Cog(1:3);
            temp0(1:3)=VectorL2G(temp1(1:3),Dis(4:6) )   !---角位移得到的节点坐标
            temp0(1:3)=temp0(1:3)-temp1(1:3)   !---角运动确定的线偏移量
            
            temp0(2)=temp0(2)+Dis(2);   !---补上船体本身的线位移(去除纵荡)      
            temp0(3)=temp0(3)+Dis(3); 
              
            !----确定瞬时位置(去除纵荡)
            temSlamShipIniRMP(i,1:3)=SlamShipIniRMP(i,1:3)+temp0(1:3)

        end do
        
        do i=1,SlamNumWaveRMP
            !----先确定波面实时位置
            temp0(1:3)=SlamWaveRMP(i,1:3)      !---平动坐标系下的节点坐标(静水面)
            
            !---瞬时位置
            temp1(1:3)=temp0(1:3)-Cog(1:3);
            temp0(1:3)=VectorL2G(temp1(1:3),Dis(4:6) )   !---角运动得到的节点线位移            
            
            temp1(1:3)=temp0(1:3)-temp1(1:3)   !---角运动确定的线偏移量
            temp1(2)=temp1(2)+Dis(2)    !---考虑船体的水平方向的线位移
            
            !----瞬时位置
            temp0(1:3)=temp0(1:3)+Cog(1:3)+Dis(1:3)   
            !---计算波面起伏
            !----计算与瞬时节点同x,y位置处的入射波面高度
            
            call Slam_WaveElevation(smtf,temp0,time,temp0feta,temWvel  )
            
            temp1(3)=temp0feta   !---垂向上的偏移量
            
            !-----确定波面上的点的瞬时位置(去除纵荡)
            temSlamWaveIniRMP(i,1:3)=SlamWaveIniRMP(i,1:3)+temp1(1:3)
 
        end do
        
        
        !-----输出(船体纵剖线上的点)
        write(67,"(f12.3,1x,\)",advance='no') time
        write(68,"(f12.3,1x,\)",advance='no') time          
          
        do i=1,SlamNumShipRMP
            write(67,"(e15.4,1x,\)",advance='no') temSlamShipIniRMP(i,1)
            write(68,"(e15.4,1x,\)",advance='no') temSlamShipIniRMP(i,3)
        end do
        write(67,"(/)")
        write(68,"(/)") 
        
        !-----输出(静水面上的点)
        write(69,"(f12.3,1x,\)",advance='no') time
        write(70,"(f12.3,1x,\)",advance='no') time          
          
        do i=1,SlamNumWaveRMP
            write(69,"(e15.4,1x,\)",advance='no') temSlamWaveIniRMP(i,1)
            write(70,"(e15.4,1x,\)",advance='no') temSlamWaveIniRMP(i,3)
        end do
        write(69,"(/)")
        write(70,"(/)")
        
        
        Waveupdown(1)=temSlamWaveIniRMP(3,1)
        Waveupdown(2)=temSlamWaveIniRMP(3,3)
        write(71,"(f12.3,1x,\)",advance='no') time
        do i=1,2
            write(71,"(e15.4,1x,\)",advance='no') Waveupdown(i)
        end do
        write(71,"(/)")
        
        
    
    end if
    
    
    

    
end subroutine SlamCase3