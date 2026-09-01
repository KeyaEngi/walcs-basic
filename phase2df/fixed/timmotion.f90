subroutine Timmotion
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ShipHullVar
    use Constant
    use PanelGeometry
    use Slamming
    use ArrayOperations
    !20231113修改
    use IrreWaveVar
    
    use lift
    use SAILPARAM_MOD, only: DP, SAIL_OK, SAIL_ERR_LOW_WIND_SPEED
    use SAILPLANNING_ADAPTER_MOD, only: ComputeSailPlanningDryRun, &
        TransformSailLoadToPlanning
    
    implicit none
    
    integer(4)::i,j,k,it,ii,jj
    

    
    real(8)::wavelscale       !----特征波长船长比
    
    integer(4)::load_ID
    character(len=4)::load_name
    
    real(8),allocatable,dimension(:)::smooth
    real(8)::smtf
    
    real(8),allocatable,dimension(:,:)::AA_mat,BB_mat,CC_mat,invAA_mat     !质量矩阵、阻尼、刚度、质量矩阵逆矩阵
    
    real(8),dimension(2*NR)::bb,zz,y,dy     !模态主坐标和一阶导，主坐标对时间的一阶导和二阶导(更新运动)
    real(8),dimension(2*NR)::Secty,Sectdy      !与y,dy对应（用于更新载荷）
    real*8,dimension(4)::a                    !龙格库塔时间步
    real(8)::tt,t

    ! Phase-1 SailModule dry-run inputs. Wind points along actual air motion
    ! in planning/global axes. R_SAIL_BODY is CE minus CG in body axes (m).
    real(DP), parameter :: SAIL_TEST_WIND_GLOBAL(3) = [0.0_DP, 10.0_DP, 0.0_DP]
    real(DP), parameter :: SAIL_TEST_R_BODY(3) = [0.0_DP, 0.0_DP, 5.0_DP]
    real(DP), parameter :: SAIL_TEST_DELTA_DEG = 0.0_DP
    logical, parameter :: SAIL_ENABLE = .TRUE.
    real(DP), parameter :: SAIL_COUPLING_TOL = 1.0E-10_DP
    real(DP) :: sail_v_cg_body(3), sail_omega_body(3), sail_v_wind_body(3)
    real(DP) :: sail_vrel_h_mag, sail_alpha_db_deg, sail_cl, sail_cd
    real(DP) :: sail_force_body(3), sail_moment_body(3)
    real(DP) :: sail_load_body(6), sail_load_planning(6)
    real(DP), dimension(NR) :: sail_force_generalized, sail_formot_before
    real(DP), dimension(NR) :: sail_formot_delta, sail_acc_before, sail_acc_after
    real(DP), dimension(NR) :: sail_acc_expected
    real(DP) :: sail_check_scale
    integer :: sail_transform_ierr
    character(len=2048) :: sail_transform_message
    integer :: sail_ierr
    character(len=2048) :: sail_message
    logical :: sail_print
    integer, parameter :: SAIL_HISTORY_UNIT = 114
    real(DP) :: sail_v_rel_body(3), sail_e_upstream_body(3)
    real(DP) :: sail_c_chord_body(3), sail_alpha_raw_deg, sail_q_dynamic
    
    !与mlm理论计算水动升力相关的变量
    real(8),dimension(6)::t0ScogM,t1ScogM        !区间前后两个时刻，船体摇荡运动的位移
    real(8)::Stime0,Stime1                       !砰击计算区间首尾时间
    real(8),dimension(NR)::SlamForce             !一个耐波性时间步长里的平均砰击载荷
    integer(4)::SlamNumT                         !一个完整的耐波性时间步里需要计算的砰击步数
    real(8),dimension(NR)::t0SlamForce           !区间前端时刻对应的砰击整体载荷
    
    
    
    
    !入射势相关变量
    real(8),dimension(3)::tmpx
    real(8)::tmpt
    real*8,dimension(Nwh,NL,3)::derPhiI
    real*8,dimension(Nwh,NL)::DtPhiI
    real(8),dimension(Nwh,NL,3)::derDtPhiI
    
    !与瞬时入射力、静水压力相关的变量
    real(8),dimension(NR)::Inst_ForceIS      !第一种方式非线性入射力
    
    
    real(8),dimension(NR)::Instbdf_ForceIS     !第二种方式非线性入射力
    
    
    !与梁振动运动有关的变量
    real(8),dimension( 1:6,1:SlamNumLine )::MEsecLoad0,MEsecLoad1
    
    

    
    
    !****前面为变量定义****
    
    !**定义方程右端合力项**
    allocate( ForMot(NR),ForceI(NR) )
    ForMot=0.0;    ForceI=0.0
    
    !**定义梁各剖面弹性振动**
    allocate( MEsecLoad(6,SlamNumLine) )
    
    
    !**机翼升力**
    allocate( air_force(NR) )
    air_force=0.0;
    
    
    !**计算网格法向
    
    
    
    
    !**给定系数矩阵**
    allocate( AA_mat(NR,NR),BB_mat(NR,NR),CC_mat(NR,NR),invAA_mat(NR,NR) )
    AA_mat=0.0;  BB_mat=0.0;   CC_mat=0.0;   invAA_mat=0.0
    
    AA_mat(1:NR,1:NR)=StruM(1:NR,1:NR)
    BB_mat(1:NR,1:NR)=StruB(1:NR,1:NR)
    
    !临界升沉、纵摇阻尼系数修正20230519
    BB_mat(3,3)=BcoefH*2.0*sqrt( StruM(3,3)*HERM(3,3) )
    BB_mat(5,5)=BcoefP*2.0*sqrt( StruM(5,5)*HERM(5,5) )
    
    
    !CC_mat(1:NR,1:NR)=StruC(1:NR,1:NR)+HERM(1:NR,1:NR)
    

    !20230519舍去静水恢复力矩阵    
    if( trim( adjustl(Non_Linear))=='LT' )then
        CC_mat(1:NR,1:NR)=StruC(1:NR,1:NR)+HERM(1:NR,1:NR)        
    else        
        CC_mat(1:NR,1:NR)=StruC(1:NR,1:NR)        
    end if
    
    
    !**计算质量矩阵逆矩阵**
    call Gauss_Jordan( NR,NR,AA_mat,invAA_mat )
    
    
   
        
        
    do load_ID=1,loadnum      !对工况数进行循环
        
    
        
        
    !**确定入射波相关参数**
    if(trim(adjustl(wavectrl))=='Calm' ) then      !静水
        amp=0.0;              !---波幅为0
        wllpp=1.0;            !---波长船长比
        wavel=Lpp;            !---认为波长是1.0倍船长
        wavek=DPi/wavel;      !---根据波长计算波浪频率
        ome=dsqrt(g0*wavek)   !---计算波浪频率
        head=Pi               !---认为静水是迎浪
        wavet=DPi/ome         !---波浪自然周期
        omee=ome-wavek*U0*dcos(head)    !---波浪遭遇频率omee
        wavete=DPi/omee                 !---波浪遭遇周期
        
        wavelscale=1.0;            !----特征波长船长比
        Nramp=int(wavet/Dtsim )    !----平滑函数作用的步数
        Nsimu=Nramp*int(time_num(load_ID))      !----总的模拟步数
        
 !20231113修改，目的是后续添加不规则模块        
    elseif( trim(adjustl(wavectrl))=='Regular' )then        !规则波
        amp=amp_num(load_ID)
        ome=ome_num(load_ID)
        wavel=DPi*g0/ome**2.0        !----波长
        wllpp=wavel/Lpp              !----波长船长比
        wavek=ome**2.0/g0            !----波数
        head=beta_num(load_ID)/180.0*Pi    !---浪向角
        wavet=DPi/ome                      !---波浪自然周期
        omee=ome-wavek*U0*dcos(head)       !---波浪遭遇频率omee
        wavete=DPi/omee                    !---波浪遭遇周期
        wavelscale=DPi*g0/omee**2.0        !----考虑航速修正的特征波长
        
        wavelscale=wavelscale/Lpp          !----特征波长船长比        
        Nramp=int(wavet/Dtsim )            !----平滑函数作用的步数
        Nsimu=int(time_num(load_ID)*wavet/Dtsim )      !---总的模拟步数
!20231113修改        
    elseif( trim(adjustl(wavectrl))=='Irregular' )then
        IrreNum=50;
        head=beta_num(load_ID)/180.*Pi    !---浪向角
        allocate( IrreAmp(1:IrreNum),IrreOme(1:IrreNum),IrreOmee(1:IrreNum) )
        allocate( Irrek(1:IrreNum),Irrepha(1:IrreNum) )
        
        call PreIrreWave(hs_num(load_ID),tz_num(load_ID))  !---自动生成不规则波
        
        do i=1,IrreNum
            Irrek(i)=Irreome(i)**2/g0
            Irreomee(i)=Irreome(i)-Irrek(i)*U0*dcos(head)
        end do
        
        wavel=g0*(tz_num(load_ID)*2.0)**2./Dpi
        ome=2.0*Pi/(tz_num(load_ID)*2.0)
        wllpp=wavel/Lpp
        
        wavelscale=0.1;
        
        Nramp=int(2.0*tz_num(load_ID)/Dtsim )     !光滑函数作用步数
        Nsimu=int(time_num(load_ID)/Dtsim )       !总的迭代步数
 
    end if
    
    
    !**平滑函数设置**
    allocate( smooth(Nramp) )
    smooth=0.0
    do it=1,Nramp
        Smooth(it)=0.5*(1.0-cos(real((it-1) )*Pi/Real(Nramp)))        
    end do
    Smooth(1)=0.0;
    
    !**CHECK**20230506
    open( 1000,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.check' )
    
    
    !**输出文件设置**
    call qinttostr(loadIDname(load_ID),load_name,4)
    
    call outputfile( load_name,load_ID )
    open(unit=SAIL_HISTORY_UNIT, &
        file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))// &
        '_LC'//trim(adjustl(load_name))//'_sail_history.csv', &
        status='replace', action='write')
    write(SAIL_HISTORY_UNIT,'(A)') &
        'time,load_ID,sail_enable,surge,sway,heave,roll,pitch,yaw,'// &
        'surge_rate,sway_rate,heave_rate,roll_rate,pitch_rate,yaw_rate,'// &
        'V_CG_BODY_X,V_CG_BODY_Y,V_CG_BODY_Z,'// &
        'OMEGA_BODY_X,OMEGA_BODY_Y,OMEGA_BODY_Z,'// &
        'V_WIND_BODY_X,V_WIND_BODY_Y,V_WIND_BODY_Z,'// &
        'V_REL_BODY_X,V_REL_BODY_Y,V_REL_BODY_Z,V_REL_H_MAG,'// &
        'E_UPSTREAM_X,E_UPSTREAM_Y,E_UPSTREAM_Z,'// &
        'C_CHORD_X,C_CHORD_Y,C_CHORD_Z,'// &
        'ALPHA_RAW_DEG,ALPHA_DB_DEG,CL,CD,Q_DYNAMIC,'// &
        'FX_BODY,FY_BODY,FZ_BODY,MX_BODY,MY_BODY,MZ_BODY,'// &
        'FX_PLANNING,FY_PLANNING,FZ_PLANNING,'// &
        'MX_PLANNING,MY_PLANNING,MZ_PLANNING'
    
    !**初始化**
    y=0.0;  dy=0.0;            !模态主坐标和一阶导，主坐标对时间的一阶导和二阶导
    t=0.0;                     !时域计算的真实时刻
    derphiI=0.0;  DtphiI=0.0;   derDtPhiI=0.0;      !入射势相关变量
    t=0.0; tt=0.0
    
    a(1)=Dtsim/2.0; a(2)=a(1); a(3)=Dtsim; a(4)=Dtsim      !荣格库塔法中时间变量
    
    zz=0.0;    bb=0.0;
    
   
    MEsecLoad=0.0;    MEsecLoad0=0.0;    MEsecLoad1=0.0

    air_force=0.0;
    
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t0ScogM=0.0;  t1ScogM=0.0;      !区间前后两个时刻，船体摇荡运动的位移
        Stime0=0.0;   Stime1=0.0;        !砰击计算区间首尾时间
        SlamForce=0.0;                  !一个耐波性时间步长里面的平均砰击载荷
        SlamNumT=5;                     !一个完整的耐波性时间步里需要计算的砰击步数
        
        t0SlamForce=0.0;                !区间前端时刻对应的砰击整体载荷
        
        SlamRelaP=0;                  !---0典型节点在波面之上，1典型节点在波面之下
        
        do i=1,SlamNumLine
            if( SlamIniType(i,1,3)<=0.0 ) SlamRelaP(i)=1
            
            SlamCase(i)=0;
            SlamIniPenetration(i)=0.0;
            SlamIniRise(i)=0.0;
   
        end do
    
        
    end if
    
    open(unit=111,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//'.che1')
    open(unit=112,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//'.che2')
    open(unit=113,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//'.che3')
    

    !**时域模拟开始**
    do it=1,Nsimu       !每个工况对总的模拟步数进行循环
        
    !**进行规则波计算时需要载荷调平**先不管20230418
        
        
    !**重新定义光滑函数**
    smtf=0.0;
    if( it<Nramp )then
        smtf=Smooth(it)
    else
        smtf=1.0
    end if
    
    !**这是用来更新运动吗**
    zz(:)=y(:)
    
    !**砰击区间前端时刻以及前端时刻对应的摇荡位移**
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t0ScogM=y(1:6)
        Stime0=t        
    end if
    
    !**计算各剖面梁振动运动MEsecLoad20230427**    
    call SectMotion( y,t )
    
    MEsecLoad0(:,:)=MEsecLoad(:,:)
    
    !**入射势给定**
    do i=1,Nwh
        do j=1,NL
            tmpx(1:3)=Node(i,j,1:3)
            
            !20231113修改
            if( IrreCtrl==0 )then
                
            tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)-omee*t
            derphiI(i,j,1)=amp*ome*exp(wavek*tmpx(3))*cos(head)* cos( tmpt ) 
            derphiI(i,j,2)=amp*ome*exp(wavek*tmpx(3))*sin(head)* cos( tmpt )
            derphiI(i,j,3)=amp*ome*exp(wavek*tmpx(3))*sin( tmpt )
            
            DtphiI(i,j)=-amp*g0*exp(wavek*tmpx(3))* cos( tmpt )
            
            derDtPhiI(i,j,1)=amp*ome**2.0*cos(head)*exp(wavek*tmpx(3))*sin( tmpt )
            derDtPhiI(i,j,2)=amp*ome**2.0*sin(head)*exp(wavek*tmpx(3))*sin( tmpt )
            derDtPhiI(i,j,3)=-amp*ome**2.0*exp(wavek*tmpx(3))*cos( tmpt )
            
            else
                
            derphiI(i,j,:)=0.0; DtphiI(i,j)=0.0; derDtPhiI(i,j,:)=0.0;
            do jj=1,IrreNum
                tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*t +Irrepha(jj)
                derphiI(i,j,1)=derphiI(i,j,1)+Irreamp(jj)*Irreome(jj)*exp(Irrek(jj)*tmpx(3))*cos(head)* cos( tmpt )
                derphiI(i,j,2)=derphiI(i,j,2)+Irreamp(jj)*Irreome(jj)*exp(Irrek(jj)*tmpx(3))*sin(head)* cos( tmpt )
                derphiI(i,j,3)=derphiI(i,j,3)+Irreamp(jj)*Irreome(jj)*exp(Irrek(jj)*tmpx(3))*sin( tmpt )  
                DtphiI(i,j)=DtphiI(i,j)-Irreamp(jj)*g0*exp(Irrek(jj)*tmpx(3))* cos( tmpt ) 
        
                derDtPhiI(i,j,1)=derDtPhiI(i,j,1)+Irreamp(jj)*Irreome(jj)**2.0*cos(head)*exp(Irrek(jj)*tmpx(3))*sin( tmpt ) 
                derDtPhiI(i,j,2)=derDtPhiI(i,j,2)+Irreamp(jj)*Irreome(jj)**2.0*sin(head)*exp(Irrek(jj)*tmpx(3))*sin( tmpt ) 
                derDtPhiI(i,j,3)=derDtPhiI(i,j,3)-Irreamp(jj)*Irreome(jj)**2.0*exp(Irrek(jj)*tmpx(3))*cos( tmpt ) 
            end do
            
            end if
                        
            
        end do        
    end do
    !光滑处理
    if( it<Nramp )then
        derphiI=derphiI*Smooth(it)
        DtphiI=DtphiI*Smooth(it)
        derDtPhiI=derDtPhiI*Smooth(it)    
    end if
    
    
    !更新运动（龙格库塔第一步）
        
    !**计算由刚度矩阵、阻尼矩阵引起的载荷**(负号是因为从方程的左端移到右端)
    sail_print = (load_ID == 1 .and. it == 1)
    call ComputeSailPlanningDryRun(U0, y, NR, t, &
        SAIL_TEST_WIND_GLOBAL, SAIL_TEST_R_BODY, SAIL_TEST_DELTA_DEG, &
        sail_v_cg_body, sail_omega_body, sail_v_wind_body, &
        sail_vrel_h_mag, sail_alpha_db_deg, sail_cl, sail_cd, &
        sail_force_body, sail_moment_body, sail_load_body, &
        sail_ierr, sail_message, 1, sail_print)
    if (.not. all(ieee_is_finite(sail_load_body))) then
        error stop 'Sail dry-run returned non-finite LOAD_6DOF during K1.'
    end if
    if (sail_ierr /= SAIL_OK .and. sail_ierr /= SAIL_ERR_LOW_WIND_SPEED) then
        write(*,'(A,I0,2A)') 'Sail dry-run K1 failed, IERR=', sail_ierr, &
            ': ', trim(sail_message)
        error stop 1
    end if

    call TransformSailLoadToPlanning(y(4:6), sail_load_body, &
        sail_load_planning, sail_transform_ierr, sail_transform_message)
    if (sail_transform_ierr /= SAIL_OK) then
        write(*,'(A,I0,2A)') 'Sail load transform K1 failed, IERR=', &
            sail_transform_ierr, ': ', trim(sail_transform_message)
        error stop 1
    end if

    ForMot(1:NR)=-matmul( CC_mat(1:NR,1:NR),y(1:NR) )-matmul( BB_mat(1:NR,1:NR),y(NR+1:2*NR) )

        
    !**非线性的入射力以及静水压力这里先不考虑20230419**
  
    if( trim( adjustl(Non_Linear))=='LT' ) then         !考虑入射势以及静水压力的非线性       
       
        call Elastic_SurfaceIntegral(DtphiI,ForceI)
        
        !**求方程右端合力**
        ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+ForceI(1:NR)+air_force(1:NR)
        
    elseif( trim( adjustl(Non_Linear))=='NL' )then
        
        if( nonlinearCtrl==0 )then
        
            Inst_ForceIS=0.0;
        
            call Instant_wetsurface( y,it,t,smtf,Inst_ForceIS )
        
            !**求方程右端合力**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Inst_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
        
        elseif( nonlinearCtrl==1 )then
            
            Instbdf_ForceIS=0.0;
            
            call Instant_bdf_wetsurface( y,it,t,smtf,Instbdf_ForceIS )
            
            !**求方程右端合力**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Instbdf_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
            
            
            write(113,'(5f16.4)')t,Instbdf_ForceIS(3),Instbdf_ForceIS(5),air_force(3),air_force(5)
            
        end if
        
        
    end if
    
    
    sail_force_generalized = 0.0_DP
    sail_force_generalized(1:6) = sail_load_planning
    sail_formot_before = ForMot
    sail_acc_before = matmul(invAA_mat, sail_formot_before)
    if (SAIL_ENABLE) ForMot(1:6) = ForMot(1:6) + sail_load_planning
    sail_formot_delta = ForMot - sail_formot_before
    sail_acc_after = matmul(invAA_mat, ForMot)
    sail_acc_expected = matmul(invAA_mat, sail_force_generalized)
    sail_check_scale = max(1.0_DP, maxval(abs(sail_load_planning)))
    if (SAIL_ENABLE) then
        if (maxval(abs(sail_formot_delta(1:6)-sail_load_planning)) > &
            SAIL_COUPLING_TOL*sail_check_scale) &
            error stop 'Sail K1 ForMot increment check failed.'
        if (NR > 6) then
            if (any(sail_formot_delta(7:NR) /= 0.0_DP)) &
                error stop 'Sail K1 modified elastic generalized loads.'
        end if
        sail_check_scale = max(1.0_DP, maxval(abs(sail_acc_expected)))
        if (maxval(abs((sail_acc_after-sail_acc_before)-sail_acc_expected)) > &
            SAIL_COUPLING_TOL*sail_check_scale) &
            error stop 'Sail K1 acceleration increment check failed.'
    else
        if (any(sail_formot_delta /= 0.0_DP)) &
            error stop 'Sail OFF unexpectedly modified ForMot during K1.'
    end if
    if (sail_print) then
        write(*,'(A,3(ES24.16,1X))') 'ANGLES = ', y(4:6)
        write(*,'(A,6(ES24.16,1X))') 'LOAD_6DOF_BODY = ', sail_load_body
        write(*,'(A,6(ES24.16,1X))') 'LOAD_6DOF_PLANNING = ', sail_load_planning
        write(*,'(A,6(ES24.16,1X))') 'ForMot_before_sail = ', sail_formot_before(1:6)
        write(*,'(A,6(ES24.16,1X))') 'ForMot_after_sail = ', ForMot(1:6)
        write(*,'(A,6(ES24.16,1X))') 'ForMot_delta = ', sail_formot_delta(1:6)
        write(*,'(A,6(ES24.16,1X))') 'expected_delta = ', sail_load_planning
        write(*,'(A,6(ES24.16,1X))') 'acc_before = ', sail_acc_before(1:6)
        write(*,'(A,6(ES24.16,1X))') 'acc_after = ', sail_acc_after(1:6)
        write(*,'(A,6(ES24.16,1X))') 'acc_expected_delta = ', sail_acc_expected(1:6)
        write(*,'(A,L1)') 'SAIL_ENABLE = ', SAIL_ENABLE
    end if

    write(111,'(7f16.4)')t,ForMot(1),ForMot(2),ForMot(3),ForMot(4),ForMot(5),ForMot(6)    
    
        
    !!**求方程右端合力**
    !ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+ForceI(1:NR)
        

                
    !**计算载荷所需要的主坐标**
    Secty(:)=y(:)
    Sectdy(1:NR)=y(NR+1:2*NR)
    Sectdy(NR+1:2*NR)=matmul( invAA_mat(1:NR,1:NR),ForMot(1:NR) )

        
    !**更新计算运动所需要的主坐标**（与载荷分开主要是与LT保持一致）
    dy(1:NR)=y(NR+1:2*NR)
    dy(NR+1:2*NR)=matmul( invAA_mat(1:NR,1:NR),ForMot(1:NR) )
        
  
    !**将纵荡、横荡、横摇、首摇强制赋值为零**
        
    if( trim(adjustl(SurgeCtr))=='NO' ) dy(NR+1)=0.0
    if( trim(adjustl(SwayCtr))=='NO' )  dy(NR+2)=0.0
    if( trim(adjustl(RollCtr))=='NO' )  dy(NR+4)=0.0
    if( trim(adjustl(YawCtr))=='NO' )   dy(NR+6)=0.0
    
        
    !**计算剖面载荷，模态叠加法（先不做）**
        
        
    call SectionLoad( Secty,t,smtf,it )
    
    

        
        
        

    
    
    bb(:)=y(:)       !当前时间步运动
    
    !龙格库塔第2，3，4步
    do k=1,3         !对龙格库塔第2，3，4步循环
        
    y(:)=zz(:)+a(k)*dy(:)
    bb(:)=bb(:)+a(k+1)*dy(:)/3.0
    
    !**更新时间**
    tt=t+a(k)

    sail_print = (load_ID == 1 .and. it == 1)
    call ComputeSailPlanningDryRun(U0, y, NR, tt, &
        SAIL_TEST_WIND_GLOBAL, SAIL_TEST_R_BODY, SAIL_TEST_DELTA_DEG, &
        sail_v_cg_body, sail_omega_body, sail_v_wind_body, &
        sail_vrel_h_mag, sail_alpha_db_deg, sail_cl, sail_cd, &
        sail_force_body, sail_moment_body, sail_load_body, &
        sail_ierr, sail_message, k+1, sail_print)
    if (.not. all(ieee_is_finite(sail_load_body))) then
        error stop 'Sail dry-run returned non-finite LOAD_6DOF during K2-K4.'
    end if
    if (sail_ierr /= SAIL_OK .and. sail_ierr /= SAIL_ERR_LOW_WIND_SPEED) then
        write(*,'(A,I0,A,I0,2A)') 'Sail dry-run K', k+1, &
            ' failed, IERR=', sail_ierr, ': ', trim(sail_message)
        error stop 1
    end if

    call TransformSailLoadToPlanning(y(4:6), sail_load_body, &
        sail_load_planning, sail_transform_ierr, sail_transform_message)
    if (sail_transform_ierr /= SAIL_OK) then
        write(*,'(A,I0,A,I0,2A)') 'Sail load transform K', k+1, &
            ' failed, IERR=', sail_transform_ierr, ': ', &
            trim(sail_transform_message)
        error stop 1
    end if

    
    !**龙格库塔第234步入射势给定**
    do i=1,Nwh
        do j=1,NL
            tmpx(1:3)=Node(i,j,1:3)
            
            !20231113修改
            if( IrreCtrl==0 )then
                
            tmpt=wavek*tmpx(1)*cos(head)+wavek*tmpx(2)*sin(head)- omee*tt
            derphiI(i,j,1)=amp*ome*exp(wavek*tmpx(3))*cos(head)* cos( tmpt )
            derphiI(i,j,2)=amp*ome*exp(wavek*tmpx(3))*sin(head)* cos( tmpt )
            derphiI(i,j,3)=amp*ome*exp(wavek*tmpx(3))*sin( tmpt )        
            DtphiI(i,j)=-amp*g0*exp(wavek*tmpx(3))* cos( tmpt )
            
            derDtPhiI(i,j,1)=amp*ome**2.0*cos(head)*exp(wavek*tmpx(3))*sin( tmpt )
            derDtPhiI(i,j,2)=amp*ome**2.0*sin(head)*exp(wavek*tmpx(3))*sin( tmpt )
            derDtPhiI(i,j,3)=-amp*ome**2.0*exp(wavek*tmpx(3))*cos( tmpt )
            
            else
                
            derphiI(i,j,:)=0.0; DtphiI(i,j)=0.0; derDtPhiI(i,j,:)=0.0;
            do jj=1,IrreNum
                tmpt=Irrek(jj)*tmpx(1)*cos(head)+Irrek(jj)*tmpx(2)*sin(head)- Irreomee(jj)*tt +Irrepha(jj)
                derphiI(i,j,1)=derphiI(i,j,1)+Irreamp(jj)*Irreome(jj)*exp(Irrek(jj)*tmpx(3))*cos(head)* cos( tmpt )
                derphiI(i,j,2)=derphiI(i,j,2)+Irreamp(jj)*Irreome(jj)*exp(Irrek(jj)*tmpx(3))*sin(head)* cos( tmpt )
                derphiI(i,j,3)=derphiI(i,j,3)+Irreamp(jj)*Irreome(jj)*exp(Irrek(jj)*tmpx(3))*sin( tmpt )  
                DtphiI(i,j)=DtphiI(i,j)-Irreamp(jj)*g0*exp(Irrek(jj)*tmpx(3))* cos( tmpt ) 
          
                derDtPhiI(i,j,1)=derDtPhiI(i,j,1)+Irreamp(jj)*Irreome(jj)**2.0*cos(head)*exp(Irrek(jj)*tmpx(3))*sin( tmpt ) 
                derDtPhiI(i,j,2)=derDtPhiI(i,j,2)+Irreamp(jj)*Irreome(jj)**2.0*sin(head)*exp(Irrek(jj)*tmpx(3))*sin( tmpt ) 
                derDtPhiI(i,j,3)=derDtPhiI(i,j,3)-Irreamp(jj)*Irreome(jj)**2.0*exp(Irrek(jj)*tmpx(3))*cos( tmpt )            
            end do
                
            end if
            
            
        end do
    end do
    
    !光滑处理
    if( it<Nramp )then
        derphiI=derphiI*Smooth(it)
        DtphiI=DtphiI*Smooth(it)
        derDtPhiI=derDtPhiI*Smooth(it)    
    end if
    
    !**计算各剖面梁振动运动MEsecLoad20230427**    
    call SectMotion( y,tt )
    MEsecLoad1(:,:)=MEsecLoad(:,:)
    
    
    !计算龙格库塔第2，3，4步气动升力
    if( trim(adjustl(Airlift))=='YES' )then
        call airforce( y(5) )
        
    end if
    
    
    
    
    !计算龙格库塔第2，3，4步的水动升力
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t1ScogM(1:6)=y(1:6)    !---t1时刻对应的船体摇荡运动位移
        Stime1=tt;             !---t1时刻
        
        if( k==1.or.k==2 )then
            i=SlamNumT/2            
        elseif( k==3 )then
            i=SlamNumT            
        end if
        
        call SlamCase4(k,smtf,i,Stime0,Stime1,t0ScogM,t1ScogM,MEsecLoad0,MEsecLoad1,t0SlamForce,SlamForce  )
        
        !SlamForce(3)=0.0;
    
    
    end if
    
    !write(111,'(2f16.6)')  SlamForce(3),SlamForce(5)
    

    
    
    !更新运动（龙格库塔第234步）
 
    !**计算由刚度矩阵、阻尼矩阵引起的载荷**(负号是因为从方程的左端移到右端)
    ForMot(1:NR)=-matmul( CC_mat(1:NR,1:NR),y(1:NR) )-matmul( BB_mat(1:NR,1:NR),y(NR+1:2*NR) )
        
    !**非线性的入射力以及静水压力这里先不考虑20230419**
    
    
    
    
    if( trim( adjustl(Non_Linear))=='LT' ) then         !考虑入射势以及静水压力的非线性       
       
        call Elastic_SurfaceIntegral(DtphiI,ForceI)
        
        !**求方程右端合力**
        ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+ForceI(1:NR)+air_force(1:NR)
        
    elseif( trim( adjustl(Non_Linear))=='NL' )then
        
        if( nonlinearCtrl==0 )then
        
            Inst_ForceIS=0.0;
        
            call Instant_wetsurface( y,it,t,smtf,Inst_ForceIS )
        
            !**求方程右端合力**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Inst_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
        
        elseif( nonlinearCtrl==1 )then
            
            Instbdf_ForceIS=0.0;
            
            call Instant_bdf_wetsurface( y,it,t,smtf,Instbdf_ForceIS )
            
            !**求方程右端合力**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Instbdf_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
            
        
        end if
        
        
    end if
        

        
        
    !**更新计算运动所需要的主坐标**（与载荷分开主要是与LT保持一致）
    sail_force_generalized = 0.0_DP
    sail_force_generalized(1:6) = sail_load_planning
    sail_formot_before = ForMot
    if (SAIL_ENABLE) ForMot(1:6) = ForMot(1:6) + sail_load_planning
    sail_formot_delta = ForMot - sail_formot_before
    sail_check_scale = max(1.0_DP, maxval(abs(sail_load_planning)))
    if (SAIL_ENABLE) then
        if (maxval(abs(sail_formot_delta(1:6)-sail_load_planning)) > &
            SAIL_COUPLING_TOL*sail_check_scale) &
            error stop 'Sail K2-K4 ForMot increment check failed.'
        if (NR > 6) then
            if (any(sail_formot_delta(7:NR) /= 0.0_DP)) &
                error stop 'Sail K2-K4 modified elastic generalized loads.'
        end if
    else
        if (any(sail_formot_delta /= 0.0_DP)) &
            error stop 'Sail OFF unexpectedly modified ForMot during K2-K4.'
    end if
    if (sail_print) then
        write(*,'(A,3(ES24.16,1X))') 'ANGLES = ', y(4:6)
        write(*,'(A,6(ES24.16,1X))') 'LOAD_6DOF_BODY = ', sail_load_body
        write(*,'(A,6(ES24.16,1X))') 'LOAD_6DOF_PLANNING = ', sail_load_planning
        write(*,'(A,6(ES24.16,1X))') 'ForMot_before_sail = ', sail_formot_before(1:6)
        write(*,'(A,6(ES24.16,1X))') 'ForMot_after_sail = ', ForMot(1:6)
        write(*,'(A,6(ES24.16,1X))') 'ForMot_delta = ', sail_formot_delta(1:6)
        write(*,'(A,6(ES24.16,1X))') 'expected_delta = ', sail_load_planning
        write(*,'(A,L1)') 'SAIL_ENABLE = ', SAIL_ENABLE
    end if
    dy(1:NR)=y(NR+1:2*NR)
    dy(NR+1:2*NR)=matmul( invAA_mat(1:NR,1:NR),ForMot(1:NR) )
        
    !**将纵荡、横荡、横摇、首摇强制赋值为零**        
    if( trim(adjustl(SurgeCtr))=='NO' ) dy(NR+1)=0.0
    if( trim(adjustl(SwayCtr))=='NO' )  dy(NR+2)=0.0
    if( trim(adjustl(RollCtr))=='NO' )  dy(NR+4)=0.0
    if( trim(adjustl(YawCtr))=='NO' )   dy(NR+6)=0.0
    
    

    
 
    end do          !对龙格库塔第2，3，4步循环
    
    !**一个完整的龙格库塔时间步**
    y(:)=bb(:)+a(1)*dy(:)/3.0
    
    
 
    !**更新时间**
    t=t+a(4)
    
    
    !**计算各剖面梁振动运动MEsecLoad20230427**    
    call SectMotion( y,t )
    
    MEsecLoad1(:,:)=MEsecLoad(:,:)
    
    
    !计算龙格库塔第一步气动升力
    if( trim(adjustl(Airlift))=='YES' )then
        call airforce( y(5) )
        
    end if
    
    write( 1000,'(2f16.4)' ) air_force(3),air_force(5)
    
    !计算龙格库塔第1步的水动升力
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t1ScogM(1:6)=y(1:6)         !---t1时刻对应的船体摇荡运动位移
        Stime1=t;             !---t1时刻
        
        call SlamCase3( smtf,SlamNumT,Stime0,Stime1,t0ScogM,t1ScogM,MEsecLoad0,MEsecLoad1,t0SlamForce,SlamForce )
        
        !SlamForce(3)=0.0;        !参考lt但要弄清楚是为啥
        
    end if
    
    write(112,'(3f16.4)') t,SlamForce(3),SlamForce(5)


    ! Recompute Sail diagnostics at the accepted full-step state. This call is
    ! output-only: the returned load is never added to ForMot here.
    call ComputeSailPlanningDryRun(U0, y, NR, t, &
        SAIL_TEST_WIND_GLOBAL, SAIL_TEST_R_BODY, SAIL_TEST_DELTA_DEG, &
        sail_v_cg_body, sail_omega_body, sail_v_wind_body, &
        sail_vrel_h_mag, sail_alpha_db_deg, sail_cl, sail_cd, &
        sail_force_body, sail_moment_body, sail_load_body, &
        sail_ierr, sail_message, 0, .false., &
        sail_v_rel_body, sail_e_upstream_body, sail_c_chord_body, &
        sail_alpha_raw_deg, sail_q_dynamic)
    if (sail_ierr /= SAIL_OK .and. sail_ierr /= SAIL_ERR_LOW_WIND_SPEED) then
        write(*,'(A,I0,2A)') 'Sail step-end diagnostic failed, IERR=', &
            sail_ierr, ': ', trim(sail_message)
        error stop 1
    end if
    call TransformSailLoadToPlanning(y(4:6), sail_load_body, &
        sail_load_planning, sail_transform_ierr, sail_transform_message)
    if (sail_transform_ierr /= SAIL_OK) then
        write(*,'(A,I0,2A)') 'Sail step-end transform failed, IERR=', &
            sail_transform_ierr, ': ', trim(sail_transform_message)
        error stop 1
    end if
    if (.not. all(ieee_is_finite([y(1:6), y(NR+1:NR+6), &
        sail_v_cg_body, sail_omega_body, sail_v_wind_body, &
        sail_v_rel_body, sail_e_upstream_body, sail_c_chord_body, &
        sail_vrel_h_mag, sail_alpha_raw_deg, sail_alpha_db_deg, &
        sail_cl, sail_cd, sail_q_dynamic, sail_load_body, &
        sail_load_planning]))) error stop 'Non-finite Sail history value.'
    write(SAIL_HISTORY_UNIT,'(*(G0,:,","))') t,load_ID,SAIL_ENABLE, &
        y(1:6),y(NR+1:NR+6),sail_v_cg_body,sail_omega_body, &
        sail_v_wind_body,sail_v_rel_body,sail_vrel_h_mag, &
        sail_e_upstream_body,sail_c_chord_body,sail_alpha_raw_deg, &
        sail_alpha_db_deg,sail_cl,sail_cd,sail_q_dynamic, &
        sail_load_body,sail_load_planning    


    
    
    
    !**输出位移、加速度**
    
    write(23,"(f12.3,1x,\)",advance='NO') t
    write(22,"(f12.3,1x,\)",advance='NO') t
    do i=1,NR
        write(23,"(e12.4,1x,\)",advance='NO') y(i)
        write(22,"(e12.4,1x,\)",advance='NO') (y(NR+i)-zz(NR+i))/Dtsim        
    end do
    write(23,"(/)")
    write(22,"(/)")
    
    
    write(*,"(A,i6,A,i6,A,i6,A,i8)") '当前工况',load_ID,' /',loadnum,'  时域模拟进度',it,'   /',Nsimu
    
    !20240418加绘图
    !call system(' start '//'C:\Users\lic\Desktop\suanli1\tempout\heave.plt')
    
  
    end do     !每个工况对总的模拟步数进行循环
    
    close(111);close(112);close(113)
    close(SAIL_HISTORY_UNIT)
    
    
    close(22);close(23);close(24);close(25);close(26);close(27);close(28);close(29)
    
    deallocate( smooth )
    
    close(1000)
    

    

    if( trim(adjustl(Slamctrl))=='YES' ) then        !平均砰击载荷水平力、垂向力
        close(65); close(66);
    end if

    !20231115修改,添加船中纵剖线各点与波面的相对位置计算
    if( trim(adjustl(SlamRelaMctr))=='YES' )then
        close(67); close(68); close(69); close(70);

    end if
    
    
    
    
    !20231113修改
    if( trim(adjustl(wavectrl))=='Irregular' )then
        deallocate( IrreAmp,IrreOme,IrreOmee,Irrek,Irrepha )        
    end if
    
    
    
    end do      !对工况数进行循环
        
    
    
    
    
    
    
    
    
end subroutine Timmotion
