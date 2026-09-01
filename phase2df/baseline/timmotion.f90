subroutine Timmotion
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use ShipHullVar
    use Constant
    use PanelGeometry
    use Slamming
    use ArrayOperations
    !20231113�޸�
    use IrreWaveVar
    
    use lift
    use SAILPARAM_MOD, only: DP, SAIL_OK, SAIL_ERR_LOW_WIND_SPEED
    use SAILPLANNING_ADAPTER_MOD, only: ComputeSailPlanningDryRun, &
        TransformSailLoadToPlanning
    
    implicit none
    
    integer(4)::i,j,k,it,ii,jj
    

    
    real(8)::wavelscale       !----��������������
    
    integer(4)::load_ID
    character(len=4)::load_name
    
    real(8),allocatable,dimension(:)::smooth
    real(8)::smtf
    
    real(8),allocatable,dimension(:,:)::AA_mat,BB_mat,CC_mat,invAA_mat     !�����������ᡢ�նȡ��������������
    
    real(8),dimension(2*NR)::bb,zz,y,dy     !ģ̬�������һ�׵����������ʱ���һ�׵��Ͷ��׵�(�����˶�)
    real(8),dimension(2*NR)::Secty,Sectdy      !��y,dy��Ӧ�����ڸ����غɣ�
    real*8,dimension(4)::a                    !�������ʱ�䲽
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
    
    !��mlm���ۼ���ˮ��������صı���
    real(8),dimension(6)::t0ScogM,t1ScogM        !����ǰ������ʱ�̣�����ҡ���˶���λ��
    real(8)::Stime0,Stime1                       !�������������βʱ��
    real(8),dimension(NR)::SlamForce             !һ���Ͳ���ʱ�䲽�����ƽ������غ�
    integer(4)::SlamNumT                         !һ���������Ͳ���ʱ�䲽����Ҫ������������
    real(8),dimension(NR)::t0SlamForce           !����ǰ��ʱ�̶�Ӧ����������غ�
    
    
    
    
    !��������ر���
    real(8),dimension(3)::tmpx
    real(8)::tmpt
    real*8,dimension(Nwh,NL,3)::derPhiI
    real*8,dimension(Nwh,NL)::DtPhiI
    real(8),dimension(Nwh,NL,3)::derDtPhiI
    
    !��˲ʱ����������ˮѹ����صı���
    real(8),dimension(NR)::Inst_ForceIS      !��һ�ַ�ʽ������������
    
    
    real(8),dimension(NR)::Instbdf_ForceIS     !�ڶ��ַ�ʽ������������
    
    
    !�������˶��йصı���
    real(8),dimension( 1:6,1:SlamNumLine )::MEsecLoad0,MEsecLoad1
    
    

    
    
    !****ǰ��Ϊ��������****
    
    !**���巽���Ҷ˺�����**
    allocate( ForMot(NR),ForceI(NR) )
    ForMot=0.0;    ForceI=0.0
    
    !**�����������浯����**
    allocate( MEsecLoad(6,SlamNumLine) )
    
    
    !**��������**
    allocate( air_force(NR) )
    air_force=0.0;
    
    
    !**����������
    
    
    
    
    !**����ϵ������**
    allocate( AA_mat(NR,NR),BB_mat(NR,NR),CC_mat(NR,NR),invAA_mat(NR,NR) )
    AA_mat=0.0;  BB_mat=0.0;   CC_mat=0.0;   invAA_mat=0.0
    
    AA_mat(1:NR,1:NR)=StruM(1:NR,1:NR)
    BB_mat(1:NR,1:NR)=StruB(1:NR,1:NR)
    
    !�ٽ���������ҡ����ϵ������20230519
    BB_mat(3,3)=BcoefH*2.0*sqrt( StruM(3,3)*HERM(3,3) )
    BB_mat(5,5)=BcoefP*2.0*sqrt( StruM(5,5)*HERM(5,5) )
    
    
    !CC_mat(1:NR,1:NR)=StruC(1:NR,1:NR)+HERM(1:NR,1:NR)
    

    !20230519��ȥ��ˮ�ָ�������    
    if( trim( adjustl(Non_Linear))=='LT' )then
        CC_mat(1:NR,1:NR)=StruC(1:NR,1:NR)+HERM(1:NR,1:NR)        
    else        
        CC_mat(1:NR,1:NR)=StruC(1:NR,1:NR)        
    end if
    
    
    !**�����������������**
    call Gauss_Jordan( NR,NR,AA_mat,invAA_mat )
    
    
   
        
        
    do load_ID=1,loadnum      !�Թ���������ѭ��
        
    
        
        
    !**ȷ�����䲨��ز���**
    if(trim(adjustl(wavectrl))=='Calm' ) then      !��ˮ
        amp=0.0;              !---����Ϊ0
        wllpp=1.0;            !---����������
        wavel=Lpp;            !---��Ϊ������1.0������
        wavek=DPi/wavel;      !---���ݲ������㲨��Ƶ��
        ome=dsqrt(g0*wavek)   !---���㲨��Ƶ��
        head=Pi               !---��Ϊ��ˮ��ӭ��
        wavet=DPi/ome         !---������Ȼ����
        omee=ome-wavek*U0*dcos(head)    !---��������Ƶ��omee
        wavete=DPi/omee                 !---������������
        
        wavelscale=1.0;            !----��������������
        Nramp=int(wavet/Dtsim )    !----ƽ���������õĲ���
        Nsimu=Nramp*int(time_num(load_ID))      !----�ܵ�ģ�ⲽ��
        
 !20231113�޸ģ�Ŀ���Ǻ�����Ӳ�����ģ��        
    elseif( trim(adjustl(wavectrl))=='Regular' )then        !����
        amp=amp_num(load_ID)
        ome=ome_num(load_ID)
        wavel=DPi*g0/ome**2.0        !----����
        wllpp=wavel/Lpp              !----����������
        wavek=ome**2.0/g0            !----����
        head=beta_num(load_ID)/180.0*Pi    !---�����
        wavet=DPi/ome                      !---������Ȼ����
        omee=ome-wavek*U0*dcos(head)       !---��������Ƶ��omee
        wavete=DPi/omee                    !---������������
        wavelscale=DPi*g0/omee**2.0        !----���Ǻ�����������������
        
        wavelscale=wavelscale/Lpp          !----��������������        
        Nramp=int(wavet/Dtsim )            !----ƽ���������õĲ���
        Nsimu=int(time_num(load_ID)*wavet/Dtsim )      !---�ܵ�ģ�ⲽ��
!20231113�޸�        
    elseif( trim(adjustl(wavectrl))=='Irregular' )then
        IrreNum=50;
        head=beta_num(load_ID)/180.*Pi    !---�����
        allocate( IrreAmp(1:IrreNum),IrreOme(1:IrreNum),IrreOmee(1:IrreNum) )
        allocate( Irrek(1:IrreNum),Irrepha(1:IrreNum) )
        
        call PreIrreWave(hs_num(load_ID),tz_num(load_ID))  !---�Զ����ɲ�����
        
        do i=1,IrreNum
            Irrek(i)=Irreome(i)**2/g0
            Irreomee(i)=Irreome(i)-Irrek(i)*U0*dcos(head)
        end do
        
        wavel=g0*(tz_num(load_ID)*2.0)**2./Dpi
        ome=2.0*Pi/(tz_num(load_ID)*2.0)
        wllpp=wavel/Lpp
        
        wavelscale=0.1;
        
        Nramp=int(2.0*tz_num(load_ID)/Dtsim )     !�⻬�������ò���
        Nsimu=int(time_num(load_ID)/Dtsim )       !�ܵĵ�������
 
    end if
    
    
    !**ƽ����������**
    allocate( smooth(Nramp) )
    smooth=0.0
    do it=1,Nramp
        Smooth(it)=0.5*(1.0-cos(real((it-1) )*Pi/Real(Nramp)))        
    end do
    Smooth(1)=0.0;
    
    !**CHECK**20230506
    open( 1000,file=trim(adjustl(InAccess))//'\'//trim(adjustl(projname))//'.check' )
    
    
    !**����ļ�����**
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
    
    !**��ʼ��**
    y=0.0;  dy=0.0;            !ģ̬�������һ�׵����������ʱ���һ�׵��Ͷ��׵�
    t=0.0;                     !ʱ��������ʵʱ��
    derphiI=0.0;  DtphiI=0.0;   derDtPhiI=0.0;      !��������ر���
    t=0.0; tt=0.0
    
    a(1)=Dtsim/2.0; a(2)=a(1); a(3)=Dtsim; a(4)=Dtsim      !�ٸ��������ʱ�����
    
    zz=0.0;    bb=0.0;
    
   
    MEsecLoad=0.0;    MEsecLoad0=0.0;    MEsecLoad1=0.0

    air_force=0.0;
    
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t0ScogM=0.0;  t1ScogM=0.0;      !����ǰ������ʱ�̣�����ҡ���˶���λ��
        Stime0=0.0;   Stime1=0.0;        !�������������βʱ��
        SlamForce=0.0;                  !һ���Ͳ���ʱ�䲽�������ƽ������غ�
        SlamNumT=5;                     !һ���������Ͳ���ʱ�䲽����Ҫ������������
        
        t0SlamForce=0.0;                !����ǰ��ʱ�̶�Ӧ����������غ�
        
        SlamRelaP=0;                  !---0���ͽڵ��ڲ���֮�ϣ�1���ͽڵ��ڲ���֮��
        
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
    

    !**ʱ��ģ�⿪ʼ**
    do it=1,Nsimu       !ÿ���������ܵ�ģ�ⲽ������ѭ��
        
    !**���й��򲨼���ʱ��Ҫ�غɵ�ƽ**�Ȳ���20230418
        
        
    !**���¶���⻬����**
    smtf=0.0;
    if( it<Nramp )then
        smtf=Smooth(it)
    else
        smtf=1.0
    end if
    
    !**�������������˶���**
    zz(:)=y(:)
    
    !**�������ǰ��ʱ���Լ�ǰ��ʱ�̶�Ӧ��ҡ��λ��**
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t0ScogM=y(1:6)
        Stime0=t        
    end if
    
    !**��������������˶�MEsecLoad20230427**    
    call SectMotion( y,t )
    
    MEsecLoad0(:,:)=MEsecLoad(:,:)
    
    !**�����Ƹ���**
    do i=1,Nwh
        do j=1,NL
            tmpx(1:3)=Node(i,j,1:3)
            
            !20231113�޸�
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
    !�⻬����
    if( it<Nramp )then
        derphiI=derphiI*Smooth(it)
        DtphiI=DtphiI*Smooth(it)
        derDtPhiI=derDtPhiI*Smooth(it)    
    end if
    
    
    !�����˶������������һ����
        
    !**�����ɸնȾ����������������غ�**(��������Ϊ�ӷ��̵�����Ƶ��Ҷ�)
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

        
    !**�����Ե��������Լ���ˮѹ�������Ȳ�����20230419**
  
    if( trim( adjustl(Non_Linear))=='LT' ) then         !�����������Լ���ˮѹ���ķ�����       
       
        call Elastic_SurfaceIntegral(DtphiI,ForceI)
        
        !**�󷽳��Ҷ˺���**
        ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+ForceI(1:NR)+air_force(1:NR)
        
    elseif( trim( adjustl(Non_Linear))=='NL' )then
        
        if( nonlinearCtrl==0 )then
        
            Inst_ForceIS=0.0;
        
            call Instant_wetsurface( y,it,t,smtf,Inst_ForceIS )
        
            !**�󷽳��Ҷ˺���**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Inst_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
        
        elseif( nonlinearCtrl==1 )then
            
            Instbdf_ForceIS=0.0;
            
            call Instant_bdf_wetsurface( y,it,t,smtf,Instbdf_ForceIS )
            
            !**�󷽳��Ҷ˺���**
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
    
        
    !!**�󷽳��Ҷ˺���**
    !ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+ForceI(1:NR)
        

                
    !**�����غ�����Ҫ��������**
    Secty(:)=y(:)
    Sectdy(1:NR)=y(NR+1:2*NR)
    Sectdy(NR+1:2*NR)=matmul( invAA_mat(1:NR,1:NR),ForMot(1:NR) )

        
    !**���¼����˶�����Ҫ��������**�����غɷֿ���Ҫ����LT����һ�£�
    dy(1:NR)=y(NR+1:2*NR)
    dy(NR+1:2*NR)=matmul( invAA_mat(1:NR,1:NR),ForMot(1:NR) )
        
  
    !**���ݵ����ᵴ����ҡ����ҡǿ�Ƹ�ֵΪ��**
        
    if( trim(adjustl(SurgeCtr))=='NO' ) dy(NR+1)=0.0
    if( trim(adjustl(SwayCtr))=='NO' )  dy(NR+2)=0.0
    if( trim(adjustl(RollCtr))=='NO' )  dy(NR+4)=0.0
    if( trim(adjustl(YawCtr))=='NO' )   dy(NR+6)=0.0
    
        
    !**���������غɣ�ģ̬���ӷ����Ȳ�����**
        
        
    call SectionLoad( Secty,t,smtf,it )
    
    

        
        
        

    
    
    bb(:)=y(:)       !��ǰʱ�䲽�˶�
    
    !���������2��3��4��
    do k=1,3         !�����������2��3��4��ѭ��
        
    y(:)=zz(:)+a(k)*dy(:)
    bb(:)=bb(:)+a(k+1)*dy(:)/3.0
    
    !**����ʱ��**
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

    
    !**���������234�������Ƹ���**
    do i=1,Nwh
        do j=1,NL
            tmpx(1:3)=Node(i,j,1:3)
            
            !20231113�޸�
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
    
    !�⻬����
    if( it<Nramp )then
        derphiI=derphiI*Smooth(it)
        DtphiI=DtphiI*Smooth(it)
        derDtPhiI=derDtPhiI*Smooth(it)    
    end if
    
    !**��������������˶�MEsecLoad20230427**    
    call SectMotion( y,tt )
    MEsecLoad1(:,:)=MEsecLoad(:,:)
    
    
    !�������������2��3��4����������
    if( trim(adjustl(Airlift))=='YES' )then
        call airforce( y(5) )
        
    end if
    
    
    
    
    !�������������2��3��4����ˮ������
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t1ScogM(1:6)=y(1:6)    !---t1ʱ�̶�Ӧ�Ĵ���ҡ���˶�λ��
        Stime1=tt;             !---t1ʱ��
        
        if( k==1.or.k==2 )then
            i=SlamNumT/2            
        elseif( k==3 )then
            i=SlamNumT            
        end if
        
        call SlamCase4(k,smtf,i,Stime0,Stime1,t0ScogM,t1ScogM,MEsecLoad0,MEsecLoad1,t0SlamForce,SlamForce  )
        
        !SlamForce(3)=0.0;
    
    
    end if
    
    !write(111,'(2f16.6)')  SlamForce(3),SlamForce(5)
    

    
    
    !�����˶������������234����
 
    !**�����ɸնȾ����������������غ�**(��������Ϊ�ӷ��̵�����Ƶ��Ҷ�)
    ForMot(1:NR)=-matmul( CC_mat(1:NR,1:NR),y(1:NR) )-matmul( BB_mat(1:NR,1:NR),y(NR+1:2*NR) )
        
    !**�����Ե��������Լ���ˮѹ�������Ȳ�����20230419**
    
    
    
    
    if( trim( adjustl(Non_Linear))=='LT' ) then         !�����������Լ���ˮѹ���ķ�����       
       
        call Elastic_SurfaceIntegral(DtphiI,ForceI)
        
        !**�󷽳��Ҷ˺���**
        ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+ForceI(1:NR)+air_force(1:NR)
        
    elseif( trim( adjustl(Non_Linear))=='NL' )then
        
        if( nonlinearCtrl==0 )then
        
            Inst_ForceIS=0.0;
        
            call Instant_wetsurface( y,it,t,smtf,Inst_ForceIS )
        
            !**�󷽳��Ҷ˺���**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Inst_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
        
        elseif( nonlinearCtrl==1 )then
            
            Instbdf_ForceIS=0.0;
            
            call Instant_bdf_wetsurface( y,it,t,smtf,Instbdf_ForceIS )
            
            !**�󷽳��Ҷ˺���**
            ForMot(1:NR)=ForMot(1:NR)+SlamForce(1:NR)+Instbdf_ForceIS(1:NR)+air_force(1:NR)+Initial2_Mg(1:NR)
            
        
        end if
        
        
    end if
        

        
        
    !**���¼����˶�����Ҫ��������**�����غɷֿ���Ҫ����LT����һ�£�
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
        
    !**���ݵ����ᵴ����ҡ����ҡǿ�Ƹ�ֵΪ��**        
    if( trim(adjustl(SurgeCtr))=='NO' ) dy(NR+1)=0.0
    if( trim(adjustl(SwayCtr))=='NO' )  dy(NR+2)=0.0
    if( trim(adjustl(RollCtr))=='NO' )  dy(NR+4)=0.0
    if( trim(adjustl(YawCtr))=='NO' )   dy(NR+6)=0.0
    
    

    
 
    end do          !�����������2��3��4��ѭ��
    
    !**һ���������������ʱ�䲽**
    y(:)=bb(:)+a(1)*dy(:)/3.0
    
    
 
    !**����ʱ��**
    t=t+a(4)
    
    
    !**��������������˶�MEsecLoad20230427**    
    call SectMotion( y,t )
    
    MEsecLoad1(:,:)=MEsecLoad(:,:)
    
    
    !�������������һ����������
    if( trim(adjustl(Airlift))=='YES' )then
        call airforce( y(5) )
        
    end if
    
    write( 1000,'(2f16.4)' ) air_force(3),air_force(5)
    
    !�������������1����ˮ������
    if( trim(adjustl(Liftctrl))=='MLM' )then
        t1ScogM(1:6)=y(1:6)         !---t1ʱ�̶�Ӧ�Ĵ���ҡ���˶�λ��
        Stime1=t;             !---t1ʱ��
        
        call SlamCase3( smtf,SlamNumT,Stime0,Stime1,t0ScogM,t1ScogM,MEsecLoad0,MEsecLoad1,t0SlamForce,SlamForce )
        
        !SlamForce(3)=0.0;        !�ο�lt��ҪŪ�����Ϊɶ
        
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


    
    
    
    !**���λ�ơ����ٶ�**
    
    write(23,"(f12.3,1x,\)",advance='NO') t
    write(22,"(f12.3,1x,\)",advance='NO') t
    do i=1,NR
        write(23,"(e12.4,1x,\)",advance='NO') y(i)
        write(22,"(e12.4,1x,\)",advance='NO') (y(NR+i)-zz(NR+i))/Dtsim        
    end do
    write(23,"(/)")
    write(22,"(/)")
    
    
    write(*,"(A,i6,A,i6,A,i6,A,i8)") '��ǰ����',load_ID,' /',loadnum,'  ʱ��ģ�����',it,'   /',Nsimu
    
    !20240418�ӻ�ͼ
    !call system(' start '//'C:\Users\lic\Desktop\suanli1\tempout\heave.plt')
    
  
    end do     !ÿ���������ܵ�ģ�ⲽ������ѭ��
    
    close(111);close(112);close(113)
    close(SAIL_HISTORY_UNIT)
    
    
    close(22);close(23);close(24);close(25);close(26);close(27);close(28);close(29)
    
    deallocate( smooth )
    
    close(1000)
    

    

    if( trim(adjustl(Slamctrl))=='YES' ) then        !ƽ������غ�ˮƽ����������
        close(65); close(66);
    end if

    !20231115�޸�,��Ӵ��������߸����벨������λ�ü���
    if( trim(adjustl(SlamRelaMctr))=='YES' )then
        close(67); close(68); close(69); close(70);

    end if
    
    
    
    
    !20231113�޸�
    if( trim(adjustl(wavectrl))=='Irregular' )then
        deallocate( IrreAmp,IrreOme,IrreOmee,Irrek,Irrepha )        
    end if
    
    
    
    end do      !�Թ���������ѭ��
        
    
    
    
    
    
    
    
    
end subroutine Timmotion
