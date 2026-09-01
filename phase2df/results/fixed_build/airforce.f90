subroutine airforce(pitch)
    use lift
    use Constant,only:pi,TrimAng,NR,U0
    use ShipHullVar,only:Cog0
    
    implicit none
    
    real(8)::pitch
    
    integer(4)::i,j,k
    integer(4)::air_it=500        !积分微元数
    
    real(8),allocatable,dimension(:)::SCL,SXcp,MCL,MXcp       !升力系数、压力中心
    real(8),allocatable,dimension(:)::SA1,SA2,MA1,MA2          !系数A1、A2
    
    real(8),allocatable,dimension(:)::dydx
    
    real(8),allocatable,dimension(:)::Sair_force,Mair_force     !升力计算临时变量
    
    real(8),allocatable,dimension(:)::temp_Sca,temp_Mca        !计算翼元面积的中间弦长
    
    real(8),allocatable,dimension(:)::Sair_area,Mair_area       !翼元面积
    
    real(8),allocatable,dimension(:)::Sair_xx,Mair_xx           !各翼元到重心处的距离
    
    !******以上为变量定义******
    
    !**计算升力系数、压力中心位置**
    allocate( SCL(Swing_num),SXcp(Swing_num),MCL(Mwing_num),MXcp(Mwing_num) )
    SCL=0.0;SXcp=0.0;MCL=0.0;MXcp=0.0;
    
    allocate( SarfaL0(Swing_num),MarfaL0(Mwing_num) )      !零升力攻角
    SarfaL0=0.0;MarfaL0=0.0
    
    allocate( SA1(Swing_num),SA2(Swing_num),MA1(Mwing_num),MA2(Mwing_num) )
    SA1=0.0;SA2=0.0;MA1=0.0;MA2=0.0
    
    !**计算瞬时攻角**
    arfa=pitch+TrimAng
    
    !计算尾翼零升力攻角以及系数A1、A2
    allocate( dydx(air_it) )
    dydx=0.0
    
    do i=1,Swing_num
        !先计算dydx
        call get_dydx( air_it,Sca(i),Spoint_num(i),Swing_node(i,:,:),dydx )
        
        !计算零升力攻角
        do j=1,air_it
            SarfaL0(i)=SarfaL0(i)-dydx(j)*(1-cos(j/real(air_it)*pi))/real(air_it)
        end do
        
        !计算系数SA1,SA2
        do j=1,air_it
            SA1(i)=SA1(i)+2*dydx(j)*cos( j/real(air_it)*pi )/real(air_it)
            
            SA2(i)=SA2(i)+2*dydx(j)*cos( 2.0*j/real(air_it)*pi )/real(air_it)
        end do
        
        !**计算压力中心**
        SXcp(i)=1.0/4*( 1.0+( SA1(i)-SA2(i) )/2.0/( arfa-SarfaL0(i) ) )
        
        !**计算升力系数**
        
        SCL(i)=2*pi*( arfa-SarfaL0(i) )
        
    
    end do
    
    
    !计算中部机翼零升力攻角以及系数A1、A2
    dydx=0.0
    do i=1,Mwing_num
        !先计算dydx
        call get_dydx( air_it,Mca(i),Mpoint_num(i),Mwing_node(i,:,:),dydx )
        
        !计算零升力攻角
        do j=1,air_it
            MarfaL0(i)=MarfaL0(i)-dydx(j)*(1-cos(j/real(air_it)*pi))/real(air_it)           
        end do
        
        !计算系数MA1,MA2
        do j=1,air_it
            MA1(i)=MA1(i)+2*dydx(j)*cos( j/real(air_it)*pi )/real(air_it)
            
            MA2(i)=MA2(i)+2*dydx(j)*cos( 2.0*j/real(air_it)*pi )/real(air_it)
        end do
        
        !**计算压力中心**
        MXcp(i)=1.0/4*( 1.0+( MA1(i)-MA2(i) )/2.0/( arfa-MarfaL0(i) ) )
        
        !**计算升力系数**
        
        MCL(i)=2*pi*( arfa-MarfaL0(i) )
        
      
    end do
    
    
    !**计算升力**
    allocate( Sair_force(NR),Mair_force(NR) )
    Sair_force=0.0;Mair_force=0.0;
    allocate( Sair_area(Swing_num),Mair_area(Mwing_num) )
    Sair_area=0.0;Mair_area=0.0;
    allocate( temp_Sca(Swing_num),temp_Mca(Mwing_num) )
    temp_Sca=0.0;temp_Mca=0.0;
    allocate( Sair_xx(Swing_num),Mair_xx(Mwing_num) )
    Sair_xx=0.0;Mair_xx=0.0;
    !尾翼    
    do i=1,Swing_num
     
        !计算翼元有效面积
        if( i<Swing_num )then            
            temp_Sca(i)=( Sca(i+1)+Sca(i) )/2.0
            Sair_area(i)=abs( Swing_point(i+1,1,2)-Swing_point(i,1,2) )*temp_Sca(i)/1000.0/1000.0
        else
            temp_Sca(i)=0.0
            Sair_area(i)=0.0
        end if
        
        !计算各翼元剖线到重心的纵向距离
        !按抬头力矩为正
        Sair_xx(i)=( Swing_point(i,Spoint_num(i),1)-Cog0(1)-SXcp(i)*Sca(i) )/1000.0
        
        
        !计算升力
        Sair_force(3)=Sair_force(3)+0.5*SCL(i)*U0*U0*Sair_area(i)
        
        Sair_force(5)=Sair_force(5)+0.5*SCL(i)*U0*U0*Sair_area(i)*Sair_xx(i)
        
    end do
    
    !中部机翼
    do i=1,Mwing_num
        !计算翼元有效面积
        if( i<Mwing_num )then            
            temp_Mca(i)=( Mca(i+1)+Mca(i) )/2.0
            Mair_area(i)=abs( Mwing_point(i+1,1,2)-Mwing_point(i,1,2) )*temp_Mca(i)/1000.0/1000.0
        else
            temp_Mca(i)=0.0
            Mair_area(i)=0.0
        end if
        
        !计算各翼元剖线到重心的纵向距离
        !按抬头力矩为正
        Mair_xx(i)=( Mwing_point(i,Mpoint_num(i),1)-Cog0(1)-MXcp(i)*Mca(i) )/1000.0
        
        !计算升力
        Mair_force(3)=Mair_force(3)+0.5*MCL(i)*U0*U0*Mair_area(i)
        
        Mair_force(5)=Mair_force(5)+0.5*MCL(i)*U0*U0*Mair_area(i)*Mair_xx(i)
        
    end do
    
    
    !计算总气动升力
    air_force=0.0;
    air_force(3)=( Sair_force(3)+Mair_force(3) )*-2.0
    air_force(5)=( Sair_force(5)+Mair_force(5) )*-2.0
    
    
    
    
    
    
    
    !write( 1000,'(2f16.4)' ) air_force(3),air_force(5)
    
    
    deallocate( SCL,SXcp,MCL,MXcp )
    
    deallocate( SarfaL0,MarfaL0 )
    
    deallocate( SA1,SA2,MA1,MA2 )
    
    deallocate( dydx )
    
    deallocate( Sair_force,Mair_force )
    
    deallocate( Sair_area,Mair_area )
    
    deallocate( temp_Sca,temp_Mca )
    
    deallocate( Sair_xx,Mair_xx )
    
end subroutine airforce 