!程序名称：三维时域高速水面飞行器水弹性分析软件
!功能：用于分析高速水面飞行器处于滑行或半滑行状态的运动及剖面载荷
!程序框架搭建：李创20230516
!特别感谢：邹健、路琳、邓保利
  
    
program main
    use SAILPARAM_MOD, only: SAIL_OK
    use SAILPLANNING_ADAPTER_MOD, only: InitializeSailPlanningAdapter, &
        FinalizeSailPlanningAdapter
    
    implicit none
    
    real*8 :: time_begin,time_end
    integer :: sail_ierr
    character(len=2048) :: sail_message
    character(len=*), parameter :: sail_database_file = &
        'D:\GitHub\walcs-basic\SailModule\sail_database.dat'
    
    call CPU_TIME(time_begin)
    
    
!文件读取    
    call readfile

!风帆模块数据库在进入时间积分之前只初始化一次
    call InitializeSailPlanningAdapter(sail_database_file, sail_ierr, sail_message)
    if (sail_ierr /= SAIL_OK) then
        write(*,'(A,I0,2A)') 'Sail adapter initialization failed, IERR=', &
            sail_ierr, ': ', trim(sail_message)
        error stop 1
    end if
    
!时域模拟    
    call Timmotion

    call FinalizeSailPlanningAdapter(sail_ierr, sail_message)
    if (sail_ierr /= SAIL_OK) then
        write(*,'(A,I0,2A)') 'Sail adapter finalization failed, IERR=', &
            sail_ierr, ': ', trim(sail_message)
        error stop 1
    end if
    
    
    call CPU_TIME(time_end)

    write(*,*) " 完整的计算时长是 ", time_end-time_begin ,"Seconds"   

    
    
    
end program main
