
!当前输出文件仅有运动、加速度和剖面载荷，输出序号22-29（20230418）
    
    
subroutine outputfile( load_name,load_ID )
    use Constant,only:OutAccess,projname,wavectrl,U0,head,Pi,amp,ome,hs_num,tz_num,Slamctrl
    use ShipHullVar,only:NBSECT
    !20231115修改
    use Slamming,only:SlamRelaMctr
    
    implicit none
    
    integer(4)::load_ID
    character(len=4)::load_name
    
    integer(4)::i,j,k
    
    character(len=1)::tampname1
    character(len=2)::tampname2
    character(len=3)::tampname3
    character(len=20)::sectname
    
    
    
!****以上是变量定义****
    
    if(trim(adjustl(wavectrl))=='Calm' ) then
        !输出摇荡位移
        open(23,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".mst")
        write(23,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(23,"(A)") '#          [运动响应时历输出文件]'
        write(23,"(A)") '#波浪型式'
        write(23,"(A)") wavectrl
        write(23,"(A)") '#航速(kn)'
        write(23,"(f12.3)") U0/0.5144
        write(23,"(A)") '#T,Surge(m),Sway(m),Heave(m),Roll(rad),Pitch(rad),Yaw(rad)'
        
        !输出摇荡加速度
        open(22,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".acc")
        write(22,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(22,"(A)") '#          [加速度时历输出文件]'
        write(22,"(A)") '#波浪型式'
        write(22,"(A)") wavectrl
        write(22,"(A)") '#航速(kn)'
        write(22,"(f12.3)") U0/0.5144
        write(22,"(A)") '#T,Surge(m/s2),Sway(m/s2),Heave(m/s2),Roll(rad/s2),Pitch(rad/s2),Yaw(rad/s2)'
        
        !输出剖面载荷
        open(24,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fx") !剖面载荷
        open(25,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fy")
        open(26,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fz")
        open(27,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Mx")
        open(28,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".My")
        open(29,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Mz")
        !fx
        write(24,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(24,"(A)") '#          [剖面载荷Fx时历输出文件]'
        write(24,"(A)") '#波浪型式'
        write(24,"(A)") wavectrl
        write(24,"(A)") '#航速(kn)'
        write(24,"(f12.3)") U0/0.5144
        write(24,"(A)") '#计算剖面数目'
        write(24,"(i8)") NBSECT
        
        write(24,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fx_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fx_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fx_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(24,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(24,"(/)")
        !fy
        write(25,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(25,"(A)") '#          [剖面载荷Fy时历输出文件]'
        write(25,"(A)") '#波浪型式'
        write(25,"(A)") wavectrl
        write(25,"(A)") '#航速(kn)'
        write(25,"(f12.3)") U0/0.5144
        write(25,"(A)") '#计算剖面数目'
        write(25,"(i8)") NBSECT
        
        write(25,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fy_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fy_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fy_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(25,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(25,"(/)")
        !fz
        write(26,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(26,"(A)") '#          [剖面载荷Fz时历输出文件]'
        write(26,"(A)") '#波浪型式'
        write(26,"(A)") wavectrl
        write(26,"(A)") '#航速(kn)'
        write(26,"(f12.3)") U0/0.5144
        write(26,"(A)") '#计算剖面数目'
        write(26,"(i8)") NBSECT
        
        write(26,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fz_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fz_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fz_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(26,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(26,"(/)")
        !mx
        write(27,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(27,"(A)") '#          [剖面载荷Mx时历输出文件]'
        write(27,"(A)") '#波浪型式'
        write(27,"(A)") wavectrl
        write(27,"(A)") '#航速(kn)'
        write(27,"(f12.3)") U0/0.5144
        write(27,"(A)") '#计算剖面数目'
        write(27,"(i8)") NBSECT
        
        write(27,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Mx_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Mx_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Mx_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(27,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(27,"(/)")
        !my
        write(28,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(28,"(A)") '#          [剖面载荷My时历输出文件]'
        write(28,"(A)") '#波浪型式'
        write(28,"(A)") wavectrl
        write(28,"(A)") '#航速(kn)'
        write(28,"(f12.3)") U0/0.5144
        write(28,"(A)") '#计算剖面数目'
        write(28,"(i8)") NBSECT
        
        write(28,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='My_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='My_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='My_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(28,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(28,"(/)")
        !mz
        write(29,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(29,"(A)") '#          [剖面载荷Mz时历输出文件]'
        write(29,"(A)") '#波浪型式'
        write(29,"(A)") wavectrl
        write(29,"(A)") '#航速(kn)'
        write(29,"(f12.3)") U0/0.5144
        write(29,"(A)") '#计算剖面数目'
        write(29,"(i8)") NBSECT
        
        write(29,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Mz_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Mz_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Mz_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(29,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(29,"(/)")
        
      
        
        
        !20231116修改，输出单条剖线上的砰击载荷
        if( trim(adjustl(Slamctrl))=='YES' ) then
            open(65,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamlineForce.Fx")
            open(66,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamlineForce.Fz")      
        end if        
        
        
        !20231115修改，添加船中纵剖线各点与波面的相对位置计算
        if( trim(adjustl(SlamRelaMctr))=='YES' )then
            open(67,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamShipLinePx.Dx")
            open(68,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamShipLinePz.Dz")
            open(69,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamWaveLinePx.Dx")
            open(70,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamWaveLinePz.Dz")
            
            open(71,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"Wave.txt")
        end if
        
        
        
        
        
    !20231113修改，目的是后续添加不规则模块    
    elseif( trim(adjustl(wavectrl))=='Regular' )then
        !输出摇荡位移
        open(23,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".mst")
        write(23,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(23,"(A)") '#          [运动响应时历输出文件]'
        write(23,"(A)") '#波浪型式'
        write(23,"(A)") wavectrl
        write(23,"(A)") '#航速(kn)'
        write(23,"(f12.3)") U0/0.5144
        write(23,"(A)") '#浪向角(°'
        write(23,"(f12.3)") 180.0-head/Pi*180.0
        write(23,"(A)") '#波幅(m)'
        write(23,"(f12.3)") amp
        write(23,"(A)") '#波浪频率(rad/s)'
        write(23,"(f12.3)") ome
        write(23,"(A)") '#T,Surge(m),Sway(m),Heave(m),Roll(rad),Pitch(rad),Yaw(rad)'
        
        !输出摇荡加速度
        open(22,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".acc") 
        write(22,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(22,"(A)") '#          [加速度时历输出文件]'
        write(22,"(A)") '#波浪型式'
        write(22,"(A)") wavectrl
        write(22,"(A)") '#航速(kn)'
        write(22,"(f12.3)") U0/0.5144
        write(22,"(A)") '#浪向角(°'
        write(22,"(f12.3)") 180.0-head/Pi*180.0
        write(22,"(A)") '#波幅(m)'
        write(22,"(f12.3)") amp
        write(22,"(A)") '#波浪频率(rad/s)'
        write(22,"(f12.3)") ome
        write(22,"(A)") '#T,Surge(m/s2),Sway(m/s2),Heave(m/s2),Roll(rad/s2),Pitch(rad/s2),Yaw(rad/s2)'
        
        !输出剖面载荷
        open(24,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fx") !剖面载荷
        open(25,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fy")
        open(26,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fz")
        open(27,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Mx")
        open(28,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".My")
        open(29,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Mz")
        !fx
        write(24,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(24,"(A)") '#          [剖面载荷Fx时历输出文件]'
        write(24,"(A)") '#波浪型式'
        write(24,"(A)") wavectrl
        write(24,"(A)") '#航速(kn)'
        write(24,"(f12.3)") U0/0.5144
        write(24,"(A)") '#浪向角(°'
        write(24,"(f12.3)") 180.0-head/Pi*180.0
        write(24,"(A)") '#波幅(m)'
        write(24,"(f12.3)") amp
        write(24,"(A)") '#波浪频率(rad/s)'
        write(24,"(f12.3)") ome
        write(24,"(A)") '#计算剖面数目'
        write(24,"(i8)") NBSECT
        
        write(24,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fx_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fx_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fx_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(24,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(24,"(/)")
        !fy
        write(25,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(25,"(A)") '#          [剖面载荷Fy时历输出文件]'
        write(25,"(A)") '#波浪型式'
        write(25,"(A)") wavectrl
        write(25,"(A)") '#航速(kn)'
        write(25,"(f12.3)") U0/0.5144
        write(25,"(A)") '#浪向角(°'
        write(25,"(f12.3)") 180.0-head/Pi*180.0
        write(25,"(A)") '#波幅(m)'
        write(25,"(f12.3)") amp
        write(25,"(A)") '#波浪频率(rad/s)'
        write(25,"(f12.3)") ome
        write(25,"(A)") '#计算剖面数目'
        write(25,"(i8)") NBSECT
        
        write(25,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fy_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fy_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fy_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(25,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(25,"(/)")
        !fz
        write(26,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(26,"(A)") '#          [剖面载荷Fz时历输出文件]'
        write(26,"(A)") '#波浪型式'
        write(26,"(A)") wavectrl
        write(26,"(A)") '#航速(kn)'
        write(26,"(f12.3)") U0/0.5144
        write(26,"(A)") '#浪向角(°'
        write(26,"(f12.3)") 180.0-head/Pi*180.0
        write(26,"(A)") '#波幅(m)'
        write(26,"(f12.3)") amp
        write(26,"(A)") '#波浪频率(rad/s)'
        write(26,"(f12.3)") ome
        write(26,"(A)") '#计算剖面数目'
        write(26,"(i8)") NBSECT
        
        write(26,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fz_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fz_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fz_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(26,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(26,"(/)")
        !mx
        write(27,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(27,"(A)") '#          [剖面载荷Mx时历输出文件]'
        write(27,"(A)") '#波浪型式'
        write(27,"(A)") wavectrl
        write(27,"(A)") '#航速(kn)'
        write(27,"(f12.3)") U0/0.5144
        write(27,"(A)") '#浪向角(°'
        write(27,"(f12.3)") 180.0-head/Pi*180.0
        write(27,"(A)") '#波幅(m)'
        write(27,"(f12.3)") amp
        write(27,"(A)") '#波浪频率(rad/s)'
        write(27,"(f12.3)") ome
        write(27,"(A)") '#计算剖面数目'
        write(27,"(i8)") NBSECT
        
        write(27,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Mx_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Mx_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Mx_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(27,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(27,"(/)")
        !my
        write(28,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(28,"(A)") '#          [剖面载荷My时历输出文件]'
        write(28,"(A)") '#波浪型式'
        write(28,"(A)") wavectrl
        write(28,"(A)") '#航速(kn)'
        write(28,"(f12.3)") U0/0.5144
        write(28,"(A)") '#浪向角(°'
        write(28,"(f12.3)") 180.0-head/Pi*180.0
        write(28,"(A)") '#波幅(m)'
        write(28,"(f12.3)") amp
        write(28,"(A)") '#波浪频率(rad/s)'
        write(28,"(f12.3)") ome
        write(28,"(A)") '#计算剖面数目'
        write(28,"(i8)") NBSECT
        
        write(28,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='My_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='My_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='My_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(28,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(28,"(/)")
        !mz
        write(29,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(29,"(A)") '#          [剖面载荷Mz时历输出文件]'
        write(29,"(A)") '#波浪型式'
        write(29,"(A)") wavectrl
        write(29,"(A)") '#航速(kn)'
        write(29,"(f12.3)") U0/0.5144
        write(29,"(A)") '#浪向角(°'
        write(29,"(f12.3)") 180.0-head/Pi*180.0
        write(29,"(A)") '#波幅(m)'
        write(29,"(f12.3)") amp
        write(29,"(A)") '#波浪频率(rad/s)'
        write(29,"(f12.3)") ome
        write(29,"(A)") '#计算剖面数目'
        write(29,"(i8)") NBSECT
        
        write(29,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Mz_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Mz_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Mz_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(29,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(29,"(/)")
        
        
        
        !20231116修改，输出单条剖线上的砰击载荷
        if( trim(adjustl(Slamctrl))=='YES' ) then
            open(65,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamlineForce.Fx")
            open(66,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamlineForce.Fz")      
        end if
        
        
        !20231115修改，添加船中纵剖线各点与波面的相对位置计算
        if( trim(adjustl(SlamRelaMctr))=='YES' )then
            open(67,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamShipLinePx.Dx")
            open(68,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamShipLinePz.Dz")
            open(69,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamWaveLinePx.Dx")
            open(70,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamWaveLinePz.Dz")
            
            open(71,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"Wave.txt")
        end if
        
        
        
    !20231113修改    
    elseif( trim(adjustl(wavectrl))=='Irregular' )then
        
        !输出摇荡运动(位移)
        open(23,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".mst") 
        write(23,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(23,"(A)") '#          [运动响应时历输出文件]'
        write(23,"(A)") '#波浪型式'
        write(23,"(A)") wavectrl
        write(23,"(A)") '#航速(kn)'
        write(23,"(f12.3)") U0/0.5144
        write(23,"(A)") '#浪向角(°'
        write(23,"(f12.3)") 180.0-head/Pi*180.0
        write(23,"(A)") '#有义波高(m)'
        write(23,"(f12.3)") hs_num(load_ID)
        write(23,"(A)") '#平均跨零周期(s)'
        write(23,"(f12.3)") tz_num(load_ID)
        write(23,"(A)") '#T,Surge(m),Sway(m),Heave(m),Roll(rad),Pitch(rad),Yaw(rad)'
        !输出摇荡加速度
        open(22,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".acc") 
        write(22,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(22,"(A)") '#          [加速度时历输出文件]'
        write(22,"(A)") '#波浪型式'
        write(22,"(A)") wavectrl
        write(22,"(A)") '#航速(kn)'
        write(22,"(f12.3)") U0/0.5144
        write(22,"(A)") '#浪向角(°'
        write(22,"(f12.3)") 180.0-head/Pi*180.0
        write(22,"(A)") '#有义波高(m)'
        write(22,"(f12.3)") hs_num(load_ID)
        write(22,"(A)") '#平均跨零周期(s)'
        write(22,"(f12.3)") tz_num(load_ID)
        write(22,"(A)") '#T,Surge(m/s2),Sway(m/s2),Heave(m/s2),Roll(rad/s2),Pitch(rad/s2),Yaw(rad/s2)'
        
        !输出剖面载荷
        open(24,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fx") !剖面载荷
        open(25,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fy")
        open(26,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Fz")
        open(27,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Mx")
        open(28,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".My")
        open(29,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//".Mz")
        !fx
        write(24,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(24,"(A)") '#          [剖面载荷Fx时历输出文件]'
        write(24,"(A)") '#波浪型式'
        write(24,"(A)") wavectrl
        write(24,"(A)") '#航速(kn)'
        write(24,"(f12.3)") U0/0.5144
        write(24,"(A)") '#浪向角(°'
        write(24,"(f12.3)") 180.0-head/Pi*180.0
        write(24,"(A)") '#有义波高(m)'
        write(24,"(f12.3)") hs_num(load_ID)
        write(24,"(A)") '#平均跨零周期(s)'
        write(24,"(f12.3)") tz_num(load_ID)
        write(24,"(A)") '#计算剖面数目'
        write(24,"(i8)") NBSECT
        
        write(24,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fx_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fx_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fx_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(24,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(24,"(/)")
        !fy
        write(25,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(25,"(A)") '#          [剖面载荷Fy时历输出文件]'
        write(25,"(A)") '#波浪型式'
        write(25,"(A)") wavectrl
        write(25,"(A)") '#航速(kn)'
        write(25,"(f12.3)") U0/0.5144
        write(25,"(A)") '#浪向角(°'
        write(25,"(f12.3)") 180.0-head/Pi*180.0
        write(25,"(A)") '#有义波高(m)'
        write(25,"(f12.3)") hs_num(load_ID)
        write(25,"(A)") '#平均跨零周期(s)'
        write(25,"(f12.3)") tz_num(load_ID)
        write(25,"(A)") '#计算剖面数目'
        write(25,"(i8)") NBSECT

        write(25,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fy_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fy_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fy_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(25,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(25,"(/)")
        !fz
        write(26,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(26,"(A)") '#          [剖面载荷Fz时历输出文件]'
        write(26,"(A)") '#波浪型式'
        write(26,"(A)") wavectrl
        write(26,"(A)") '#航速(kn)'
        write(26,"(f12.3)") U0/0.5144
        write(26,"(A)") '#浪向角(°'
        write(26,"(f12.3)") 180.0-head/Pi*180.0
        write(26,"(A)") '#有义波高(m)'
        write(26,"(f12.3)") hs_num(load_ID)
        write(26,"(A)") '#平均跨零周期(s)'
        write(26,"(f12.3)") tz_num(load_ID)
        write(26,"(A)") '#计算剖面数目'
        write(26,"(i8)") NBSECT

        write(26,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Fz_SEC'//trim(adjustl(tampname1))//'(kN)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Fz_SEC'//trim(adjustl(tampname2))//'(kN)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Fz_SEC'//trim(adjustl(tampname3))//'(kN)'
            end if
            write(26,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(26,"(/)")
        !mx
        write(27,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(27,"(A)") '#          [剖面载荷Mx时历输出文件]'
        write(27,"(A)") '#波浪型式'
        write(27,"(A)") wavectrl
        write(27,"(A)") '#航速(kn)'
        write(27,"(f12.3)") U0/0.5144
        write(27,"(A)") '#浪向角(°'
        write(27,"(f12.3)") 180.0-head/Pi*180.0
        write(27,"(A)") '#有义波高(m)'
        write(27,"(f12.3)") hs_num(load_ID)
        write(27,"(A)") '#平均跨零周期(s)'
        write(27,"(f12.3)") tz_num(load_ID)
        write(27,"(A)") '#计算剖面数目'
        write(27,"(i8)") NBSECT

        write(27,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Mx_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Mx_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Mx_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(27,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(27,"(/)")
        
        !my
        write(28,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(28,"(A)") '#          [剖面载荷My时历输出文件]'
        write(28,"(A)") '#波浪型式'
        write(28,"(A)") wavectrl
        write(28,"(A)") '#航速(kn)'
        write(28,"(f12.3)") U0/0.5144
        write(28,"(A)") '#浪向角(°'
        write(28,"(f12.3)") 180.0-head/Pi*180.0
        write(28,"(A)") '#有义波高(m)'
        write(28,"(f12.3)") hs_num(load_ID)
        write(28,"(A)") '#平均跨零周期(s)'
        write(28,"(f12.3)") tz_num(load_ID)
        write(28,"(A)") '#计算剖面数目'
        write(28,"(i8)") NBSECT
 
        write(28,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='My_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='My_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='My_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(28,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(28,"(/)")
        
        !mz
        write(29,"(A)") '#[ COMPASS-WALCS-NL V1.0 ]'
        write(29,"(A)") '#          [剖面载荷Mz时历输出文件]'
        write(29,"(A)") '#波浪型式'
        write(29,"(A)") wavectrl
        write(29,"(A)") '#航速(kn)'
        write(29,"(f12.3)") U0/0.5144
        write(29,"(A)") '#浪向角(°'
        write(29,"(f12.3)") 180.0-head/Pi*180.0
        write(29,"(A)") '#有义波高(m)'
        write(29,"(f12.3)") hs_num(load_ID)
        write(29,"(A)") '#平均跨零周期(s)'
        write(29,"(f12.3)") tz_num(load_ID)
        write(29,"(A)") '#计算剖面数目'
        write(29,"(i8)") NBSECT

        write(29,"(A,\)",advance='NO') '#   T(s)     '
        do i=1,NBSECT
            if(i<10) then
                call qinttostr(i,tampname1,1)
                sectname='Mz_SEC'//trim(adjustl(tampname1))//'(kN.m)'
            elseif(i<100) then
                call qinttostr(i,tampname2,2)
                sectname='Mz_SEC'//trim(adjustl(tampname2))//'(kN.m)'
            elseif(i<1000) then
                call qinttostr(i,tampname3,3)
                sectname='Mz_SEC'//trim(adjustl(tampname3))//'(kN.m)'
            end if
            write(29,"(A,4x,\)",advance='NO') trim(adjustl(sectname))
        end do
        write(29,"(/)")
        


        !20231116修改，输出单条剖线上的砰击载荷
        if( trim(adjustl(Slamctrl))=='YES' ) then
            open(65,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamlineForce.Fx")
            open(66,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamlineForce.Fz")      
        end if
        
        !20231115修改，添加船中纵剖线各点与波面的相对位置计算
        if( trim(adjustl(SlamRelaMctr))=='YES' )then
            open(67,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamShipLinePx.Dx")
            open(68,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamShipLinePz.Dz")
            open(69,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamWaveLinePx.Dx")
            open(70,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"SlamWaveLinePz.Dz")
            
            open(71,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(projname))//"_LC"//trim(adjustl(load_name))//"Wave.txt")
        end if
        
        
        

        
    end if
    
    
    
    
    
end subroutine outputfile