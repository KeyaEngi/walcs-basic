!这个模块用来实现CPU时间的计算。可用于不同程序的测试比较。
Module CPUTime_MOD

USE IFPORT

Implicit None

integer(2),Save::StartTime(4),EndTime(4),CPUTime(4)
Private StartTime,EndTime

Contains

   SubRoutine GetStartTime()
      call gettim(StartTime(1),StartTime(2),StartTime(3),StartTime(4))
   EndSubroutine    

   Subroutine GetEndTime()
     Integer(2) i
       call gettim(EndTime(1),EndTime(2),EndTime(3),EndTime(4))
       do i=1,4
        if(Endtime(i)>=StartTime(i))then
            CPUTime(i)=Endtime(i)-StartTime(i)
        else
            CPUTime(i-1)=CPUTime(i-1)-1
            if(i/=4)then
                CPUTime(i)=Endtime(i)+60-Starttime(i)
            else
                CPUTime(i)=EndTime(i)+100-StartTime(i)
            endif
        endif
      EndDo
	  CPUTIME(4)= CPUTIME(1)*3600+CPUTIME(2)*60+CPUTIME(3)  !---------
	  write(*,*)
      write(*,'(A,2x,i4,A)')" 所用时间:",CPUTime(4)," s" !---------
   EndSubroutine    

End Module CPUTime_MOD

!Use
!  M_CPUTime
!...  
!Call GetStartTime()
!...
!Call GetEndTime()