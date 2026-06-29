!这个模块用来实现CPU时间的计算。可用于不同程序的测试比较。

Module M_CPUTime

USE IFPORT

Implicit None

  Integer(2),Save::StartTime(4),EndTime(4),CPUTime(4)
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
      write(*,'(2x,3(i2,a1),i2)') CPUTime(1),':',CPUTime(2),':',CPUTime(3),':',CPUTime(4) 
   EndSubroutine    

End Module M_CPUTime
