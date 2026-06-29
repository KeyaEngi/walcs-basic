program main
!------------------------------------------------------------------------------------------
!本程序集合WALCS无对称性和含有对称性的组合程序。
!王川  2012.09.16
!------------------------------------------------------------------------------------------
USE DFLIB !?????
USE MESH_MOD,ONLY:NPANEL,SYMTYPE !加钢盖的总数？？？,X or XY
USE READFILE_MOD
USE GETMESH_MOD
USE GETM_MOD
USE INOUTACCESS_MOD,ONLY:PROJNAME
USE ENVIRONMENT_MOD,ONLY:H
USE CAL_MOD,ONLY:WATERDEPTH
USE MAINCAL_MOD
USE PRINT_MOD,ONLY:PRTMOTION

IMPLICIT NONE
!-------------------------------------
INTERFACE !定义子程序的使用接口
	SUBROUTINE	IUTMP(IU)
	INTEGER	::	IU 
	END SUBROUTINE	IUTMP
END INTERFACE
!-------------------------------------
INTEGER(4)::FREQUENCY=300 !频率
INTEGER(4)::DURATION=200  !持续的时间
INTEGER(4) ::  IU_TMP !???
CHARACTER(LEN=14)  ::  RF_TMP
CHARACTER(LEN=2)   ::  NO

LOGICAL::EXF
!-----------------------------------------------------------------
!-----------------------------------------------------------------------------
DO	IU_TMP=1,	99
	CALL QINTTOSTR(IU_TMP,NO,2)			!此时NO为长度为5的字符串
	RF_TMP='WaLCS_NL'//NO//'.TMP'
	INQUIRE(FILE=RF_TMP,EXIST=EXF)
	IF(EXF==.true.)	then
		OPEN(IU_TMP,FILE=RF_TMP,STATUS='OLD',DISPOSE='DELETE')
	endif
END DO
CALL READFILE()
call READMESH()

write(*,*)"三维波浪载荷计算系统COMPASS-WALCS-BASIC正在计算中......"
write(*,*) 
CALL GETMATRIX()

write(*,'(a,a)')" 工程名称: ",projname
write(*,*) 
!by shiyuyun
if(trim(adjustl(SYMTYPE))=="NO")then
write(*,'(a,i4)')" 水动力网格总数目: ",NPanel," 计算水动力网格数目: ",NPanel
else if(trim(adjustl(SYMTYPE))=="X")then
write(*,'(a,i4)')" 水动力网格总数目: ",2*NPanel," 计算水动力网格数目: ",NPanel
else if(trim(adjustl(SYMTYPE))=="XY")then
write(*,'(a,i4)')" 水动力网格总数目: ",4*NPanel," 计算水动力网格数目: ",NPanel
end if
write(*,*) 
IF (trim(adjustl(WATERDEPTH))=="INFINITE")THEN
	write(*,*)"工作水深: 无限水深 "
else
	write(*,"(a,f6.2,a)")" 工作水深: 有限水深，水深为 ",H," 米"
endif
write(*,*) 


CALL MAINCAL()




CALL BEEPQQ(FREQUENCY,DURATION)


pause
!--------------------------------------------------
end program MAIN