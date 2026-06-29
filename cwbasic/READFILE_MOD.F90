!*********************************************************************************************************************
!读取路径变量(ACCESS.DAT)、总体设置文件(config.hyp)、计算控制参数(ProjName.cfg)、环境参数参数(ProjName.dhc)
!*********************************************************************************************************************
MODULE  READFILE_MOD

USE INOUTACCESS_MOD
USE ENVIRONMENT_MOD

USE CAL_MOD
USE PRINT_MOD

IMPLICIT NONE
PRIVATE
PUBLIC::READFILE

CONTAINS 
!-------------------------------------------------------------------------------------
SUBROUTINE READFILE()    
IMPLICIT NONE

character(LEN=300)	::	InputAccessHYP
CHARACTER(LEN=300)::CFGINPUT,DHCINPUT,CHEOUTPUT
INTEGER::I
!---------------------读取路径变量(ACCESS.DAT)---------------------------
open(1,file='Access.dat',STATUS='OLD')
call IUTMP(1)
READ(1,'(A)')	InAccess
READ(1,'(A)')	OutAccess
close(1)
!----------Open Config File总体设置文件(config.hyp)---------------------
InputAccessHYP=trim(adjustl(InAccess))//'\'//'config.hyp'
open(10,file=InputAccessHYP,status='old')
CALL IUTMP(10)
READ(10,*)	projname !工程名称
close(10)
!-------------------------读入计算控制参数(ProjName.cfg)-------------------------
CFGINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.cfg'
OPEN (20,FILE=CFGINPUT,STATUS="OLD") !20=cfg 10=hyp
CALL IUTMP(20)
!读取输入控制参数
READ(20,"(A)") WATERDEPTH !无限水深or有限水深
READ(20,"(A)") INPUTSTATIC !船舶静力特性
!READ(20,"(A)") INPUTINERTIA
READ(20,"(A)") MoorLine !锚链效应
READ(20,*)     ROllDAMPING  !横摇阻尼
READ(20,"(A)") HULLTYPE !船体类型：单体or多体
READ(20,"(A)") OMETYPE !波浪频率输入方式：自然频率or周期

!READ(20,"(A)") MATRIXSOLVE
! READ(20,"(A)") WAVEFI
! READ(20,*) POINTPRESS
!READ(20,"(A)") ModiDeflect

!读取打印输出控制参数
READ(20,"(A)") PRTSTATIC
READ(20,"(A)") PRTMOTION
READ(20,"(A)") PRTpPRESS
READ(20,"(A)") PRTSECTLOAD
READ(20,"(A)") PRTva
READ(20,"(A)") ansMethod
READ(20,"(A)") PRTFEPRESS
!READ(20,"(A)") FLOATATION
!READ(20,"(A)") BILGEFORCE
CLOSE(20)
!--------------------读取环境参数参数(ProjName.dhc)--------------------------------------
NBU=1;NBHEAD=1;NBOME=1

DHCINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.dhc'
OPEN (30,FILE=DHCINPUT,STATUS="OLD") !30=dhc
CALL IUTMP(30)

READ(30,*)ROU,H !水密度 水深
READ(30,*)NBHEAD,NBOME !浪向角数目，规则波数目（频率数目）

ALLOCATE(UKN(1:NBU),HEAD(1:NBHEAD),OME0(0:NBOME),AMP(0:NBOME),PHASE(0:NBOME),OME(0:NBOME),OME1(0:NBOME))

READ(30,*)(UKN(i),i=1,NBU)		  !读入航速
READ(30,*)(HEAD(i),i=1,NBHEAD)	  !读入浪向角
READ(30,*)(OME0(i),i=1,NBOME)	  !读入波浪频率
if(TRIM(ADJUSTL(ansMethod))=="YES")	then
    READ(30,*)(AMP(i),i=1,NBOME)	  !读入波幅
    READ(30,*)(PHASE(i),i=1,NBOME)    !读入相位
else
    AMP=1
    PHASE=0
end if
CLOSE(30)

if(trim(adjustl(OMETYPE))=="PERIOD")	then	   !若输入为波浪周期时，换算为波浪频率
      OME0(:)=2.0*PI/OME0(:)
ENDIF
OME(0)=0.15
OME1(0)=0.15
OME0(0)=0.15

!---------------------------check file---------------------------------------------
CHEOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.che'
OPEN(12,FILE=CHEOUTPUT)
write(12,'(a)')'#              [ WALCS Version 2.0 ]'
write(12,*)
write(12,'(a)')'#              [ 输入数据验证文件 ]'
write(12,*)
write(12,'(a)')'#[ 控制参数输入文件 ]'
write(12,*)
if(trim(adjustl(INPUTSTATIC))=="YES")THEN
write(12,'(a)')"  输入船舶静力特性"
else
write(12,"(a)")" 不输入船舶静力特性"
endif
write(12,*)


if(trim(adjustl(MoorLine))=="YES")THEN
write(12,'(a)')"  考虑锚链效应"
else
write(12,"(a)")" 不考虑锚链效应"
endif
write(12,*)

select case(ROllDAMPING)
case(0)
write(12,'(a)')" 横摇修正方法：横摇不修正"
case(1)
write(12,'(a)')" 横摇修正方法：Rma,Rmb法"
case(2)
write(12,'(a)')" 横摇修正方法：Rma,Rmb,ek0法"
case(3)
write(12,'(a)')" 横摇修正方法：米勒法"
case(4)
write(12,'(a)')" 横摇修正方法：驳船经验公式"
case(5)
write(12,'(a)')" 横摇修正方法：临界阻尼法"
endselect
write(12,*)""

if(trim(adjustl(HULLTYPE))=="SINGLE")THEN
write(12,'(a)')" 单体"
else
write(12,"(a)")" 双体"
endif
write(12,*)""

if(trim(adjustl(OMETYPE))=="FREQUENCY")THEN
write(12,'(a)')" 自然频率"
else
write(12,"(a)")" 周期"
endif
write(12,*)""
write(12,"(a)")"#[ 输出控制参数 ]"
write(12,*)""
if(trim(adjustl(PRTSTATIC))=="YES")THEN
write(12,'(a)')" 输出静水载荷计算结果"
else
write(12,"(a)")" 不输出静水载荷计算结果"
endif
write(12,*)

if(trim(adjustl(PRTMOTION))=="YES")THEN
write(12,'(a)')" 输出运动响应计算结果"
else
write(12,"(a)")" 不输出运动响应计算结果"
endif
write(12,*)

if(trim(adjustl(PRTpPRESS))=="YES")THEN
write(12,'(a)')" 输出压力计算点响应结果"
else
write(12,"(a)")" 不输出压力计算点响应结果"
endif
write(12,*)""

if(trim(adjustl(PRTSECTLOAD))=="YES")THEN
write(12,'(a)')" 输出剖面载荷响应结果"
else
write(12,"(a)")" 不输出剖面载荷响应结果"
endif
write(12,*)""

if(trim(adjustl(PRTva))=="YES")THEN
write(12,'(a)')" 输出加速度点响应结果"
else
write(12,"(a)")" 不输出加速度点响应结果"
endif
write(12,*)
if(trim(adjustl(ansMethod))=="NO")THEN
write(12,'(a)')" 谱分析方法"
else
write(12,"(a)")" 设计波方法"
endif
write(12,*)

write(12,"(A)")"#------------------------------------------------------------"
write(12,*)
write(12,'(a)')"#[ 波浪参数输入文件 ] "
write(12,*)
write(12,'(a)')"[海水密度(ton/m^3),工作水深(m)]"
write(12,*)
write(12,'(3f8.3)')rou,g0,H
write(12,*)
write(12,'(a)')"[浪向角数目,波浪频率数目]"
write(12,*)
write(12,'(3i8)') NBHEAD,NBOME
write(12,*)
write(12,"(A)")"#------------------------------------------------------------"
write(12,*)
!--------------------------------------------------------------------------------

END SUBROUTINE READFILE
!------------------------------------------------------------------------------------------
ENDMODULE  READFILE_MOD 
!**************************************************************************************************