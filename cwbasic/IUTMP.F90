!******************************************************************************************************************
!用于略去输入数据文件中的#注释内容:#必须在第一行，注释内容以字符"#"打头,该行在#后的内容都被认为是注释内容
!默认读取数据的列数不能大于150，否则就要加大字符串C的声明长度
!程序编制： 单鹏昊
!******************************************************************************************************************
SUBROUTINE IUTMP(IU)     ! 完全版修改版IUTMP  
IMPLICIT NONE 
INTEGER::IU                             !输入文件号
!-----以下说明的变量仅用于本子程序段--------------
INTEGER::IU_TMP                     !临时文件号
LOGICAL::EX,EXF
CHARACTER(LEN=14)::RF_TMP    !临时文件名
CHARACTER(LEN=2)::NO
CHARACTER(LEN=500)::C
INTEGER::PP
!------以下为判断文件号IU_TMP是否已经被使用的代码-----------------------------
IU_TMP=1
DO
	INQUIRE(IU_TMP,OPENED=EX)    !判断文件号IU_TMP是否已经被使用
	CALL QINTTOSTR(IU_TMP,NO,2)    !此时NO为长度为2的字符串
	RF_TMP='WaLCS_NL'//NO//'.TMP'
	INQUIRE(FILE=RF_TMP,EXIST=EXF)
	IF(EX==.FALSE..AND.EXF==.FALSE.)EXIT
	IU_TMP=IU_TMP+1
END DO
!------以上为判断文件号IU_TMP是否已经被使用的代码-----------------------------
OPEN(IU_TMP,FILE=RF_TMP)        !RF_TMP用于临时保存不含注释的文件内容
DO WHILE( .TRUE. )
	    READ(IU,'(A)',IOSTAT=PP)C   
	    IF  (PP/=0) EXIT                     !判断文件是否正常读取内容，否，则退出
	    IF (C /=" ") THEN                   !如果读到的不是空行
               IF(C(1:1)=='#')THEN      
                    CYCLE
               ELSE
                   WRITE(IU_TMP,'(A)')C
               END IF
	    ELSE                                        !空行
    	
	    END IF
END DO
CLOSE(IU);CLOSE(IU_TMP)
OPEN(IU,FILE=RF_TMP,STATUS='OLD',DISPOSE='DELETE')  !输入文件号重定向数据文件，文件调用完删除
END SUBROUTINE IUTMP
!--------------------------------------------------------------------------------------
SUBROUTINE QINTTOSTR(M,STR,N)
!功    能：利用文件号创建临时的文件名;
!将一个整型数据转换为包含N个字符的字符串(取末尾的几个数字)，1--01；2--02；...8--08
!使用说明：N不能超过10(若需要更长的字符串可通过调整程序中数组INC的大小实现)
IMPLICIT NONE
INTEGER::M,N
CHARACTER(LEN=N)::STR
!-----以下说明的变量仅用于本子程序段--------------
INTEGER::N1,I
INTEGER,DIMENSION(1:10)::INC
!-----以下为程序代码------------------------------
N1=M
DO I=1,N
	INC(I)=N1-INT(N1/10)*10
	N1=INT(N1/10)
END DO
STR=CHAR(48+INC(N))
DO I=2,N
	STR=TRIM(STR)//CHAR(48+INC(N-I+1))
END DO
END SUBROUTINE QINTTOSTR
!-----------------------------IUTMP.F90/(2001/11)----------------------------
