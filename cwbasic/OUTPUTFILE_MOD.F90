MODULE OUTPUTFILE_MOD

USE INOUTACCESS_MOD
USE ENVIRONMENT_MOD,ONLY:NBU,NBHEAD,NBOME,UKN,OME0,HEAD,NBSECT,g0,OME,PI,AMP,PHASE,OME1
USE MESH_MOD,ONLY:NPANEL,SL
USE CAL_MOD,ONLY:MASSSOLVE,WATERDEPTH
USE PRINT_MOD,ONLY:PRTSECTLOAD,ANSMETHOD,PRTMOTION,PRTva,PRTPPRESS,BILGEFORCE
USE PVPOINT_MOD,ONLY:NPVA,NPPOINT,Npvf
USE OUTPUT_MOD

IMPLICIT NONE
PRIVATE
PUBLIC::OUTPUTFILE,m_putouts


CONTAINS
!----------------------------------------------------------------------
SUBROUTINE OUTPUTFILE() 
IMPLICIT NONE

INTEGER::I
!------------------------------------------------------------------------------------------
OUTPUT3HC=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.3hc'
OPEN (5,FILE=OUTPUT3HC)
write(5,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
write(5,*)
write(5,'(a)')"#              [ 三维水动力系数结果输出文件 ]"
write(5,*)
write(5,'(a)')"#[工程名称]"
write(5,*)
write(5,'(a)')projname
write(5,*)
write(5,'(a)')"#[水深参数]"
write(5,*)
if(trim(adjustl(WATERDEPTH))=="INFINITE")write(5,'(a)')"无限水深"
write(5,*)
write(5,'(a)')"#[水动力网格数目]"
write(5,*)
write(5,'(i5)')NPanel
write(5,*)
write(5,'(a)')'#[浪向角数目,波浪频率数目]'
write(5,*)
write(5,'(2(i4,1x))')NBHEAD,NBOME
write(5,*)
write(5,'(a)')"#[航速(kn)]"
write(5,*)
write(5,'(10f8.3)')(UKN(i),i=1,NBU)
write(5,*)
write(5,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
write(5,*)
write(5,'(50f8.2)')(HEAD(i),i=1,NBHEAD)
write(5,*)
write(5,'(a)')"#[频率(rad/s)Ome0(i),i=1,NBOME]"
write(5,*)
write(5,'(100f8.3)')(Ome0(i),i=1,NBOME)
write(5,*)
!-----------------------------------------------------------------------------------------------------
if(trim(adjustl(ANSMETHOD))=="NO")	then
    SLTOUTPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.slt'
 	open(33,file=SLTOUTPUT) 
    write(33,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
    write(33,*)
	write(33,'(a)')"#              [ 用于长、短期统计预报的输入文件]"
	write(33,*)   	
	write(33,'(a)')"#[重力加速度,垂线间长,航速数目,浪向角数目,自然频率数目]" 
	write(33,*)
	write(33,'(2x,f8.2,2x,f8.3,2x,I3,2x,I3,2x,I3)')g0,sl,NBU,NBHEAD,NBOME
endif
!------------------------------------------------------------------------------------------------------
if(trim(adjustl(PRTMOTION))=="YES")	then
	if(trim(adjustl(ANSMETHOD))=="YES")	then    
        MONOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.mst'
        open(6,file=MONOUTPUT)
        write(6,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(6,*)
		write(6,'(a)')"#              [ 运动响应传递函数计算结果文件]"
		write(6,*)
		write(6,'(a)')"#[航速数目,浪向角数目,波浪频率数目]"
		write(6,*)
		write(6,'(3i8)')NBU,NBHEAD,NBOME
		write(6,*)
		write(6,'(a)')"#[航速(kn)]"
		write(6,*)
		write(6,'(5f8.3)')(UKN(i),i=1,NBU)
		write(6,*)
        write(6,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(6,*)
        write(6,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        write(6,*)
        write(6,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(6,*)
        write(6,'(100f8.3)')(Ome0(i),i=1,NBOME)
!----------------------------------------------------
        ACEOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.acc'
        open(15,file=ACEOUTPUT)
        write(15,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(15,*)
        write(15,'(a)')"#              [ 浮体重心处加速度响应结果文件]"
        write(15,*)
        write(15,'(a)')"#[浪向角数目,波浪频率数目]"
        write(15,*)
        write(15,'(2i8)')NBHEAD,NBOME
        write(15,*)
        write(15,'(a)')"#[航速(kn)]"
        write(15,*)
        write(15,'(5f8.3)')(UKN(i),i=1,NBU)
        write(15,*)
        write(15,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(15,*)
        write(15,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        !---------------------2013.9.20lizhifu-------------------------------        
        GSTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.gst'
        open(116,file=GSTOUTPUT)
        write(116,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(116,*)
		write(116,'(a)')"#              [ 整船加速度响应传递函数计算结果文件]"
		write(116,*)
		write(116,'(a)')"#[航速数目,浪向角数目,波浪频率数目]"
		write(116,*)
		write(116,'(3i8)')NBU,NBHEAD,NBOME
		write(116,*)
		write(116,'(a)')"#[航速(kn)]"
		write(116,*)
		write(116,'(5f8.3)')(UKN(i),i=1,NBU)
		write(116,*)
        write(116,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(116,*)
        write(116,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        write(116,*)
        write(116,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(116,*)
        write(116,'(100f8.3)')(Ome0(i),i=1,NBOME)   
         !---------------------2014.4.20LiHui-------------------------------        
        MXTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.m4x'
        open(55,file=MXTOUTPUT)
        write(55,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(55,*)
		write(55,'(a)')"#              [ 横摇阻尼力矩响应]"
		write(55,*)
		write(55,'(a)')"#[航速数目,浪向角数目,波浪频率数目]"
		write(55,*)
		write(55,'(3i8)')NBU,NBHEAD,NBOME
		write(55,*)
		write(55,'(a)')"#[航速(kn)]"
		write(55,*)
		write(55,'(5f8.3)')(UKN(i),i=1,NBU)
		write(55,*)
        write(55,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(55,*)
        write(55,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        write(55,*)
        write(55,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(55,*)
        write(55,'(100f8.3)')(Ome0(i),i=1,NBOME)         
        !----------------------------------------------------------
	ELSE if(trim(adjustl(ANSMETHOD))=="NO")	then
	    MSTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.mst'
        open(131,file=MSTOUTPUT)
        write(131,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(131,*)
		write(131,'(a)')"#              [ 运动响应传递函数计算结果文件]"
		write(131,*)
		write(131,'(a)')"#[航速数目,浪向角数目,波浪频率数目]"
		write(131,*)
		write(131,'(3i8)')NBU,NBHEAD,NBOME
		write(131,*)
		write(131,'(a)')"#[航速(kn)]"
		write(131,*)
		write(131,'(5f8.3)')(UKN(i),i=1,NBU)
		write(131,*)
        write(131,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(131,*)
        write(131,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        write(131,*)
        write(131,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(131,*)
        write(131,'(100f8.3)')(Ome0(i),i=1,NBOME)
		
		ACEOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.acc'
        open(15,file=ACEOUTPUT)
        write(15,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(15,*)
        write(15,'(a)')"#              [ 浮体重心处加速度响应结果文件]"
        write(15,*)
        write(15,'(a)')"#[浪向角数目,波浪频率数目]"
        write(15,*)
        write(15,'(2i8)')NBHEAD,NBOME
        write(15,*)
        write(15,'(a)')"#[航速(kn)]"
        write(15,*)
        write(15,'(5f8.3)')(UKN(i),i=1,NBU)
        write(15,*)
        write(15,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(15,*)
        write(15,'(5f8.3)')(HEAD(i),i=1,NBHEAD)
        !------------------------------2013.9.20lizhifu---------------------------------------        
        GSTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.gst'
        open(231,file=GSTOUTPUT)
        write(231,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(231,*)
		write(231,'(a)')"#              [ 整船加速度响应传递函数计算结果文件]"
		write(231,*)
		write(231,'(a)')"#[航速数目,浪向角数目,波浪频率数目]"
		write(231,*)
		write(231,'(3i8)')NBU,NBHEAD,NBOME
		write(231,*)
		write(231,'(a)')"#[航速(kn)]"
		write(231,*)
		write(231,'(5f8.3)')(UKN(i),i=1,NBU)
		write(231,*)
        write(231,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(231,*)
        write(231,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        write(231,*)
        write(231,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(231,*)
        write(231,'(100f8.3)')(Ome0(i),i=1,NBOME)
                 !---------------------2014.4.20LiHui-------------------------------        
        MXTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.m4x'
        open(55,file=MXTOUTPUT)
        write(55,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(55,*)
		write(55,'(a)')"#              [ 横摇阻尼力矩响应]"
		write(55,*)
		write(55,'(a)')"#[航速数目,浪向角数目,波浪频率数目]"
		write(55,*)
		write(55,'(3i8)')NBU,NBHEAD,NBOME
		write(55,*)
		write(55,'(a)')"#[航速(kn)]"
		write(55,*)
		write(55,'(5f8.3)')(UKN(i),i=1,NBU)
		write(55,*)
        write(55,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
        write(55,*)
        write(55,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
        write(55,*)
        write(55,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(55,*)
        write(55,'(100f8.3)')(Ome0(i),i=1,NBOME)         
        
	endif
endif

!------------------------------------------------------------------------------------------------------
IF(trim(adjustl(PRTSECTLOAD))=="YES")THEN
    SSTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.sst'
    OPEN (132,FILE=SSTOUTPUT)
    write(132,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
    write(132,*)
    if(trim(adjustl(ANSMETHOD))=="YES")	then  
    write(132,'(a)')"#              [ 剖面载荷响应计算结果文件]"
    else
    write(132,'(a)')"#              [ 剖面载荷响应传递函数计算结果文件]"
    end if
    write(132,*)
    write(132,'(a)')"#[航速数目,浪向角数目,波浪频率数目,剖面数目]"
    write(132,*)
    if(trim(adjustl(MASSSOLVE))=="SECT") then 
	    write(132,'(4i8)')NBU,NBHEAD,NBOME,Nbsect
    else
	    write(132,'(4i8)')NBU,NBHEAD,NBOME,Nbsect
    end if
    write(132,*)
    write(132,'(a)')"#[航速(kn)]"
    write(132,*)
    write(132,'(5f8.3)')(UKN(i),i=1,NBU)
    write(132,*)
    write(132,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
    write(132,*)
    write(132,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
    write(132,*)
    write(132,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
    write(132,*)
    write(132,'(100f8.3)')(Ome0(i),i=1,NBOME)
    write(132,*)
    write(132,'(a)')"#[计算剖面号]"
    write(132,*)
    write(132,'(100i8)')(i,i=1,NBsect)
    write(132,*)
ENDIF
!------------------------------------------------------------------------------------------------------------
if(trim(adjustl(PRTva))=="YES")	then
    if(trim(adjustl(ANSMETHOD))=="YES")	then    
	    VACOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.ast'
        OPEN(19,FILE=VACOUTPUT)
	    write(19,"(a)")"#              [ COMPASS-WALCS-BASIC V1.0 ]"
	    write(19,*)
	    write(19,"(a)")"#              [ 计算点速度加速度响应计算结果文件]"
	    write(19,*)
	    write(19,'(a)')"#[航速数目,浪向角数目,波浪频率数目,计算点数目]"
	    write(19,*)
	    write(19,'(4i8)')NBU,NBHEAD,NBOME,NPVA
	    write(19,*)
	    write(19,'(a)')"#[航速(kn)]"
	    write(19,*)
	    write(19,'(5f8.3)')(UKN(i),i=1,NBU)
	    write(19,*)
	    write(19,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
	    write(19,*)
	    write(19,'(50f6.1)')(HEAD(i),i=1,NBHEAD)
	    write(19,*)
	    write(19,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
	    write(19,*)
	    write(19,'(100f8.3)')(Ome0(i),i=1,NBOME)
	    write(19,*)
	ELSE if(trim(adjustl(ANSMETHOD))=="NO") then
        ASTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.ast'
        OPEN(134,FILE=ASTOUTPUT)
		write(134,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
		write(134,*)
		write(134,'(a)')"#              [ 计算点速度加速度响应传递函数计算结果文件]"
		write(134,*)
		write(134,'(a)')"#[航速数目,浪向角数目,波浪频率数目,计算点数目]"
		write(134,*)
		write(134,'(4i8)')NBU,NBHEAD,NBOME,NPVA
		write(134,*)
		write(134,'(a)')"#[航速(kn)]"
		write(134,*)
		write(134,'(5f8.3)')(UKN(i),i=1,NBU)
		write(134,*)
		write(134,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
		write(134,*)
		write(134,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
		write(134,*)
	    write(134,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
	    write(134,*)
	    write(134,'(100f8.3)')(Ome0(i),i=1,NBOME)
	    write(134,*)
	endif
endif
if(trim(adjustl(PRTMOTION))=="YES")	then
        PREOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.pre'
        OPEN (3,FILE=PREOUTPUT)
        write(3,"(a)")"#              [ COMPASS-WALCS-BASIC V1.0 ]"
        write(3,*)""
        write(3,"(a)")"#              [ 压力响应计算结果输出文件 ]"
         write(3,*)""
	    write(3,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
	    write(3,*)
	    write(3,'(50f8.2)')(HEAD(i),i=1,NBHEAD)
	    write(3,*)
	    write(3,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
	    write(3,*)
	    write(3,'(100f8.3)')(Ome0(i),i=1,NBOME)
	    write(3,*)
end if

!------------------------------------------------------------------------------------------
if(trim(adjustl(PRTPPRESS))=="YES")then
    if(trim(adjustl(ANSMETHOD))=="YES")	then
	    PPROUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.pst' 
		OPEN(18,FILE=PPROUTPUT)  
	    write(18,*)"#              [ COMPASS-WALCS-BASIC V1.0 ]"
	    write(18,*)
	    write(18,*)"#              [ 计算点压力响应结果文件]"
	    write(18,*)
	    write(18,'(a)')"#[航速数目,浪向角数目,波浪频率数目,计算点数目]"
	    write(18,*)
	    write(18,'(4i8)')NBU,NBHEAD,NBOME,NPPOINT
	    write(18,*)
	    write(18,'(a)')"#[航速(kn)]"
	    write(18,*)
	    write(18,'(5f8.3)')(UKN(i),i=1,NBU)
	    write(18,*)
	    write(18,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
	    write(18,*)
	    write(18,'(50f8.2)')(HEAD(i),i=1,NBHEAD)
	    write(18,*)
	    write(18,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
	    write(18,*)
	    write(18,'(100f8.3)')(Ome0(i),i=1,NBOME)
	    write(18,*)
    ELSE if(trim(adjustl(ANSMETHOD))=="NO") then
	    PSTOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.pst'
		OPEN(133,FILE=PSTOUTPUT)
		write(133,'(a)')"#              [ COMPASS-WALCS-BASIC V1.0 ]"
		write(133,*)
		write(133,'(a)')"#              [ 计算点压力响应传递函数计算结果文件]"
		write(133,*)
		write(133,'(a)')"#[航速数目,浪向角数目,波浪频率数目,压力计算点数目]"
		write(133,*)
		write(133,'(4i8)')NBU,NBHEAD,NBOME,NPPoint
		write(133,*)
		write(133,'(a)')"#[航速(kn)]"
		write(133,*)
		write(133,'(5f8.3)')(UKN(i),i=1,NBU)
		write(133,*)
		write(133,'(a)')"#[浪向角(deg)HEAD(i),i=1,NBHEAD]"
		write(133,*)
		write(133,'(50f8.3)')(HEAD(i),i=1,NBHEAD)
		write(133,*)
		write(133,'(a)')"#[波浪频率(rad/s)Ome0(i),i=1,NBOME]"
        write(133,*)
        write(133,'(100f8.3)')(Ome0(i),i=1,NBOME)
        write(133,*)
	endif
endif
!------------------------------------------------------------------------------------------------------
IF(trim(adjustl(BILGEFORCE))=="YES")THEN
    BVFOUTPUT=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.bvf'
    OPEN (119,FILE=BVFOUTPUT)
    write(119,*)"#              [ COMPASS-WALCS-BASIC V1.0 ]"
    write(119,*) 
    write(119,'(a)')"#[VISCOUS ROLL DAMPING FOR FORCES ON BILGE KEELS]"
    write(119,*)
    write(119,'(a)')"#[NBU,NBHEAD,NBOME,Npvf]"
    write(119,*)
    write(119,'(4i8)')NBU,NBHEAD,NBOME,Npvf
    write(119,*)
    write(119,'(a)')"#[UKN(i),i=1,NBU]"
    write(119,*)
    write(119,'(5f8.3)')(UKN(i),i=1,NBU)
    write(119,*)
    write(119,'(a)')"#[HEAD(i),i=1,NBHEAD]"
    write(119,*)
    write(119,'(5f8.3)')(HEAD(i),i=1,NBHEAD)
    write(119,*)
    write(119,'(a)')"#[Ome0(i),i=1,NBOME]"
    write(119,*)
    write(119,'(10f8.3)')(Ome0(i),i=1,NBOME)
    write(119,*)
ENDIF
!------------------------------------------------------------------------------------------------------------
END SUBROUTINE OUTPUTFILE



!------------------------------------------------------------------------------------------
!输出运动结果
!------------------------------------------------------------------------------------------
subroutine m_putouts(IB,motion,M4x)
implicit none

INTEGER,INTENT(IN)::IB
real(8),intent(in),dimension(1:2,1:6,0:NBOME)::motion
real(8),intent(in),dimension(1:2,0:NBOME)::M4x
real(8),dimension(1:2,1:6)::temp
real(8),dimension(1:2,1:2)    ::  Acc
real(8),dimension(0:NBOME) :: phai,theta
integer(4)::i,j,k

if(trim(adjustl(ANSMETHOD))=="NO")then
!-------------运动-------------------

    WRITE(131,*)
    WRITE(131,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(131,'(a,4x,a,5x,a,6x,a,7x,a,6x,a,6x,a,6x,a,7x,a,6x,a,6x,a,6x,a,8x,a,6x,a,1x)')"#   Ome0","Ome","Surge","PHA","Sway","PHA","Heave","PHA","Roll","PHA","Pitch","PHA","Yaw","PHA"
	do i=1,NBOME
		do j=1,6
			call comp_to_ampha(motion(1,J,I),motion(2,J,I),temp(1,J),temp(2,J))
		enddo
		temp(1,:)=temp(1,:)/(amp(i))
		write(131,'(2f8.3,6(e12.4,f8.3))')ome0(i),ome1(i),((temp(K,J),k=1,2),j=1,6)
	enddo
!-----------加速度----------------------------------------------------------    	
    WRITE(15,*)
	WRITE(15,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(15,'(a,3x,a,1x,12(4x,a,3x))')"#  Ome0","Ome","SurAC","SurAS","SwaAC","SwaAS","HeaAC","HeaAS","RolAC","RolAS","PitAC","PitAS","YawAC","YawAS"
		
	do i=1,NBOME
	
        do j=1,6
	        temp(1,J)=-ome(i)*ome(i)*motion(1,J,I)
	        temp(2,J)=-ome(i)*ome(i)*motion(2,J,I)
        enddo
        
        temp(1,1)=temp(1,1)-g0*motion(1,5,i)
	    temp(2,1)=temp(2,1)-g0*motion(2,5,i)
	    
	    temp(1,2)=temp(1,2)+g0*motion(1,4,i) 
	    temp(2,2)=temp(2,2)+g0*motion(2,4,i) 
	    
		write(15,'(2f7.3,12e12.4)')ome0(i),ome1(i),((temp(J,K),j=1,2),k=1,6)
		
	enddo	
	
	
	!---------GST-------------2013.9.20lizhifu---------------------------------
	WRITE(231,*)
    WRITE(231,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(231,'(a,4x,a,5x,a,4x,a,6x,a,6x,a,4x,a,5x,a,6x,a,6x,a,4x,a,4x,a,7x,a,6x,a,1x)')"#   Ome0","Ome","Surge_G","PHA","Sway_G","PHA","Heave_G","PHA","Roll_G","PHA","Pitch_G","PHA","Yaw_G","PHA"
	do i=1,NBOME
	    do j=1,6
			call comp_to_ampha(-ome(i)*ome(i)*motion(1,J,I),-ome(i)*ome(i)*motion(2,J,I),temp(1,J),temp(2,J))
		enddo
	!--------------------------------2013.11.17shiyuyun-------------------	
	
!	        temp(2,:)  =	temp(2,:)+180		
!			
!			if(temp(2,J)>360)then
!			
!			    temp(2,J)=temp(2,J)-360
!			    
!			else if(temp(2,J)<-360)then
!			
!			    temp(2,J)=temp(2,J)+360 
!			    
!			end if 			
!		enddo
!		temp(1,:)=temp(1,:)/(amp(i))
		
!----------------重心处加速度：考虑重力的贡献-By LiHui 2014.04.20------------------------------	

!		temp(1,:)=ome(i)*ome(i)*temp(1,:)
		
		Acc(1,1) =  -ome(i)*ome(i)*motion(1,1,i)-g0*motion(1,5,i)
		
		Acc(2,1) =  -ome(i)*ome(i)*motion(2,1,i)-g0*motion(2,5,i)
		
		
		Acc(1:2,2) =  -ome(i)*ome(i)*motion(1:2,2,i)+g0*motion(1:2,4,i) 
		
		do  j=  1,   2 
		 
		    call comp_to_ampha(Acc(1,j),Acc(2,j),temp(1,j),temp(2,j))
		    
		enddo    
		
		write(231,'(2f8.3,6(e12.4,f8.3))')ome0(i),ome1(i),((temp(K,J),k=1,2),j=1,6)
	enddo
!--------------------横摇阻尼力矩------------------------------------------------------
    WRITE(55,*)
	WRITE(55,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(55,'(a,3x,a,1x,12(4x,a,3x))')"#  Ome0","Ome","M4xC","M4xS"
	
    do i=1,NBOME
	    write(55,'(2f7.3,2e12.4)')ome0(i),ome1(i),M4x(1,i),M4x(2,i)
	end do
!---------------------谱分析Ok----------------------------------------------------------
!---------------------------------------------------------------------------------------
ELSE if(trim(adjustl(ANSMETHOD))=="YES")then
    WRITE(6,*)
    WRITE(6,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(6,'(a,4x,a,5x,a,6x,a,7x,a,6x,a,6x,a,6x,a,7x,a,6x,a,6x,a,6x,a,8x,a,6x,a,1x)')"#   Ome0","Ome","Surge","PHA","Sway","PHA","Heave","PHA","Roll","PHA","Pitch","PHA","Yaw","PHA"

    do i=1,NBOME
		do j=1,6
		
			call comp_to_ampha(motion(1,J,I),motion(2,J,I),temp(1,J),temp(2,J))

		enddo	
				
		temp(2,:)=temp(2,:)+Phase(i)    !考虑入射波相位
				
		temp(1,:)=temp(1,:)*cos(temp(2,:)*PI/180.0)
		
		theta(i)   =   temp(1,4)   !横摇角		
		phai(i)    =   temp(1,5)   !纵摇角

		write(6,'(2f8.3,6(e12.4,f8.3))')ome0(i),ome1(i),((temp(K,J),k=1,2),j=1,6)
		
	enddo


!----------------------------------------------------------------------------    
    WRITE(15,*)
	WRITE(15,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(15,*)"#Given incident wave phase values(rad for angle):"
	write(15,'(3x,a,3x,a,2x,3x,a,2x,3x,a,3x,3x,a,2x,3x,a,3x,3x,a,2x,4x,a,3x)')"#   Ome0","Ome","Surge_a","Sway_a","Heave_a","Roll_a","Pitch_a","Yaw_a"
!-----------------------------------------------------------------------------	
	WRITE(116,*)
    WRITE(116,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(116,'(a,4x,a,5x,a,4x,a,6x,a,6x,a,4x,a,5x,a,6x,a,6x,a,4x,a,4x,a,7x,a,6x,a,1x)')"#   Ome0","Ome","Surge_G","PHA","Sway_G","PHA","Heave_G","PHA","Roll_G","PHA","Pitch_G","PHA","Yaw_G","PHA"
	
    do i=1,NBOME
    
		do j=1,6
		
			call comp_to_ampha(-ome(i)*ome(i)*motion(1,J,I),-ome(i)*ome(i)*motion(2,J,I),temp(1,J),temp(2,J))

		enddo	
				
		temp(2,:)=temp(2,:)+Phase(i)    !考虑入射波相位
		
		temp(1,:)=temp(1,:)*cos(temp(2,:)*PI/180.0)
		
		temp(1,1)=temp(1,1)-g0*sin(phai(i))
		
		temp(1,2)=temp(1,2)+g0*sin(theta(i))

	    write(15,'(2f7.3,6e12.4)')ome0(i),ome1(i),(temp(1,J),j=1,6)	
	    
		write(116,'(2f8.3,6(e12.4,f8.3))')ome0(i),ome1(i),((temp(K,J),k=1,2),j=1,6)   
		 
	enddo

  
  !------------------------------------------------------------------------------------------
    WRITE(55,*)
	WRITE(55,'(A,I4,A,F5.1)')"#IB=",IB,"     HEAD= ",HEAD(IB)*180.0/PI
	write(55,'(a,3x,a,1x,12(4x,a,3x))')"#  Ome0","Ome","M4xC","M4xS"
	
	do i=1,NBOME
	
	    call comp_to_ampha(M4x(1,i),M4x(2,i),temp(1,i),temp(2,i))
	    
		temp(2,i)=temp(2,i)+Phase(i)    !考虑入射波相位
			
		if(temp(2,i)>360)then
			
		    temp(2,i)=temp(2,i)-360
			    
		else if(temp(2,i)<-360)then
			
		    temp(2,i)=temp(2,i)+360 
			    
		end if	
		
		temp(1,i)=temp(1,i)*cos(temp(2,i)*PI/180.0)	
			
	    write(55,'(2f7.3,2e12.4)')ome0(i),ome1(i),temp(1,i),temp(2,i)
	end do


!----------------------------------------------------------------------------------
!原程序输出.ace文件时此处有误
endif
end subroutine m_putouts
!------------------------------------------------------


!----------------------------------------------------------------------
END MODULE OUTPUTFILE_MOD