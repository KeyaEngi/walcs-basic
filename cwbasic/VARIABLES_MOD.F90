!*************************************************************************************
!*************************************************************************************
!路径变量  ACCESS.DAT  config.hyp
MODULE INOUTACCESS_MOD  
IMPLICIT NONE
SAVE
character(100)	::	PROJNAME
character(200)	::	InAccess,OutAccess

END MODULE INOUTACCESS_MOD
!*******************************************************************************
!水动力网格信息变量
MODULE MESH_MOD   
IMPLICIT NONE
SAVE
CHARACTER(LEN=15)::MESHTYPE    !BDF  MES
CHARACTER(LEN=15)::SYMTYPE    !是否考虑对称性 
INTEGER(4)::TarNum,Space  !是对称性的修正系数  edited by 王川2012.09.16
INTEGER(4)::Factor  !是对称性的修正系数  edited by 王川2012.09.16
INTEGER(4)::NPanel !加刚盖的总数 加内部自由面网格 总网格数
INTEGER(4)::NB !湿表面网格
INTEGER(4)::ArcNum !每剖面网格数
! INTEGER(4),ALLOCATABLE,DIMENSION(:)::IEC !面元类型
REAL(8),ALLOCATABLE,dimension(:,:,:)::xn !网格顶点坐标
! REAL(8),ALLOCATABLE,dimension(:,:,:)::TEMPxn !网格顶点坐标
REAL(8),ALLOCATABLE,dimension(:,:)::xav   !网格中心点坐标
REAL(8),ALLOCATABLE,dimension(:,:,:)::XQ   !网格顶点在局部坐标
REAL(8),ALLOCATABLE,dimension(:,:,:)::e   
REAL(8),ALLOCATABLE,dimension(:)::EA 
real(8),	dimension(1:3,1:3)	::	Ty				!纵倾坐标变换矩阵
real(8)::vol,XB,YB,ZB,deltax  !网格体积，浮心坐标
real(8)::vol_x,vol_y,vol_z
real(8)::SL,SB,ST,Ta,TrimAng
REAL(8),ALLOCATABLE,dimension(:)::sectEA
REAL(8)::WETEA
real(8),dimension(1:6,1:6)::HRM  !hydrostatic restoring matrix静水恢复力矩阵
real(8),ALLOCATABLE,dimension(:,:,:)::SECTHRM !截面的静水恢复力矩阵

END MODULE MESH_MOD
!***************************************************************************************
!环境和尺度信息变量
MODULE ENVIRONMENT_MOD   
IMPLICIT NONE
SAVE

REAL(8)::ROU, H 
REAL(8),parameter::	PI=3.14159265358979
REAL(8),parameter::   g0=9.81

REAL(8)::TTMASS
REAL(8)::ZG ,XG ,YG         !重心在网格输入坐标系（z=0.在静水面处）下的z坐标,重心在模型坐标系中的x坐标
INTEGER(4)::NITEM , NBSECT	
REAL(8),ALLOCATABLE,dimension(:,:)::SECTCOG,SECTPSN  !SECTCOG:截断剖面到重心系的坐标；SECTPSN各剖面在初始网格坐标系下的位置
REAL(8),dimension(1:6,1:6)::MASSMATRIX
REAL(8),ALLOCATABLE,dimension(:,:,:)::SECTMASSMATRIX
real(8),allocatable::XM(:),YM(:),ZM(:),MM(:),IXX(:),IYY(:),IZZ(:),IXY(:),IXZ(:),IYZ(:),ZSC(:)
integer(4),allocatable::SECTTYPE(:)   !SECTTYPE=1横剖面 SECTTYPE=2纵剖面

real(8),allocatable,dimension(:,:,:)::SML   !单根锚链刚度阵
INTEGER::MLN1    !锚链根数
real(8),allocatable,dimension(:)::XML,YML,ZML  !XL-转塔与系泊油船艉间的距离,ZL-着链点与系泊油船基线间的距离

INTEGER::NBU,NBHEAD,NBOME   !number of 航速，浪向，频率
REAL(8),ALLOCATABLE,dimension(:)::UKN   
REAL(8),ALLOCATABLE,dimension(:)::HEAD 
REAL(8),ALLOCATABLE,dimension(:)::OME0
REAL(8),ALLOCATABLE,dimension(:)::OME 
REAL(8),ALLOCATABLE,dimension(:)::OME1
REAL(8),ALLOCATABLE,dimension(:)::AMP
REAL(8),ALLOCATABLE,dimension(:)::PHASE

END MODULE  ENVIRONMENT_MOD
!*******************************************************************************************
!有限元网格信息参数
MODULE PVPOINT_MOD 
IMPLICIT NONE
SAVE
!-----------------------------------------------------
CHARACTER(LEN=15)::PRESSNAME  !YES or NO
INTEGER::NLC
!----------------------------------------------------
INTEGER(4)::NFEPOINT
REAL(8),ALLOCATABLE,DIMENSION(:,:)::FEMPOINT,NORMAL
integer(4),allocatable,dimension(:)::WetEleID
character(len=100)::Endchar
character(len=20),allocatable,dimension(:,:,:)::LCName
!-----------------------------------------------------
INTEGER::NPVA
REAL(8),ALLOCATABLE,DIMENSION(:,:)::VAPOINT
!-----------------------------------------------------
REAL(8),ALLOCATABLE,DIMENSION(:,:)::PRESSPOINT
INTEGER(4)::NPPOINT
!-----------------------------------------------------
INTEGER::Npvf   !舭龙骨附近等效节点个数
real(8),allocatable,dimension(:,:)::vForceNode    !等效节点坐标
!REAL(8)::PGL   !计算中间变量
integer::FEM_N  !有限元模型的法向标志==0向内，==1向外
END MODULE PVPOINT_MOD
!*****************************************************************************************************
!剖面载荷
MODULE WAVEFORCE_MOD  
IMPLICIT NONE
SAVE
REAL(8),ALLOCATABLE,dimension(:,:,:)::SECTLOAD

END  MODULE WAVEFORCE_MOD
!******************************************************************************************************
!输入计算参数选择变量
MODULE CAL_MOD
USE INOUTACCESS_MOD   
IMPLICIT NONE
SAVE

CHARACTER(LEN=15)::WATERDEPTH  !FINITE  INFINITE
!INTEGER::infiFreq    !是否计算频率为无穷时的附加质量(0否 1是)
!INTEGER::zeroFreq    !是否计算频率为0时的附加质量(0否 1是)
CHARACTER(LEN=15)::INPUTSTATIC !YES or NO 是否输入船舶静力特性Aw,Sy,hx,hy
!CHARACTER(LEN=15)::INPUTINERTIA !YES or NO 是否输入船舶质量特性I11，I22，I33，I13
integer::ROllDAMPING 
!横摇修正方法	
!	0---横摇不修正，
!	1---ram,rmb法,
!	2---ram,rmb,ek0法
!	3---米勒法(需输入舭龙骨数据)
!	4---驳船经验公式
!   5---临界阻尼系数  
!CHARACTER(LEN=15)::MATRIXSOLVE !GAUSS GMRES
CHARACTER(LEN=15):: MASSSOLVE !WHOLE  SECT
CHARACTER(LEN=15):: SECTLOADSOLVE !YES  NO

! CHARACTER(LEN=15)::WAVEFI        !QUAS  PRECISE(波浪主干扰力计算方法：中心点乘面积 精确积分法)
INTEGER::POINTPRESS    !对角线坐标变换插值0 贴网格法1 直接计算2
CHARACTER(LEN=15)::MoorLine      !YES  NO
CHARACTER(LEN=15)::ModiDeflect !YES  NO
CHARACTER(LEN=15)::HULLTYPE   !单体or多体（SINGLE or MULTI）
CHARACTER(LEN=15)::OMETYPE    !波浪频率输入方式FREQUENCY or PERIOD
CHARACTER(LEN=15)::FLOATATION !用户自定义or质量计算(USER or FMAS or NO)

!---------------------------------------------
ENDMODULE CAL_MOD
!*******************************************************************************************************
!输出参数选择变量
MODULE PRINT_MOD              
IMPLICIT  NONE
SAVE
CHARACTER(LEN=15)::PRTSTATIC !YES NO
CHARACTER(LEN=15)::PRTMOTION !YES NO
CHARACTER(LEN=15)::PRTPPRESS   !YES NO
CHARACTER(LEN=15)::PRTSECTLOAD !YES NO
CHARACTER(LEN=15)::PRTFEPRESS  !NO PATRAN ANSYS
CHARACTER(LEN=15)::PRTva !YES NO
CHARACTER(LEN=15)::ansMethod   !SPECTRUM DESIGNWAVE YES OR NO
CHARACTER(LEN=15)::BILGEFORCE !YES  NO
END  MODULE PRINT_MOD
!*********************************************************************************************************
!输出文件
MODULE OUTPUT_MOD  
IMPLICIT NONE
SAVE

CHARACTER(LEN=300)::SSTOUTPUT    !剖面载荷输出文件
CHARACTER(LEN=300)::OUTPUT3HC    !水动力系数及波浪干扰力输出文件
CHARACTER(LEN=300)::SLTOUTPUT    
CHARACTER(LEN=300)::MONOUTPUT
CHARACTER(LEN=300)::MSTOUTPUT
CHARACTER(LEN=300)::ACEOUTPUT
CHARACTER(LEN=300)::VACOUTPUT
CHARACTER(LEN=300)::ASTOUTPUT
CHARACTER(LEN=300)::PPROUTPUT
CHARACTER(LEN=300)::PSTOUTPUT
CHARACTER(LEN=300)::BVFOUTPUT
CHARACTER(LEN=300)::PREOUTPUT
CHARACTER(LEN=300)::GSTOUTPUT    !重心处加速度输出文件 byshiyuyun
CHARACTER(LEN=300)::MXTOUTPUT    !横摇阻尼力矩输出文件 LiHui

END  MODULE OUTPUT_MOD
!*********************************************************************************************************