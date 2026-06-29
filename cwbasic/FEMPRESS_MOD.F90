MODULE FEMPRESS_MOD

USE INOUTACCESS_MOD
USE ENVIRONMENT_MOD,ONLY:NBU,NBHEAD,NBOME,PI
USE MESH_MOD,ONLY:Ta,ST,SB,DELTAX,ArcNum,NB,EA,E,XAV,XQ,TY,Factor,TarNum,Space
USE CAL_MOD,ONLY:POINTPRESS
USE PVPOINT_MOD,ONLY:NFEPOINT,FEMPOINT,NORMAL,WetEleID,PRESSNAME,NLC,Endchar,LCname,FEM_N
USE PRINT_MOD,ONLY:PRTFEpRESS,ansMethod

IMPLICIT NONE
PRIVATE
PUBLIC::FEPPOINT,FEPRESS


CONTAINS
!------------------------------------------------------------------------------------
SUBROUTINE FEPPOINT(RATIO)

real(8),allocatable,dimension(:,:)::PressPoint0,normal0
integer(4),	allocatable,dimension(:)::WetEleID0
real(8),dimension(1:3)		::	tempcoor
REAL,INTENT(OUT)::ratio
!integer   ::  FEM_N
CHARACTER(LEN=300)::FEMINPUT,BDFINPUT,DATINPUT
! character(len=20),INTENT(OUT),dimension(:,:,:)	::LCName
character(len=100)::BDFNAME,ANSNAME
! character(len=100)::Endchar
INTEGER::I, J, K
INTEGER::TEMP
!--------------------------------------------------------------------------------------
FEMINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.FEM'

if(TRIM(ADJUSTL(PRTFEPRESS))=="YES")PRTFEPRESS="PATRAN"

if(TRIM(ADJUSTL(PRTFEPRESS))=="PATRAN")	then
	open(31,file=FEMINPUT,status='old')
	call IUTMP(31)
	read(31,*)	bdfname
	read(31,*)(tempcoor(i),i=1,3),ratio,FEM_N  !压力计算点所在坐标系在初始坐标系下坐标
	tempcoor(1)=deltax-tempcoor(1)/ratio          !随船平动坐标系在压力计算点坐标系下的坐标
	tempcoor(2)=-tempcoor(2)/ratio
	tempcoor(3)=ST-tempcoor(3)/ratio
!	tempcoor(1)=deltax-tempcoor(1)/ratio         !随船平动坐标系在压力计算点坐标系下的坐标
!	tempcoor(2)=-tempcoor(2)/ratio
!	tempcoor(3)=Ta-tempcoor(3)/ratio
	read(31,*)	PRESSNAME
	if(TRIM(ADJUSTL(PRESSNAME))=="YES")	then
		if(TRIM(ADJUSTL(ansMethod))=="NO")	then
			NLC=2*NBOME
			allocate(Lcname(1:NLC,1:NBHEAD,1:NBU))
			do	I=1,	NBU
				do	J=1,	NBHEAD
					do	K=1,	NLC
						read(31,*)temp,Lcname(K,J,I)
					end	do
				end	do
			end	do
		elseif(TRIM(ADJUSTL(ansMethod))=="YES")	then
			NLC=NBOME
			allocate(Lcname(1:NLC,1:NBHEAD,1:NBU))
			do	I=1,	NBU
				do	J=1,	NBHEAD
					do	K=1,	NLC
						read(31,*)temp,Lcname(K,J,I)
					end	do
				end	do
			end	do
		endif
	endif
	BDFINPUT=trim(adjustl(InAccess))//'\'//trim(adjustl(bdfname))
	allocate(PressPoint0(1:3,1:1000000),normal0(1:3,1:1000000),WetEleID0(1:1000000))

	call	GetPpoint(BDFINPUT,PressPoint0,WetEleID0,tempcoor,ratio,normal0)   !获得patran输出的bdf文件的有限元网格信息
	
	
    ALLOCATE(FEMPOINT(1:3,1:NFEPOINT),NORMAL(1:3,1:NFEPOINT),WetEleID(1:NFEPOINT))
	DO I=1,NFEPOINT
	    FEMPOINT(:,I)=PressPoint0(:,I)
        normal(:,I)=normal0(:,I)
		WetEleID(I)=WetEleID0(I)
	ENDDO
	DEALLOCATE(PressPoint0,normal0,WetEleID0)
	!----------纵倾坐标变换----------
!	do i=1,NFEPOINT
!		FEMPOINT(1:3,i)=Matmul(Ty,FEMPOINT(1:3,i))
!	enddo
	close(31)
elseif(TRIM(ADJUSTL(PRTFEPRESS))=="ANSYS")	then

	open(31,file=FEMINPUT,status='old')
	call IUTMP(31)
	read(31,*)	Ansname
	read(31,*)(tempcoor(i),i=1,3),ratio           !压力计算点所在坐标系在初始坐标系下坐标
	tempcoor(1)=deltax-tempcoor(1)/ratio          !随船平动坐标系在压力计算点坐标系下的坐标
	tempcoor(2)=-tempcoor(2)/ratio
	tempcoor(3)=ST-tempcoor(3)/ratio
!    tempcoor(1)=deltax-tempcoor(1)/ratio          !随船平动坐标系在压力计算点坐标系下的坐标
!	tempcoor(2)=-tempcoor(2)/ratio
!	tempcoor(3)=Ta-tempcoor(3)/ratio
	
	read(31,*)	PRESSNAME
	
	if(TRIM(ADJUSTL(PRESSNAME))=="YES")	then
	
		if(TRIM(ADJUSTL(ansMethod))=="NO")	then
			NLC=2*NBOME
			allocate(Lcname(1:NLC,1:NBHEAD,1:NBU))
			do	I=1,	NBU
				do	J=1,	NBHEAD
					do	K=1,	NLC
						read(31,*)temp,Lcname(K,J,I)
					end	do
				end	do
			end	do
		elseif(TRIM(ADJUSTL(ansMethod))=="YES")	then
			NLC=NBOME
			allocate(Lcname(1:NLC,1:NBHEAD,1:NBU))
			do	I=1,	NBU
				do	J=1,	NBHEAD
					do	K=1,	NLC
						read(31,*)temp,Lcname(K,J,I)
					end	do
				end	do
			end	do
		endif
	endif
	
	DATINPUT=trim(adjustl(InAccess))//'\'//trim(adjustl(Ansname))
	
	allocate(PressPoint0(1:3,1:1000000),normal0(1:3,1:1000000),WetEleID0(1:1000000))

	call	GetAnsPoint(DATINPUT,PressPoint0,WetEleID0,tempcoor,ratio,normal0)   !获得ANSYS输出的dat文件的有限元网格信息
	
    ALLOCATE(FEMPOINT(1:3,1:NFEPOINT),NORMAL(1:3,1:NFEPOINT),WetEleID(1:NFEPOINT))
    
	DO I=1,NFEPOINT
	
	    FEMPOINT(:,I)=PressPoint0(:,I)
        normal(:,I)=normal0(:,I)
		WetEleID(I)=WetEleID0(I)
		
	ENDDO
	
	DEALLOCATE(PressPoint0,normal0,WetEleID0)
	!----------纵倾坐标变换----------
!	do i=1,NFEPOINT
!	
!		FEMPOINT(1:3,i)=Matmul(Ty,FEMPOINT(1:3,i))
!		
!	enddo
	
	close(31)
	
endif
END SUBROUTINE FEPPOINT
!-------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE FEPRESS(IU,IB,IP,BETA,OMEI0,wavephase,PS,PW,PD,PR,PALL,RATIO)
IMPLICIT NONE

REAL(8),INTENT(IN)::BETA,OMEI0,wavephase
INTEGER,intent(in)::IU,IB,IP
real(8),intent(in),dimension(1:2,1:Factor,1:NB)::ps,pw,pd,pr,pall
REAL,INTENT(IN)::RATIO
! character(len=20),intent(in),dimension(1:NLC)::LCName
real(8),allocatable,dimension(:,:)::pps,ppw,ppd,ppr,ppall	!船舶计算点表面压力
real(8),allocatable,dimension(:,:)::zps,zpw,zpd,zpr,zpall
integer::i,j
allocate(pps(1:2,1:NFEPoint),ppw(1:2,1:NFEPoint),ppd(1:2,1:NFEPoint),ppr(1:2,1:NFEPoint),ppall(1:2,1:NFEPoint))

if(POINTPRESS==0)	then
     allocate(zps(1:2,1:Factor*NB),zpw(1:2,1:Factor*NB),zpd(1:2,1:Factor*NB),zpr(1:2,1:Factor*NB),zpall(1:2,1:Factor*NB))
    !x坐标赋值  满足原来插值需求 edited by 王川
    if(Factor==1)then
        zps(:,:)=ps(:,1,:)
        zpw(:,:)=pw(:,1,:)
        zpd(:,:)=pd(:,1,:)
        zpr(:,:)=pr(:,1,:)
        zpall(:,:)=pall(:,1,:)
    else if(Factor==2)then
        do i=1,2*NB/ArcNum
            do j=1,ArcNum/2
                !zps赋值
                zps(:,j+(i-1)*ArcNum)=ps(:,1,j+(i-1)*ArcNum/2)
                zps(:,j+(i-1)*ArcNum+ArcNum/2)=ps(:,2,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                
                !zpw赋值
                zpw(:,j+(i-1)*ArcNum)=pw(:,1,j+(i-1)*ArcNum/2)
                zpw(:,j+(i-1)*ArcNum+ArcNum/2)=pw(:,2,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                
                !zpd赋值
                zpd(:,j+(i-1)*ArcNum)=pd(:,1,j+(i-1)*ArcNum/2)
                zpd(:,j+(i-1)*ArcNum+ArcNum/2)=pd(:,2,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                
                !zpr赋值
                zpr(:,j+(i-1)*ArcNum)=pr(:,1,j+(i-1)*ArcNum/2)
                zpr(:,j+(i-1)*ArcNum+ArcNum/2)=pr(:,2,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                    
                !zpall赋值
                zpall(:,j+(i-1)*ArcNum)=pall(:,1,j+(i-1)*ArcNum/2)
                zpall(:,j+(i-1)*ArcNum+ArcNum/2)=pall(:,2,(i-1)*ArcNum/2+ArcNum/2+1-j) 
            end do
        end do   
    else if(Factor==4)then
        do i=1,2*NB/ArcNum
            do j=1,ArcNum/2
                !zps赋值
                zps(:,j+(i-1)*ArcNum)=ps(:,2,j+(i-1)*ArcNum/2)
                zps(:,j+(i-1)*ArcNum+ArcNum/2)=ps(:,3,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                zps(:,2*NB+j+(i-1)*ArcNum)=ps(:,1,j+NB-i*ArcNum/2)
                zps(:,2*NB+(i-1)*ArcNum+ArcNum/2+j)=ps(:,4,NB-(i-1)*ArcNum/2+1-j)
                
                !zpw赋值
                zpw(:,j+(i-1)*ArcNum)=pw(:,2,j+(i-1)*ArcNum/2)
                zpw(:,j+(i-1)*ArcNum+ArcNum/2)=pw(:,3,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                zpw(:,2*NB+j+(i-1)*ArcNum)=pw(:,1,j+NB-i*ArcNum/2)
                zpw(:,2*NB+(i-1)*ArcNum+ArcNum/2+j)=pw(:,4,NB-(i-1)*ArcNum/2+1-j)
                
                !zpd赋值
                zpd(:,j+(i-1)*ArcNum)=pd(:,2,j+(i-1)*ArcNum/2)
                zpd(:,j+(i-1)*ArcNum+ArcNum/2)=pd(:,3,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                zpd(:,2*NB+j+(i-1)*ArcNum)=pd(:,1,j+NB-i*ArcNum/2)
                zpd(:,2*NB+(i-1)*ArcNum+ArcNum/2+j)=pd(:,4,NB-(i-1)*ArcNum/2+1-j)
                
                !zpr赋值
                zpr(:,j+(i-1)*ArcNum)=pr(:,2,j+(i-1)*ArcNum/2)
                zpr(:,j+(i-1)*ArcNum+ArcNum/2)=pr(:,3,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                zpr(:,2*NB+j+(i-1)*ArcNum)=pr(:,1,j+NB-i*ArcNum/2)
                zpr(:,2*NB+(i-1)*ArcNum+ArcNum/2+j)=pr(:,4,NB-(i-1)*ArcNum/2+1-j)
                        
                !zpall赋值
                zpall(:,j+(i-1)*ArcNum)=pall(:,2,j+(i-1)*ArcNum/2)
                zpall(:,j+(i-1)*ArcNum+ArcNum/2)=pall(:,3,(i-1)*ArcNum/2+ArcNum/2+1-j) 
                zpall(:,2*NB+j+(i-1)*ArcNum)=pall(:,1,j+NB-i*ArcNum/2)
                zpall(:,2*NB+(i-1)*ArcNum+ArcNum/2+j)=pall(:,4,NB-(i-1)*ArcNum/2+1-j)
            end do
        end do
    end if    
    call tranpress3 (xav(:,1:NB),zps,zpw,zpd,zpr,zpall,pps,ppw,ppd,ppr,ppall)  !分别沿X,Z方向插值传递压力(hydro to FEM)
    
elseIF(POINTPRESS==1) THEN	
		
 	call Precal	(ps,pw,pd,pr,pall,pps,ppw,ppd,ppr,ppall)    !贴体网格法传递压力(hydro to FEM)
 	allocate(zps(1:2,1:1),zpw(1:2,1:1),zpd(1:2,1:1),zpr(1:2,1:1),zpall(1:2,1:1))
 	
endif

if(IP>0)   then

    call ep_putouts	(IU,IB,IP,omeI0,ppall,ratio,beta,wavephase)  !输出有限元网格压力

endif

DEALLOCATE(pps,ppw,ppd,ppr,ppall,zps,zpw,zpd,zpr,zpall)

END SUBROUTINE FEPRESS
!--------------------------------------------------------------------------------------------------




!----------------------------------------------------------------------------------
!输出用于有限元加载的压力结果
!----------------------------------------------------------------------------------
subroutine ep_putouts(IU,IB,IP,omeI0,ppall,ratio,BETA,wavephase)
implicit none

integer,intent(in)::IU,IB,IP
REAL(8),INTENT(IN)::BETA,OMEI0,wavephase
real,INTENT(IN)::ratio
!real(8),intent(in),dimension(1:3,1:NFEPoint)::FEMPOINT
real(8),intent(in),dimension(1:2,1:NFEPoint)::ppall
! character(len=20),INTENT(IN),dimension(1:NLC)	::	LCName
! character(len=100),INTENT(IN)::	Endchar

real(8),dimension(1:2,1:NFEPoint)::ppresure
real(8),dimension(1:2,1:NFEpoint)::temp
integer(4)::LoadSetC,LoadSets
integer(4)::i
! real(8)::PI
character::B*2,O*2
character(len=10),dimension(1:NFEpoint)	::	CharP

! PI=4.0*atan(1.0)
call QINTTOSTR(IB,B,2)
! do j=1,NBOME
	if(TRIM(ADJUSTL(ansMethod))=="NO")	then
		LoadSetC=(((IU-1)*NBHEAD+(IB-1))*NBOME+IP)*2-1
		LoadSetS=(((IU-1)*NBHEAD+(IB-1))*NBOME+IP)*2
		if(TRIM(ADJUSTL(PRTFEPRESS))=="PATRAN")	then
			if(TRIM(ADJUSTL(PRESSNAME))=="NO")	then
				call QINTTOSTR(IP,O,2)			
				open(333,file=trim(adjustl(OutAccess))//'\'//'EPBeta'//B//'Ome'//O//'C.bdf')
				open(334,file=trim(adjustl(OutAccess))//'\'//'EPBeta'//B//'Ome'//O//'S.bdf')
			else IF(TRIM(ADJUSTL(PRESSNAME))=="YES") THEN
				open(333,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(LCName(2*IP-1,IB,IU)))//'.bdf')
				open(334,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(LCName(2*IP,IB,IU)))//'.bdf')
			endif
			write(333,'(a,I3,2x,a)')"$ LoadCase : ",LoadSetC,"cos"
			write(334,'(a,I3,2x,a)')"$ LoadCase : ",LoadSetS,"sin"
			if(ratio==1.0)	then
				write(333,'(a)')"$ Unit     :   Pa"
				write(334,'(a)')"$ Unit     :   Pa"
			else
				write(333,'(a)')"$ Unit     :   MPa"
				write(334,'(a)')"$ Unit     :   MPa"
			endif
			write(333,'(a,1x,f6.2)')"$ Beta     :",beta/pi*180.0
			write(334,'(a,1x,f6.2)')"$ Beta     :",beta/pi*180.0
			write(333,'(a,1x,f6.2)')"$ FREQ     :",OmeI0
			write(334,'(a,1x,f6.2)')"$ FREQ     :",OmeI0
		elseif(TRIM(ADJUSTL(PRTFEPRESS))=="ANSYS")	then
			if(TRIM(ADJUSTL(PRESSNAME))=="NO")	then
				call QINTTOSTR(IP,O,2)			
				open(333,file=trim(adjustl(OutAccess))//'\'//'EPBeta'//B//'Ome'//O//'C.dat')
				open(334,file=trim(adjustl(OutAccess))//'\'//'EPBeta'//B//'Ome'//O//'S.dat')
			elseIF(TRIM(ADJUSTL(PRESSNAME))=="YES") THEN
				open(333,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(LCName(2*IP-1,IB,IU)))//'.dat')
				open(334,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(LCName(2*IP,IB,IU)))//'.dat')
			endif
		endif
		!修改压力符号：by LiHui 2014.04.18
		do i=1,NFEPoint
		    if(FEM_N==0)    then
			    ppresure(1,I)=ppall(1,I)*1000/(ratio*ratio)  !*(-1.0) !去掉*(-1.0)  By LiHui 2014.04.20
			    ppresure(2,I)=ppall(2,I)*1000/(ratio*ratio)  !*(-1.0) !去掉*(-1.0)  By LiHui 2014.04.20
			else
			    ppresure(1,I)=-ppall(1,I)*1000/(ratio*ratio)  !*(-1.0) !去掉*(-1.0)  By LiHui 2014.04.20
			    ppresure(2,I)=-ppall(2,I)*1000/(ratio*ratio)  !*(-1.0) !去掉*(-1.0)  By LiHui 2014.04.20		
			endif		
		end do
		if(TRIM(ADJUSTL(PRTFEPRESS))=="PATRAN")	then
			call	FormatTrans(NFEpoint,ppresure(1,:),CharP)  !去掉压力前面的负号  By LiHui 2014.04.20
			do	i=1,	NFEPoint
				WRITE(333,'(A6,2X,3x,I3,2x,A8,1x,I7)')'PLOAD2',LoadSetC,charp(i),WetEleID(i)
!				WRITE(333,'(A6,2X,3x,I3,2x,f8.2,I8)')'PLOAD2',LoadSetC,ppresure(1,I),WetEleID(i)
			end	do
			call	FormatTrans(NFEpoint,ppresure(2,:),CharP)  !去掉压力前面的负号  By LiHui 2014.04.20
			do	i=1,	NFEPoint
				WRITE(334,'(A6,2X,3x,I3,2x,A8,1x,I7)')'PLOAD2',LoadSetS,charp(i),WetEleID(i)
!				WRITE(334,'(A6,2X,3x,I3,2x,f8.2,I8)')'PLOAD2',LoadSetS,ppresure(2,I),WetEleID(i)
			end	do
			
			WRITE(333,'(a)')	trim(adjustl(Endchar))
			WRITE(334,'(a)')	trim(adjustl(Endchar))
			close(333);close(334)
		elseif(TRIM(ADJUSTL(PRTFEPRESS))=="ANSYS")	then
			write(333,'(1x,i8)')	NFEpoint
			do	i=1,	NFEPoint
				WRITE(333,'(2X,I8,2x,e12.4)')WetEleID(i),ppresure(1,I)
			end	do
			write(334,'(1x,i8)')	NFEpoint
			do	i=1,	NFEPoint
				WRITE(334,'(2X,I8,2x,e12.4)')WetEleID(i),ppresure(2,I)
			end	do
			WRITE(333,'(a)')'\End_of_File\'
			WRITE(334,'(a)')'\End_of_File\'
			close(333);close(334)
		endif
	elseif(TRIM(ADJUSTL(ansMethod))=="YES")	then
		LoadSetC=((IU-1)*NBHEAD+(IB-1))*NBOME+IP
		if(TRIM(ADJUSTL(PRTFEPRESS))=="PATRAN")	then
			if(TRIM(ADJUSTL(PRESSNAME))=="NO")	then
				call QINTTOSTR(IP,O,2)
				open(333,file=trim(adjustl(OutAccess))//'\'//'EPBeta'//B//'Ome'//O//'.bdf')	
			elseIF(TRIM(ADJUSTL(PRESSNAME))=="YES") THEN
				open(333,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(LCName(IP,IB,IU)))//'.bdf')			
			endif
			write(333,'(a,I3,2x)')"$ LoadCase : ",LoadSetC
			if(ratio==1.0)	then
				write(333,'(a)')"$ Unit     :   Pa"
			else
				write(333,'(a)')"$ Unit     :   MPa"
			endif
			write(333,'(a,1x,f6.2)')"$ Beta     :",beta/pi*180.0
			write(333,'(a,1x,f6.2)')"$ FREQ     :",OmeI0
		elseif(TRIM(ADJUSTL(PRTFEPRESS))=="ANSYS")	then
			if(TRIM(ADJUSTL(PRESSNAME))=="NO")	then
				call QINTTOSTR(IP,O,2)
				open(333,file=trim(adjustl(OutAccess))//'\'//'EPBeta'//B//'Ome'//O//'.dat')	
			else IF(TRIM(ADJUSTL(PRESSNAME))=="YES") THEN
				open(333,file=trim(adjustl(OutAccess))//'\'//trim(adjustl(LCName(IP,IB,IU)))//'.dat')			
			endif		
		endif
		do i=1,NFEPoint
		
			call	comp_to_ampha(ppall(1,I),ppall(2,I),temp(1,I),temp(2,I))
			
			if(FEM_N==0)    then
			    temp(1,I)=temp(1,I)*1000/(ratio*ratio)
			else
			    temp(1,I)=-temp(1,I)*1000/(ratio*ratio)
			endif
			
			temp(1,I)=temp(1,I)*cos((temp(2,I)+WavePhase)/180.0*PI)  !*(-1.0)	!去掉*(-1.0)  By LiHui 2014.04.20	
				
		end do
		if(TRIM(ADJUSTL(PRTFEPRESS))=="PATRAN")	then
			call	FormatTrans(NFEpoint,temp(1,:),CharP)
			do	i=1,	NFEPoint
				WRITE(333,'(A6,2X,3x,I3,2x,A8,1x,I7)')'PLOAD2',LoadSetC,charp(i),WetEleID(i)
			end	do
			WRITE(333,'(a)')	trim(adjustl(Endchar))
			close(333)
		elseif(TRIM(ADJUSTL(PRTFEPRESS))=="ANSYS")	then
			write(333,'(1x,I8)')NFEpoint
			do	i=1,	NFEPoint
				WRITE(333,'(2X,I8,2x,e12.4)')WetEleID(i),temp(1,I)
			end	do
			WRITE(333,'(a)')'\End_of_File\'
			close(333)
		endif
	endif
! enddo
end subroutine ep_putouts
!----------------------------------------------------------------------

!------------------------------------------------------
subroutine	FormatTrans(N,P,CharP)
implicit	none
integer	::	N
real(8),	dimension(1:N)	::	P
character(len=10),	dimension(1:N)	::	CharP
character(len=10)	::	tempchar
integer(4)	::	i,j,k

open(301,file='tempp.dat')
do	i=1,	N
	write(301,'(E10.3)')	P(i)
end	do
close(301)

open(301,file='tempp.dat')
do	i=1,	N
	read(301,'(A10)')	tempchar
	do	j=1,	8	
		if(tempchar(j:j)=='E')	then
			if(tempchar(j+1:j+1)/='-')	then
				tempchar(j:j)='+'
				if(tempchar(j+2:j+2)=='0')	then
					do	k=j+1,8
						tempchar(k:k)=tempchar(k+2:k+2)	
					end	do
					tempchar(9:10)=' '				
				else
					do	k=j+1,9
						tempchar(k:k)=tempchar(k+1:k+1)	
					end	do
				endif
			else
				tempchar(j:j)='-'
				if(tempchar(j+2:j+2)=='0')	then
					do	k=j+1,8
						tempchar(k:k)=tempchar(k+2:k+2)	
					end	do
					tempchar(10:10)=' '
				else
					do	k=j+1,9
						tempchar(k:k)=tempchar(k+1:k+1)	
					end	do
				end	if
			end	if
		end	if
	end	do
	CharP(i)=tempchar
end	do
close(301)
end	subroutine FormatTrans
!----------------------------------------------------------------



!************************************************************************************
!************************************************************************************
!程序功能：由面元控制点处压力线性插值计算给定点压力
!注:首尾处剖面由于横坐标不同，插值有一些误差      
!************************************************************************************
subroutine tranpress3(xav0,ps,pw,pd,pr,pall,cpps,cppw,cppd,cppr,cppall)
implicit none

real(8),intent(in),dimension(1:3,1:NB)::xav0
real(8),intent(in),dimension(1:2,1:Factor*NB)::ps,pw,pd,pr,pall
real(8),intent(out),dimension(1:2,1:NFEPoint)::cpps,cppw,cppd,cppr,cppall
real(8),dimension(1:2,1:2)::tempps,temppw,temppd,temppr,temppall
real(8),allocatable,dimension(:,:)::x00
integer::i,j,IP
allocate(x00(1:3,1:Factor*NB))
if(Factor==1)then
 x00=xav0
 else if(Factor==2)then

 do i=1,2*NB/ArcNum
    do j=1,ArcNum/2

            !x坐标赋值
            x00(1,j+(i-1)*ArcNum)=xav0(1,j+(i-1)*ArcNum/2)
            x00(1,j+(i-1)*ArcNum+ArcNum/2)=xav0(1,(i-1)*ArcNum/2+ArcNum/2+1-j)

            !y坐标赋值
            x00(2,j+(i-1)*ArcNum)=xav0(2,j+(i-1)*ArcNum/2)
            x00(2,j+(i-1)*ArcNum+ArcNum/2)=-xav0(2,(i-1)*ArcNum/2+ArcNum/2+1-j)
            !z坐标赋值
            x00(3,j+(i-1)*ArcNum)=xav0(3,j+(i-1)*ArcNum/2)
            x00(3,j+(i-1)*ArcNum+ArcNum/2)=xav0(3,(i-1)*ArcNum/2+ArcNum/2+1-j)
end do
end do
elseif(Factor==4)then
do i=1,2*NB/ArcNum
    do j=1,ArcNum/2
            !x坐标赋值
            x00(1,j+(i-1)*ArcNum)=xav(1,j+(i-1)*ArcNum/2)
            x00(1,j+(i-1)*ArcNum+ArcNum/2)=-xav(1,(i-1)*ArcNum/2+ArcNum/2+1-j)
            x00(1,2*NB+j+(i-1)*ArcNum)=xav(1,j+NB-i*ArcNum/2)
            x00(1,2*NB+(i-1)*ArcNum+ArcNum/2+j)=xav(1,NB-(i-1)*ArcNum/2+1-j)

            !y坐标赋值
            x00(2,j+(i-1)*ArcNum)=xav(2,j+(i-1)*ArcNum/2)
            x00(2,j+(i-1)*ArcNum+ArcNum/2)=-xav(2,(i-1)*ArcNum/2+ArcNum/2+1-j)
            x00(2,2*NB+j+(i-1)*ArcNum)=xav(2,j+NB-i*ArcNum/2)
            x00(2,2*NB+(i-1)*ArcNum+ArcNum/2+j)=-xav(2,NB-(i-1)*ArcNum/2+1-j)

            !z坐标赋值
            x00(3,j+(i-1)*ArcNum)=xav(3,j+(i-1)*ArcNum/2)
            x00(3,j+(i-1)*ArcNum+ArcNum/2)=xav(3,(i-1)*ArcNum/2+ArcNum/2+1-j)
            x00(3,2*NB+j+(i-1)*ArcNum)=xav(3,j+NB-i*ArcNum/2)
            x00(3,2*NB+(i-1)*ArcNum+ArcNum/2+j)=xav(3,NB-(i-1)*ArcNum/2+1-j)
end do
end do
end if
do i=1,NFEPoint
		IP=0
		do j=1,Factor*NB/ArcNum
			if(FEMPOINT(1,i)<x00(1,(j-1)*ArcNum+1))then
				IP=j
				exit
			endif
		enddo
		if(IP==0)IP=Factor*NB/ArcNum+1
		if(IP/=1.and.IP/=Factor*NB/ArcNum+1)then
			if(FEMPOINT(2,i)>0.0)then !左舷
				!IP-2剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pall(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,1),1,pall(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									ps(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,1),1,ps(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pw(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,1),1,pw(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pd(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,1),1,pd(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pr(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,1),1,pr(j,(IP-2)*ArcNum+ArcNum/2+1))
				enddo
				!IP-1剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pall(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,2),1,pall(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									ps(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,2),1,ps(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pw(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,2),1,pw(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pd(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,2),1,pd(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pr(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,2),1,pr(j,(IP-1)*ArcNum+ArcNum/2+1))
				enddo
			elseif(FEMPOINT(2,i)<=0.0)then !右舷
				!IP-2剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pall(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,1),0,pall(j,(IP-2)*ArcNum+ArcNum/2+1-1))

					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									ps(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,1),0,ps(j,(IP-2)*ArcNum+ArcNum/2+1-1))

					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pw(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,1),0,pw(j,(IP-2)*ArcNum+ArcNum/2+1-1))

					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pd(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,1),0,pd(j,(IP-2)*ArcNum+ArcNum/2+1-1))

					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pr(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,1),0,pr(j,(IP-2)*ArcNum+ArcNum/2+1-1))

				enddo
				!IP-1剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pall(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,2),0,pall(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									ps(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,2),0,ps(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pw(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,2),0,pw(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pd(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,2),0,pd(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pr(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,2),0,pr(j,(IP-1)*ArcNum+ArcNum/2+1-1))
				enddo
			endif
			cpps(1:2,i)=(tempps(1:2,2)-tempps(1:2,1))*(FEMPOINT(1,i)-x00(1,(IP-2)*ArcNum+1))/(x00(1,(IP-1)*ArcNum+1)-x00(1,(IP-2)*ArcNum+1))+tempps(1:2,1)
			cppw(1:2,i)=(temppw(1:2,2)-temppw(1:2,1))*(FEMPOINT(1,i)-x00(1,(IP-2)*ArcNum+1))/(x00(1,(IP-1)*ArcNum+1)-x00(1,(IP-2)*ArcNum+1))+temppw(1:2,1)
			cppd(1:2,i)=(temppd(1:2,2)-temppd(1:2,1))*(FEMPOINT(1,i)-x00(1,(IP-2)*ArcNum+1))/(x00(1,(IP-1)*ArcNum+1)-x00(1,(IP-2)*ArcNum+1))+temppd(1:2,1)
			cppr(1:2,i)=(temppr(1:2,2)-temppr(1:2,1))*(FEMPOINT(1,i)-x00(1,(IP-2)*ArcNum+1))/(x00(1,(IP-1)*ArcNum+1)-x00(1,(IP-2)*ArcNum+1))+temppr(1:2,1)
			cppall(1:2,i)=(temppall(1:2,2)-temppall(1:2,1))*(FEMPOINT(1,i)-x00(1,(IP-2)*ArcNum+1))/(x00(1,(IP-1)*ArcNum+1)-x00(1,(IP-2)*ArcNum+1))+temppall(1:2,1)
		elseif(IP==1)then
			if(FEMPOINT(2,i)>0.0)then !左舷
				!IP-1剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pall(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,2),1,pall(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									ps(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,2),1,ps(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pw(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,2),1,pw(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pd(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,2),1,pd(j,(IP-1)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),x00(3,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),&
									pr(j,(IP-1)*ArcNum+1:(IP-1)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,2),1,pr(j,(IP-1)*ArcNum+ArcNum/2+1))
				enddo
			elseif(FEMPOINT(2,i)<=0.0)then !右舷
				!IP-1剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pall(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,2),0,pall(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									ps(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,2),0,ps(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pw(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,2),0,pw(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pd(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,2),0,pd(j,(IP-1)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),x00(3,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),&
									pr(j,(IP-1)*ArcNum+ArcNum/2+1:IP*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,2),0,pr(j,(IP-1)*ArcNum+ArcNum/2+1-1))
				enddo
			endif
			cpps(1:2,i)=tempps(1:2,2)
			cppw(1:2,i)=temppw(1:2,2)
			cppd(1:2,i)=temppd(1:2,2)
			cppr(1:2,i)=temppr(1:2,2)
			cppall(1:2,i)=temppall(1:2,2)
		elseif(IP==Factor*NB/ArcNum+1)then
			if(FEMPOINT(2,i)>0.0)then !左舷
				!IP-2剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pall(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,1),1,pall(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									ps(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,1),1,ps(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pw(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,1),1,pw(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pd(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,1),1,pd(j,(IP-2)*ArcNum+ArcNum/2+1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),x00(3,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),&
									pr(j,(IP-2)*ArcNum+1:(IP-2)*ArcNum+ArcNum/2),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,1),1,pr(j,(IP-2)*ArcNum+ArcNum/2+1))
				enddo
			elseif(FEMPOINT(2,i)<=0.0)then !右舷
				!IP-2剖面
				do j=1,2
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pall(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppall(j,1),0,pall(j,(IP-2)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									ps(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),tempps(j,1),0,ps(j,(IP-2)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pw(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppw(j,1),0,pw(j,(IP-2)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pd(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppd(j,1),0,pd(j,(IP-2)*ArcNum+ArcNum/2+1-1))
					call LinDiagPola(ArcNum/2,x00(2,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),x00(3,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),&
									pr(j,(IP-2)*ArcNum+ArcNum/2+1:(IP-1)*ArcNum),FEMPOINT(2,i),FEMPOINT(3,i),temppr(j,1),0,pr(j,(IP-2)*ArcNum+ArcNum/2+1-1))
				enddo
			endif
			cpps(1:2,i)=tempps(1:2,1)
			cppw(1:2,i)=temppw(1:2,1)
			cppd(1:2,i)=temppd(1:2,1)
			cppr(1:2,i)=temppr(1:2,1)
			cppall(1:2,i)=temppall(1:2,1)
		endif
enddo
deallocate(x00)
end subroutine tranpress3
!---------------------------------------------------------------------------------------------------------------------------


!************************************************************************************
!程序功能：剖面压力局部对角线坐标系线性插值
!输入参数：N为船体半剖面点数，y0(1:n),z0(1:n)为半剖面坐标，pa(1:n)为压力值
!          y,z为压力计算点坐标
!输出参数：pa为计算点压力值
!************************************************************************************
subroutine LinDiagPola(n,y0,z0,pa0,y,z,pa,IS_Stern,pa_IS)
implicit none
integer(4),intent(in)::n,IS_Stern
real(8),intent(in),dimension(1:n)::y0,z0,pa0
real(8),intent(in)::y,z,pa_IS
real(8),intent(out)::pa
real(8),dimension(1:3,1:3)::ee !方向余弦
real(8),dimension(1:3)::xyz0,xyz00,tran_xyz
real(8),allocatable,dimension(:,:)::tran_xyz0
integer(4)::i
ee=0.0
ee(1,1)=1.0
ee(2,2)=y0(n)-y0(1);ee(3,2)=z0(n)-z0(1);ee(:,2)=ee(:,2)/sqrt(sum(ee(:,2)**2))
ee(2,3)=ee(3,1)*ee(1,2)-ee(1,1)*ee(3,2)
ee(3,3)=ee(1,1)*ee(2,2)-ee(2,1)*ee(1,2)
ee(:,3)=ee(:,3)/sqrt(sum(ee(:,3)**2))
allocate(tran_xyz0(1:3,1:n))
xyz0=0.0
xyz00=0.0;xyz00(2)=y0(1);xyz00(3)=z0(1)
do i=1,n
	xyz0(2)=y0(i)
	xyz0(3)=z0(i)
	tran_xyz0(2,i)=dot_product(ee(:,2),xyz0-xyz00)
	tran_xyz0(3,i)=dot_product(ee(:,3),xyz0-xyz00)
enddo
xyz0=0.0;xyz0(2)=y;xyz0(3)=z
tran_xyz(2)=dot_product(ee(:,2),xyz0-xyz00)
tran_xyz(3)=dot_product(ee(:,3),xyz0-xyz00)
!插值
do i=1,n
	if(tran_xyz(2)<tran_xyz0(2,i))then
		if(i==1)then
			if(IS_Stern==0)then
				pa=(pa0(1)+pa_IS)/2.0
			else
				pa=(pa0(2)-pa0(1))*(tran_xyz(2)-tran_xyz0(2,1))/(tran_xyz0(2,2)-tran_xyz0(2,1))+pa0(1)
			endif
			exit
		else
			pa=(pa0(i)-pa0(i-1))*(tran_xyz(2)-tran_xyz0(2,i-1))/(tran_xyz0(2,i)-tran_xyz0(2,i-1))+pa0(i-1)
			exit
		endif
	endif
	if(i==n)then
		if(IS_Stern==1)then
			pa=(pa0(n)+pa_IS)/2.0
		else
			pa=(pa0(n)-pa0(n-1))*(tran_xyz(2)-tran_xyz0(2,n-1))/(tran_xyz0(2,n)-tran_xyz0(2,n-1))+pa0(n-1)
		endif
	endif
enddo
deallocate(tran_xyz0)
end subroutine LinDiagPola
!----------------------------------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------------------------
subroutine	Precal	(pps,ppw,ppd,ppr,ppall,cpps,cppw,cppd,cppr,cppall)	
implicit	none
real(8),	intent(in)	::	PPs(1:2,1:Factor,1:NB),PPw(1:2,1:Factor,1:NB),PPd(1:2,1:Factor,1:NB),PPr(1:2,1:Factor,1:NB),PPall(1:2,1:Factor,1:NB)
real(8),	intent(out)	::	CPPs(1:2,1:NFEpoint),CPPw(1:2,1:NFEpoint),CPPd(1:2,1:NFEpoint),CPPr(1:2,1:NFEpoint),CPPall(1:2,1:NFEpoint)
integer	::	position(1:NFEpoint),L(1:NFEpoint)
integer	::	i
L=0
do	i=1,	NFEPoint

	call	PointPosition(FEMPOINT(:,i),position(i),L(i),normal(:,i))

	if(position(i)>0)	then
	       
        if(abs(Factor-2)<0.01.and.(L(i)>2)) L(i)=2
   

		cppS(:,I)=PPS(:,L(i),Position(i))
		cppW(:,I)=PPW(:,L(i),Position(i))
		cppD(:,I)=PPD(:,L(i),Position(i))
		cppR(:,I)=PPR(:,L(i),Position(i))
		cppALL(:,I)=PPALL(:,L(i),Position(i))	
	else
		cppS(:,I)=0.0
        cppW(:,I)=0.0
		cppD(:,I)=0.0
		cppR(:,I)=0.0
        cppALL(:,I)=0.0
        
	endif	
end	do
!----------------test------------------
!open(3,file='position.dat')
!do	i=1,NFEpoint
!	write(3,*)i,position(i),L(i)
!enddo
!close(3)
!--------------------------------------
end	subroutine PreCal
!----------------------------------------------------------------------------------------------------------------------



!---------------------------------------------------------------------------------
subroutine PointPosition(xyz0,position,AreaNum,normal0)
implicit	none

real(8),	intent(in)	::	xyz0(1:3)
real(8),	intent(in)	::	normal0(1:3)
integer,INTENT(OUT)	::	position
real(8)	::	xq0(1:2,1:NB)			!面元顶点局部坐标
! real(8)	::	Pxy(1:2)
real(8),allocatable::x00(:,:),ee(:,:,:),xqq(:,:,:)
real(8)	::	xyz(1:3)
real(8)	::	x1(1:2,1:6)
real(8)	::	ss(1:5)
real(8)	::	test1(1:NB)
real(8)::  LP(1:NB)
real(8)	::	esp,esp1,esp2,test,t
real(8)	::	min
real(8)	::	inclination(1:NB),d(1:NB)
real(8)	::	p1,p2,w1,w2,ea1,ea2,ea0
real(8)	::	temp
integer	::	i,j,k,AreaNum,L
!---------------------------------------------------------------
allocate(x00(1:3,1:NB),ee(1:3,1:3,1:NB),xqq(1:2,1:4,1:NB))
!=========================================================================
esp     =   1.0e-1
esp1    =   1.0e-6
esp2    =   1.0e-1
position=   0

do L=1,TarNum,Space

    if(AreaNum>0)   exit
    
    select case(L)
        case    (1)
            x00=xav
            ee=e
            xqq=xq            
        case    (2)
        
            x00(1,:)    =   -xav(1,:)
            x00(2,:)    =    xav(2,:)   
            ee(1,3,:)   =   -e(1,3,:)
            ee(2,3,:)   =    e(2,3,:)
            ee(3,3,:)   =    e(3,3,:)    
            ee(1,1,:)   =   -e(1,1,:)
            ee(2,1,:)   =    e(2,1,:)
            ee(3,1,:)   =    e(3,1,:)
            ee(1,2,:)   =    e(1,2,:)
            ee(2,2,:)   =   -e(2,2,:)
            ee(3,2,:)   =   -e(3,2,:)
            xqq(1,:,:)  =    xq(1,:,:)
            xqq(2,:,:)  =   -xq(2,:,:)
            
       case     (3)
       
            x00(1,:)    =   -xav(1,:)
            x00(2,:)    =   -xav(2,:)
            ee(1,3,:)   =   -e(1,3,:)
            ee(2,3,:)   =   -e(2,3,:)
            ee(3,3,:)   =    e(3,3,:)
            ee(1,1,:)   =   -e(1,1,:)
            ee(2,1,:)   =   -e(2,1,:)
            ee(3,1,:)   =    e(3,1,:)
            ee(1,2,:)   =   -e(1,2,:)
            ee(2,2,:)   =   -e(2,2,:)
            ee(3,2,:)   =    e(3,2,:)
            xqq(1,:,:)  =    xq(1,:,:)
            xqq(2,:,:)  =    xq(2,:,:)
            
        case    (4)
        
            x00(1,:)    =    xav(1,:)
            x00(2,:)    =   -xav(2,:)
            ee(1,3,:)   =    e(1,3,:)
            ee(2,3,:)   =   -e(2,3,:)
            ee(3,3,:)   =    e(3,3,:)
            ee(1,1,:)   =    e(1,1,:)
            ee(2,1,:)   =   -e(2,1,:)
            ee(3,1,:)   =    e(3,1,:)
            ee(1,2,:)   =   -e(1,2,:)
            ee(2,2,:)   =    e(2,2,:)
            ee(3,2,:)   =   -e(3,2,:)
            xqq(1,:,:)  =    xq(1,:,:)
            xqq(2,:,:)  =   -xq(2,:,:)
    end select
    do	i=1,	NB
    
	    test=abs(xyz0(1)-x00(1,i))+abs(xyz0(2)-x00(2,i))+abs(xyz0(3)-xav(3,i))
	    	
	    if(abs(test)<=ESP)	then
	    
		    position=i 
		    AreaNum=L
		    exit
		    
	    else
	    
		    t=-(ee(1,3,i)*(xyz0(1)-x00(1,i))+&
			    ee(2,3,i)*(xyz0(2)-x00(2,i))+&
			    ee(3,3,i)*(xyz0(3)-xav(3,i)))
		    xyz(:)=ee(:,3,i)*t+xyz0(:)
		    
		    do	j=1,	2
			    xq0(j,i)=dot_product(ee(:,j,i),xyz(:)-x00(:,i))
		    end	do
		    
		    
!		    LP(i)   =   sqrt(sum((xqq(:,1,i)-xqq(:,2,i))**2))   +   &
!		                sqrt(sum((xqq(:,2,i)-xqq(:,3,i))**2))   +   &
!		                sqrt(sum((xqq(:,3,i)-xqq(:,4,i))**2))   +   &
!		                sqrt(sum((xqq(:,4,i)-xqq(:,1,i))**2))
!		                
!		    LP(i)   =   LP(i)/4.0/sqrt(NFEPoint*1.0/NB)    
		  
		    
		    d(i)    =	sqrt((xyz0(1)-xyz(1))**2+(xyz0(2)-xyz(2))**2+(xyz0(3)-xyz(3))**2)    
		    	
		    	
	        !计算面元面积
		    ea0=0.0
		    
		    do j=1,4
			    if(j==1)	then
				    ss(1)=sqrt(sum((xqq(:,1,i)-xq0(:,i))**2))
				    ss(2)=sqrt(sum((xqq(:,2,i)-xq0(:,i))**2))
				    ss(3)=sqrt(sum((xqq(:,1,i)-xqq(:,2,i))**2))
			    elseif(j==2)	then
				    ss(1)=sqrt(sum((xqq(:,2,i)-xq0(:,i))**2))
				    ss(2)=sqrt(sum((xqq(:,3,i)-xq0(:,i))**2))
				    ss(3)=sqrt(sum((xqq(:,2,i)-xqq(:,3,i))**2))
			    elseif(j==3)	then
				    ss(1)=sqrt(sum((xqq(:,3,i)-xq0(:,i))**2))
				    ss(2)=sqrt(sum((xqq(:,4,i)-xq0(:,i))**2))
				    ss(3)=sqrt(sum((xqq(:,3,i)-xqq(:,4,i))**2))
			    else
				    ss(1)=sqrt(sum((xqq(:,4,i)-xq0(:,i))**2))
				    ss(2)=sqrt(sum((xqq(:,1,i)-xq0(:,i))**2))
				    ss(3)=sqrt(sum((xqq(:,4,i)-xqq(:,1,i))**2))
			    endif
			    
			    p1=0.5*(ss(1)+ss(2)+ss(3))
			    w1=p1*(p1-ss(1))*(p1-ss(2))*(p1-ss(3))
			    ea1=0.0	
			    		
			    if(w1>esp1)	then
				    ea1=sqrt(w1)
			    else
				    ea1=0.0	
			    endif		
			    ea0=ea0+ea1
			    
		    enddo
	    
		    
		    test1(i)=abs(1-ea0/ea(i))
		    
		    temp=dot_product(normal0(1:3),ee(:,3,i))
		    
		    if(temp<-1.0)	then
			    temp=-1.0
		    elseif(temp>1.0)	then
			    temp=1.0
		    endif

			inclination=acos(temp)/PI*180
			
!		    if(test1(i)<=esp.and.abs(inclination(i))<=80.0.and.(d(i)/LP(i))<=10*esp2)	then 	
	    
		    if(test1(i)<=5.0e-2.and.abs(inclination(i))<=80.0.and.d(i)<=SB/20.0)	then		    		    		    
			    position=i	
			    AreaNum=L
			    exit			    		
		    endif		
	    end	if		
    end	do
    

    
end do
! if(position>0)	Pxy(1:2)=xq0(:,position)
return
end subroutine PointPosition
!------------------------------------------------------------------



!---------------------------------------------------------------------------------------------------------
!NPOINT有限元网格数；WetEleID(1:NFEPOINT)网格单元号；coor(1:3,1:NFEPOINT)随船平动坐标系下有限元网格中心点坐标
!normal(1:3,1:NFEPOINT)有限元网格单位法向量(指向船体外部)
!---------------------------------------------------------------------------------------------------------
subroutine	GetPpoint(BDFINPUT,coor,WetEleID0,tempcoor,ratio,normal)		
implicit	none

character(len=300),INTENT(IN)::BDFINPUT
real(8),INTENT(IN),dimension(1:3)::tempcoor
real,INTENT(IN)::ratio
!integer   ::  FEM_N
! INTEGER(4),INTENT(OUT)::NwElm
integer(4),intent(out),dimension(1:1000000)::WetEleID0
! character(len=100),INTENT(OUT)::Endchar
real(8),INTENT(OUT),dimension(1:3,1:100000)	::	coor,Normal

character(len=100)	::	char1,char2
character(len=6)	::	temp1,temp2
character(len=6)	::	temp3,temp4
character(len=1),	dimension(1:100)	::	Tempchar,CTempchar
integer(4),	dimension(1:10)		::	position1,position2
integer(4),	dimension(1:1000000)		::	EleID,NodeID,WetNodeID	
integer(4),	dimension(1:1000000,1:4)	::	EleNode
real(8),	dimension(1:1000000,1:3)	::	tempNodecoor,Nodecoor
real(8),	allocatable,	dimension(:,:)	::	xavv,ee
real(8),	dimension(1:3)	::	r13,r24
integer(4)	::	NC,NG,NT,lenth1,lenth2,StartNO,EndNO,flag,Nlap,NW
integer(4)	::	i,j,k,m	,jj
!------以下为程序代码-----------------------------
NC=0;NG=0;NT=0

open(1001,file=BDFINPUT)                    	!打开给定的工程
open(1002,file='TEMPC.DAT')						!存放四边形单元的节点号
open(1003,file='TempT.dat')						!存放三角形单元的节点号
open(1004,file='TEMPG.DAT')						!存放网格节点的坐标
do
	read(1001,'(A)')	char1
	if(char1(1:7)=='ENDDATA')	then
		Endchar=char1
		exit     !读到文件末尾，跳出循环
	endif
	if(char1(1:8)=='CQUAD4  ')		then
		lenth1=len(trim(adjustl(char1)))   !所读文件此行的长度
        DO I=1, 7
              write(1002,'(a)')char1(((i-1)*8+1):i*8)
        ENDDO
		NC=NC+1    !四边形单元数
	elseif(char1(1:8)=='CQUAD4* ')	then
		lenth1=len(trim(adjustl(char1)))
		Nlap=int(lenth1/8)
		do	i=1,	Nlap
			if(char1(((i-1)*8+1):i*8)/="        ")	then
				write(1002,'(a)')char1(((i-1)*8+1):i*8)
			endif
		end	do
		write(1002,'(a)')char1((Nlap*8+1):lenth1)
		read(1001,'(1x,a)')	char2
		write(1002,'(a)')	char2		
		NC=NC+1    !四边形单元的个数
	elseif(char1(1:8)=='CTRIA3  ')	then
		lenth1=len(trim(adjustl(char1)))
        DO I=1, 6
              write(1003,'(a)') char1(((i-1)*8+1):i*8)
        ENDDO     
		NT=NT+1	 !三角形单元数
	elseif(char1(1:8)=='CTRIA3* ')	then
		lenth1=len(trim(adjustl(char1)))
		Nlap=int(lenth1/8)
		do	i=1,	Nlap
			if(char1(((i-1)*8+1):i*8)/="        ")	then
				write(1003,'(a)')char1(((i-1)*8+1):i*8)
			endif
		end	do
		write(1003,'(a)')char1((Nlap*8+1):lenth1)
		read(1001,'(1x,a)')	char2
		write(1003,'(a)')	char2		
		NT=NT+1	 !三角形单元的个数
	elseif(char1(1:5)=='GRID    ')	then
		lenth1=len(trim(adjustl(char1)))
		Nlap=int(lenth1/8)
		NW=0
		do	i=1,	Nlap
			if(NW<=5)	then
				if(char1(((i-1)*8+1):i*8)/="        ")	then
					jj=0
					do	j=2,	8
						if(char1((i-1)*8+j:(i-1)*8+j)=='-'.and.char1((i-1)*8+j-1:(i-1)*8+j-1)/=' ') jj=j    !科学计数中的“-”jj为每八个字符中的第几个为"-" 如-9.6-11
					enddo
					if(jj/=0)	then
						tempchar(:)=' '
						do	j=1,	jj-1
							tempchar(j)=char1((i-1)*8+j:(i-1)*8+j)    !把科学计数中的小数部分输出到tempchar(1:jj-1)
						end	do
						tempchar(jj)="E"   !科学计数E
						do	j=jj,	8
							tempchar(j+1)=char1((i-1)*8+j:(i-1)*8+j)
						end	do
						write(1004,'(9(a))')(tempchar(j),j=1,9)    !科学计数，如-9.6E-11
						NW=NW+1
					else
						write(1004,'(a)')char1(((i-1)*8+1):i*8)
						NW=NW+1
					endif
				endif
			endif
		end	do
		if(NW<5)	then
			jj=0
			if(lenth1-8*Nlap>=1)	then
				do	j=2,	lenth1-8*Nlap
					if(char1(8*Nlap+j:8*Nlap+j)=="-".and.char1(8*Nlap+j-1:8*Nlap+j-1)/=' ') jj=j
				enddo
				if(jj/=0)	then
					tempchar(:)=' '
					do	j=1,	jj-1
						tempchar(j)=char1((i-1)*8+j:(i-1)*8+j)
					enddo
					tempchar(jj)="E"
					do	j=jj,	lenth1-8*Nlap
						tempchar(j+1)=char1((i-1)*8+j:(i-1)*8+j)
					enddo
					write(1004,*)(tempchar(j),j=1,9)
				else
					write(1004,'(a)')char1((Nlap*8+1):lenth1)
				endif		
			end	if
		end	if
		NG=NG+1   !点node的个数
	elseif(char1(1:5)=='GRID*   ')	then
		read(1001,'(1x,a)')char2
		char2=trim(adjustl(char2))
		lenth1=len(trim(adjustl(char1)))
		lenth2=len(trim(adjustl(char2)))
		Nlap=int(lenth1/8)
		NW=0
	
		do	i=1,	Nlap
			if(NW<=5)	then
				if(char1(((i-1)*8+1):i*8)/="        ")	then
					jj=0
					do	j=2,	8
						if(char1((i-1)*8+j:(i-1)*8+j)=="-".and.char1((i-1)*8+j-1:(i-1)*8+j-1)/=" ") jj=j
					enddo
					if(jj/=0)	then
						tempchar(:)=' '
						do	j=1,	jj-1
							tempchar(j)=char1((i-1)*8+j:(i-1)*8+j)
						enddo
						tempchar(jj)="E"
						do	j=jj,	8
							tempchar(j+1)=char1((i-1)*8+j:(i-1)*8+j)
						enddo
						write(1004,'(9(a))')(tempchar(j),j=1,9)
						NW=NW+1
					else
						write(1004,'(a)')char1(((i-1)*8+1):i*8)
						NW=NW+1
					endif
				endif
			endif
		end	do
		jj=0
		if(NW<5)	then
			if(lenth1-8*Nlap>=1)	then
				do	j=2,	lenth1-8*Nlap
					if(char1(8*Nlap+j:8*Nlap+j)=="-".and.char1(8*Nlap+j-1:8*Nlap+j-1)/=" ") jj=j
				end	do
				if(jj/=0)	then
					tempchar(:)=' '
					do	j=1,	jj-1
						tempchar(j)=char1((i-1)*8+j:(i-1)*8+j)
					end	do
					tempchar(jj)="E"
					do	j=jj,	lenth1-8*Nlap
						tempchar(j+1)=char1((i-1)*8+j:(i-1)*8+j)
					end	do
					write(1004,*)(tempchar(j),j=1,9)
				else
					write(1004,'(a)')char1((Nlap*8+1):lenth1)
				endif
			end	if
		end	if
		jj=0
		do	j=2,	lenth2
			if(char2(j:j)=="-".and.char2(j-1:j-1)/=' ') jj=j
		enddo
		if(jj/=0)	then
			tempchar(:)=' '
			do	j=1,	jj-1
				tempchar(j)=char2(j:j)
			end	do
			tempchar(jj)="E"
			do	j=jj,	lenth2
				tempchar(j+1)=char2(j:j)
			end	do
			write(1004,*)(tempchar(j),j=1,lenth2+1)
		else
			write(1004,'(a)')char2(1:lenth2)
		endif
		NG=NG+1
	endif
end	do
close(1001)
close(1002)
close(1003)
close(1004)
open(1001,file='TEMPC.DAT')
call	IUTMP(1001)
do	i=1,	NC   !四边形单元号EleID(1:NC)
	read(1001,*)	temp1
	read(1001,*)	EleID(I)
	read(1001,*)	temp2
	do	j=1,	4
		read(1001,*)	EleNode(i,j)   !EleNode(1:NC,1:4)四边形单元的四个节点号
	end	do
end	do
close(1001)
open(1002,file='TempT.dat')
call	IUTMP(1002) 
do	i=NC+1,	NC+NT
	read(1002,*)	temp1
	read(1002,*)	EleID(I)   !EleID(NC+1:NC+NT)三角形单元号
	read(1002,*)	temp2
	do	j=1,	3
		read(1002,*)	EleNode(i,j)    !EleNode(NC+1:NC+NT,1:3)三角形单元的三个节点号
	enddo
	EleNode(i,4)=EleNode(i,3)
end	do
close(1002)
open(1003,file='TEMPG.DAT')
call	IUTMP(1003)
do	i=1,	NG
	read(1003,*)	temp3
	read(1003,*)	NodeID(i)  !节点号
	do	j=1,	3
		read(1003,*)	tempNodecoor(i,j)   !每个网格节点的坐标
	end	do
end	do
close(1003)
do	i=1,	NG
	Nodecoor(NodeID(i),:)=tempNodecoor(i,:)/ratio-tempcoor(:)   !节点坐标移动到随船平动坐标系下
	Nodecoor(NodeID(i),:)=Matmul(Ty,Nodecoor(NodeID(i),:))
end	do

allocate(xavv(1:NC+NT,1:3),ee(1:3,1:NC+NT))
xavv=0.0;coor=0.0
NFEPOINT=0
do	i=1,	NC+NT
	if(i<=NC)	then
		xavv(i,:)=(Nodecoor(EleNode(i,1),:)+Nodecoor(EleNode(i,2),:)+Nodecoor(EleNode(i,3),:)+Nodecoor(EleNode(i,4),:))/4.0   !四边形单元网格中心点坐标
	else
		xavv(i,:)=(Nodecoor(EleNode(i,1),:)+Nodecoor(EleNode(i,2),:)+Nodecoor(EleNode(i,3),:))/3.0    !三角形单元网格中心点坐标
	endif
	
	r13=Nodecoor(EleNode(i,3),:)-Nodecoor(EleNode(i,1),:)
	r24=Nodecoor(EleNode(i,4),:)-Nodecoor(EleNode(i,2),:)
	ee(1,i)=r13(2)*r24(3)-r13(3)*r24(2)
	ee(2,i)=r13(3)*r24(1)-r13(1)*r24(3)
	ee(3,i)=r13(1)*r24(2)-r13(2)*r24(1)	
	ee(:,i)=ee(:,i)/(sqrt(sum(ee(:,i)**2)))    !法向量单位化
	
	if(xavv(i,3)<=0.001)	then		 !平均湿表面
	
		NFEPOINT=NFEPOINT+1   !计算点的个数Npoint（计算的有限元网格数）
		WetEleID0(NFEPOINT)=EleID(i)   !网格单元号
		wetNodeID(NFEPOINT)=i   
		coor(:,NFEPOINT)=xavv(i,:)     !计算点坐标（计算的有限元网格中心点坐标-随船平动坐标系下）
		
		if(FEM_N==0)    then
		    normal(:,NFEPOINT)=ee(:,i)   !有限元计算网格的法向量-指向船体内部
		else
		    normal(:,NFEPOINT)=-ee(:,i)   !有限元计算网格的法向量-指向船体外部
		endif		
	end	if
end	do
close(1001)
close(1002)
close(1003)
close(1004)
deallocate(xavv)
open(1001,file='TEMPC.DAT',status='old',dispose='delete')
open(1002,file='TEMPG.DAT',status='old',dispose='delete')
open(1003,file='TEMPT.DAT',status='old',dispose='delete')

end	subroutine getppoint
!-----------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
subroutine	GetAnsPoint(DATINPUT,coor,WetEleID0,tempcoor,ratio,normal)	
!------------------------------------------------------------------------------------	
implicit	none

real,INTENT(IN)::ratio
character(len=300),INTENT(IN)::DATINPUT
real(8),INTENT(IN),dimension(1:3)::tempcoor
! INTEGER(4),INTENT(OUT)::NwElm
integer(4),intent(out),dimension(1:1000000)::WetEleID0
real(8),INTENT(OUT),dimension(1:3,1:100000)	::	coor,Normal

integer	::	NwElm
integer,	allocatable	::	EleID(:)
integer,	allocatable	::	EleType(:)
integer,	allocatable	::	NodeID(:,:)
real(8),	allocatable	::	Nodecoor(:,:,:)
real(8),	allocatable,	dimension(:,:)	::	xavv,ee
real(8),	dimension(1:3)	::	r13,r24
integer(4)	::	i,j,k
!------以下为程序代码-----------------------------
open(1001,file=DATINPUT)	!打开给定的工程
call	IUTMP(1001)
read(1001,*)NwElm
allocate(EleID(1:NwElm),EleType(1:NwElm),NodeID(1:NwElm,1:4),Nodecoor(1:NwElm,1:4,1:3))
EleType=0
do	i=1,	NwElm
	read(1001,*)	EleID(I),(NodeID(i,j),j=1,4),(Nodecoor(i,1,j),j=1,3),(Nodecoor(i,2,j),j=1,3),	&
					(Nodecoor(i,3,j),j=1,3),(Nodecoor(i,4,j),j=1,3)
	if(sum(abs(Nodecoor(i,3,1:3)-Nodecoor(i,4,1:3)))<=1.0e-6)	EleType(i)=1
end	do
close(1001)
do	i=1,	NwElm
	do	j=1,	4
		Nodecoor(i,j,:)=Nodecoor(i,j,:)/ratio-tempcoor(:)
	end	do
end	do
allocate(xavv(1:NwElm,1:3),ee(1:3,1:NwElm))
xavv=0.0;coor=0.0
NFEPOINT=0
do	i=1,	NwElm
	if(EleType(i)==0)	then
		xavv(i,:)=(Nodecoor(i,1,:)+Nodecoor(i,2,:)+Nodecoor(i,3,:)+Nodecoor(i,4,:))/4.0
	else
		xavv(i,:)=(Nodecoor(i,1,:)+Nodecoor(i,2,:)+Nodecoor(i,3,:))/3.0
	endif
	r13=Nodecoor(i,1,:)-Nodecoor(i,3,:)
	r24=Nodecoor(i,2,:)-Nodecoor(i,4,:)
	ee(1,i)=r13(2)*r24(3)-r13(3)*r24(2)
	ee(2,i)=r13(3)*r24(1)-r13(1)*r24(3)
	ee(3,i)=r13(1)*r24(2)-r13(2)*r24(1)	
	ee(:,i)=ee(:,i)/(sqrt(sum(ee(:,i)**2))) !单位化
	if(xavv(i,3)<=0.0)	then		
		NFEPOINT=NFEPOINT+1
		WetEleID0(NFEPOINT)=EleID(i)
		coor(:,NFEPOINT)=xavv(i,:)
		normal(:,NFEPOINT)=ee(:,i)
	end	if
end	do
close(1001)
deallocate(xavv)

end	subroutine GetAnsPoint
!==============================================================================================================================




!-----------------------------------------------------------------------
END MODULE FEMPRESS_MOD