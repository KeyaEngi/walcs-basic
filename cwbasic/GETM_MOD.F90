!************************************************************************************
!程序功能：计算船舶质量矩阵(包括整船和到各站剖面)
!BY SUN：加入修正重力浮力差，重心与浮心纵向位置不一致部分
!***********************************************************************************
MODULE GETM_MOD

USE INOUTACCESS_MOD !文件路径变量
USE CAL_MOD,ONLY: MASSSOLVE, MODIDeflect, MOORLINE, INPUTSTATIC,HULLTYPE
USE ENVIRONMENT_MOD,ONLY:TTMASS,MASSMATRIX,SECTMASSMATRIX,SECTCOG,SECTPSN, ZG ,XG ,YG ,NITEM ,NBSECT,rou,G0,H,MLN1,SML,XML,YML,ZML
USE ENVIRONMENT_MOD,ONLY:XM ,YM ,ZM ,MM ,IXX ,IYY ,IZZ ,IXY,IXZ,IYZ ,ZSC ,SECTTYPE 
USE MESH_MOD
USE PRINT_MOD,ONLY:PRTSTATIC,PRTSECTLOAD

IMPLICIT NONE
PRIVATE
TYPE::MS !分段质量模型信息
   REAL(8):: MI , X1 , X2 , XGG , YR , ZGG !分段质量，起始坐标、终止坐标、重心纵向坐标、横摇惯性半径，重心垂向坐标
END TYPE MS
PUBLIC::GETMATRIX
CONTAINS

!------------------------------------------------------------------------------------------------------
SUBROUTINE GETMATRIX()
IMPLICIT NONE

REAL(8)::TEMP,RX,RY,RZ,RXZ
real(8)::Aw,Sy,hx,hy,I11,I22,I33,I13			!船舶湿面积，静矩，横、纵稳心高，惯性矩
CHARACTER(LEN=300) ::DMFINPUT
! CHARACTER(LEN=15)::MASSSOLVE
INTEGER::MASSTYPE
INTEGER::I,J
TYPE(MS),ALLOCATABLE,DIMENSION(:) :: PMASS           ! EACH MASS ITEM
INTEGER::STATMM
REAL(8)::X0,XALL
ModiDeflect="NO"
MASSMATRIX(:,:)=0.
!write(*,*)ModiDeflect
DMFINPUT=TRIM(ADJUSTL(INACCESS))//'\'//TRIM(ADJUSTL(PROJNAME))//'.DMF'

OPEN(50,FILE=DMFINPUT,STATUS='OLD')
CALL IUTMP(50)

IF(HULLTYPE=="SINGLE")THEN
    READ(50,"(A)") MASSSOLVE 
    !-----------------------读入整船质量数据----------------------------------------------
IF (trim(adjustl(MASSSOLVE))=="WHOLE")THEN  
    if(trim(adjustl(PRTSECTLOAD))=='YES')Then
    write(*,*)"如果需要计算船舶(平台)剖面载荷，请选择船舶(平台)分段质量模型或多体质量模型"
    write(*,*)
    write(*,*)"感谢您使用三维波浪载荷计算系统(WALCS)."
    write(*,*)
    stop
    end if 
    if(trim(adjustl(PRTSTATIC))=='YES')Then
    write(*,*)"如果需要计算船舶(平台)静水载荷，请选择船舶(平台)分段质量模型。"
    write(*,*)
    write(*,*)"感谢您使用三维波浪载荷计算系统(WALCS)."
    write(*,*)
    stop
    end if 

        READ(50,*) TTMASS, XG, YG, ZG !总质量 重心坐标
        
        IF (trim(adjustl(INPUTSTATIC))=="YES")THEN   !是否输入静力特性
              READ(50,*) Aw,Sy,hx,hy
	    ELSE
		      READ(50,*) temp,temp,temp,temp
	    ENDIF
    	
    !    IF (trim(adjustl(INPUTINERTIA))=="YES")THEN  !是否输入惯性矩或惯性半径
              READ(50,*) MASSTYPE    !MASStype=0..输入惯性半径 or MASStype=1..输入惯性矩
              IF(MASSTYPE==0) THEN
                    READ(50,*) RX, RY, RZ, RXZ
                    I11=TTMASS*RX**2
                    I22=TTMASS*RY**2
                    I33=TTMASS*RZ**2
                    I13=TTMASS*RXZ**2
              ELSEIF(MASSTYPE==1)THEN
                    READ(50,*) I11, I22, I33, I13
              ENDIF        
    !	ELSE
    !		  READ(50,*) temp,temp,temp,temp
    !	ENDIF
!	    ZG=ZG/cos(TrimAng)-ST    !重心垂向坐标相对参考点(吃水)的位置
	    ZG=ZG-ST	
	    MASSMATRIX(1,1)=TTMASS; MASSMATRIX(2,2)=TTMASS; MASSMATRIX(3,3)=TTMASS 
	    MASSMATRIX(4,4)=I11; MASSMATRIX(5,5)=I22; MASSMATRIX(6,6)=I33
	    MASSMATRIX(4,6)=I13; MASSMATRIX(6,4)=I13	  
	    !----------------------------check file-------------
	    WRITE(12,'(a)')'[ Global Inertia Matrix : ]'
        DO J=1,6
            WRITE(12,'(6E14.6)') (MASSMATRIX(I,J),I=1,6)
        ENDDO
        write(12,*)    

    !-----------------------读分段计算数据----------------------------------------------
    ELSEIF(TRIM(ADJUSTL(MASSSOLVE))=="SECT")THEN 
        MODIDeflect="NO"  !默认修正重力（心）
        READ(50,*) TTMASS
        READ(50,*) XG, YG, ZG
        READ(50,*) NITEM , NBSECT !质量段数目，计算剖面数目
        ALLOCATE(PMASS(1:NITEM));ALLOCATE(SECTPSN(1:3,1:NBSECT));ALLOCATE(SECTCOG(1:3,1:NBSECT))
        ALLOCATE(ZSC(1:NBSECT))
        DO I=1,NITEM
           READ(50,*) TEMP,PMASS(I)%MI, PMASS(I)%X1, PMASS(I)%X2 ,PMASS(I)%XGG, PMASS(I)%YR ,PMASS(I)%ZGG 
        ENDDO
       
       X0=PMASS(1)%X1    !重量起始点
       XALL=PMASS(NITEM)%X2   !重量终止点
       
        DO I=1,NBSECT
               READ(50,*)TEMP,(SECTPSN(J,I),J=1,3),ZSC(I)
        ENDDO
        
        IF (trim(adjustl(INPUTSTATIC))=="YES")THEN   !是否输入静力特性
              READ(50,*) Aw,Sy,hx,hy
	    ELSE
		      READ(50,*) temp,temp,temp,temp
	    ENDIF  	
    !    IF (trim(adjustl(INPUTINERTIA))=="YES")THEN  !是否输入惯性矩或惯性半径
    !          READ(50,*) MASSTYPE    !MASStype=0..输入惯性半径 or MASStype=1..输入惯性矩
    !          IF(MASSTYPE==0) THEN
    !                READ(50,*) RX, RY, RZ, RXZ
    !                I11=TTMASS*RX**2
    !                I22=TTMASS*RY**2
    !                I33=TTMASS*RZ**2
    !                I13=TTMASS*RXZ**2
    !          ELSEIF(MASSTYPE==1)THEN
    !                READ(50,*) I11, I22, I33, I13
    !          ENDIF        
    !	ELSE
    !		  READ(50,*) temp,temp,temp,temp
    !	ENDIF	           
        CLOSE(50)     
        ALLOCATE(SECTMASSMATRIX(1:6,1:6,1:NBSECT))
        !----------------------------check file--------------------------------------------------
        WRITE(12,'(a)')'[重心坐标(m),( 未修正 ) ]'
        write(12,*)
        WRITE(12,'(3F8.3)') XG, YG, ZG
        write(12,*)
        WRITE(12,'(a)')'[ 浮体总质量(ton), ( 未修正 )  ]'
        write(12,*)
        WRITE(12,'(F12.3)') TTMASS
        write(12,*)
        
        IF(TRIM(ADJUSTL(ModiDeflect))=="YES")THEN   !BY SUN,对重心纵向坐标进行修正，使其与浮心一致
	             XG=deltax
	            !----------------------------check file--------------------------------------------------
	            WRITE(12,'(a)')'[重心坐标(m),( 修正 ) ]'
                WRITE(12,*)
                WRITE(12,'(3F8.3)') XG, YG, ZG
                WRITE(12,*)
	    ENDIF   
        
        CALL SECTMASS(PMASS)   !求分段质量矩阵SECTMASSMATRIX(1:6,1:6,1:nbsect)
        CALL GLBMASS(PMASS)    !求整船质量矩阵MASSMATRIX(1:6,1:6)
        
        IF(TRIM(ADJUSTL(ModiDeflect))=="YES")THEN
	            WRITE(12,'(a)')'[ 总质量矩阵,(修正) ]'
                DO J=1,6
                    WRITE(12,'(6E14.6)') (MASSMATRIX(I,J),I=1,6)
                ENDDO
                write(12,*)
        ELSE
	            WRITE(12,'(a)')'[总质量矩阵,(未修正)  ]'
	            write(12,*)
                DO J=1,6
                    WRITE(12,'(6E14.6)') (MASSMATRIX(I,J),I=1,6)
                ENDDO
                write(12,*)             
        ENDIF
        
	    ZG=ZG-ST    !重心垂向坐标相对参考点(吃水)的位置
!	    ZG=ZG/cos(TrimAng)-ST
	    ZSC(:)=ZSC(:)-ST
	    sectcog(1,:) = sectcog(1,:)-XG
	    sectpsn(1,:) = sectpsn(1,:)-XG
	    sectpsn(3,:) = sectpsn(3,:)-ST
	    
	    IF(TRIM(ADJUSTL(PRTSTATIC))=="YES")THEN
	           CALL getHYDROSTATIC(X0,XALL)
	    ENDIF   
	    
    END IF
!-----------------------读入多体质量信息文件----------------------------------------------       
ELSEIF(TRIM(ADJUSTL(HULLTYPE))=="MULTI")	    THEN
    if(trim(adjustl(PRTSTATIC))=='YES')Then
    write(*,*)"如果需要计算船舶(平台)静水载荷，请选择船舶(平台)分段质量模型。"
    write(*,*)
    write(*,*)"感谢您使用三维波浪载荷计算系统(WALCS)."
    write(*,*)
    stop
    end if 
        READ(50,"(A)") MASSSOLVE
        READ(50,*)   NITEM , NBSECT
       
        ALLOCATE(XM(1:NITEM),YM(1:NITEM),ZM(1:NITEM),MM(1:NITEM),IXX(1:NITEM),IYY(1:NITEM),IZZ(1:NITEM),IXY(1:NITEM),IXZ(1:NITEM),IYZ(1:NITEM))
        ALLOCATE(SECTPSN(1:3,1:NBSECT),ZSC(1:NBSECT))
        ALLOCATE(SECTTYPE(1:NBSECT))
        ALLOCATE(SECTMASSMATRIX(1:6,1:6,1:NBSECT))
         
        DO I=1,   NITEM
            READ(50,*) TEMP,  XM(I), YM(I),  ZM(I),  MM(I), IXX(I),IYY(I),IZZ(I),IXY(I),IXZ(I),IYZ(I)
        ENDDO
    	
        DO  I=1,  NBSECT
            READ(50,*) TEMP,    SECTPSN(1,I),    SECTPSN(2,I), SECTPSN(3,I),  ZSC(I), SECTTYPE(I)   
        END   DO
       
        CALL MULTIMASS  ()
        DEALLOCATE(XM,YM,ZM,MM,IXX,IYY,IZZ,IXY,IXZ,IYZ)
        
        IF (trim(adjustl(INPUTSTATIC))=="YES")THEN   !是否输入静力特性
              READ(50,*) Aw,Sy,hx,hy
        ELSE
		      READ(50,*) temp,temp,temp,temp
        ENDIF
        
        ZG=ZG-ST
        ZSC(:)=ZSC(:)-ST
        CALL  MASSMATRIXPUTOUT()
END IF

if(trim(adjustl(HULLTYPE))=="SINGLE")then
   CALL RESTOMATRIX(Aw,Sy,hx,hy)    !求整体恢复力矩阵HRM(1:6,1:6)
ELSEif(trim(adjustl(HULLTYPE))=="MULTI")then
   CALL RESTOMATRIX_SECT(Aw,Sy,hx,hy)  !计算整体及分块恢复力矩阵和湿表面积
endif
Ta=ST*cos(TrimAng)-Xg*tan(TrimAng)
END SUBROUTINE GETMATRIX

!----------------------------------------------------------------------------------------------------
SUBROUTINE SECTMASS(PMASS)
IMPLICIT NONE

REAL(8)::ttDeflect
REAL(8),dimension(1:nitem)::Deflect
INTEGER :: I , J 
TYPE(MS),INTENT(INOUT),DIMENSION(1:NITEM) :: PMASS           ! EACH MASS ITEM
REAL(8):: M , MIX , MIZ , I44 , I55 , I66 , I45 , I46 , I56
REAL(8) :: XGC , YGC , ZGC
REAL(8) :: T1 ,TEMP
REAL(8),DIMENSION(1:3):: COG ,INCOG           ! CENTRE GRAVITY OF THE BODY
INCOG (1)=XG; INCOG(2)=YG; INCOG(3)=ZG

COG = 0 ; PMASS(:)%MI = PMASS(:)%MI * TTMASS/SUM(PMASS(:)%MI)

Deflect=0.
!-------------重力与浮力差加权平均到各质量段上-----------------------------------
IF(TRIM(ADJUSTL(ModiDeflect))=="YES")THEN
     ttDeflect=rou*vol-SUM(PMASS(:)%MI)   !by sun
     do i=1,nitem
        Deflect(i)=(PMASS(i)%X2-PMASS(i)%X1)/(PMASS(Nitem)%X2-PMASS(1)%X1)*ttDeflect
        PMASS(I)%MI=PMASS(I)%MI+Deflect(i)
     end do
end if
!---------------------------------------------------------------------------------
!DO I = 1 , NITEM
!   COG(1) = COG(1) + PMASS(I)%MI*PMASS(I)%XGG
!   COG(3) = COG(3) + PMASS(I)%MI*PMASS(I)%ZGG
!   COG(2) = COG(2) + PMASS(I)%MI
!END DO
!COG(:) = COG(:) / COG(2) ; COG(2) = 0
!
!PMASS(:)%X1 = PMASS(:)%X1 - ( COG(1)-INCOG(1) )
!PMASS(:)%X2 = PMASS(:)%X2 - ( COG(1)-INCOG(1) )
!PMASS(:)%XGG = PMASS(:)%XGG - ( COG(1)-INCOG(1) )
!PMASS(:)%ZGG = PMASS(:)%ZGG - ( COG(3)-INCOG(3) )
!COG = INCOG

SECTMASSMATRIX = 0. ; SECTCOG = 0.
DO I = 1 , NBSECT
   M = 0. ; MIX = 0 ; MIZ = 0 ; I44 = 0 ;  I55 = 0 ; I66 = 0 ; I45 = 0 ; I46 = 0 ; I56 = 0
   XGC = 0 ; YGC = 0 ; ZGC = 0
   DO J = 1 , NITEM
      IF(PMASS(J)%X2<=SECTPSN(1,I)) THEN
	    M = M + PMASS(J)%MI
		MIX = MIX + PMASS(J)%MI * PMASS(J)%XGG
        MIZ = MIZ + PMASS(J)%MI * PMASS(J)%ZGG
	  ELSE IF (PMASS(J)%X2>SECTPSN(1,I).AND.PMASS(J)%X1<SECTPSN(1,I)) THEN
	    T1 = PMASS(J)%MI/(PMASS(J)%X2-PMASS(J)%X1)*(SECTPSN(1,I)-PMASS(J)%X1)
		M = M + T1
		MIX = MIX + T1 * ((PMASS(J)%X1+SECTPSN(1,I))*0.5)
		MIZ = MIZ + T1 * PMASS(J)%ZGG
	  END IF
   END DO
   if(M<0.0000001)then
   write(*,"(a,i3,a)")"剖面",i," 设置无效，请重新设置。"
   stop
   endif
   SECTCOG(1,I) = MIX/M ; SECTCOG(2,I) = 0 ; SECTCOG(3,I) = MIZ/M
   XGC = SECTCOG(1,I) - INCOG(1)
   YGC = SECTPSN(2,I) - INCOG(2)
   ZGC = SECTCOG(3,I) - INCOG(3) 

   DO J = 1 , NITEM
      IF(PMASS(J)%X2<=SECTPSN(1,I)) THEN

		I44 = I44 + PMASS(J)%MI*PMASS(J)%YR**2 + PMASS(J)%MI*(INCOG(3)-PMASS(J)%ZGG)**2
		I55 = I55 + PMASS(J)%MI*((PMASS(J)%X2-PMASS(J)%X1)/4)**2 + & 
		       PMASS(J)%MI*((INCOG(1)-PMASS(J)%XGG)**2 + (INCOG(3)-PMASS(J)%ZGG)**2)        
		I66 = I66 + PMASS(J)%MI*((PMASS(J)%X2-PMASS(J)%X1)/4)**2 + &
		       PMASS(J)%MI*(INCOG(1)-PMASS(J)%XGG)**2  + PMASS(J)%MI*(INCOG(2)-PMASS(J)%YR)**2
		I46 = I46 + PMASS(J)%MI*(-INCOG(1)+PMASS(J)%XGG)*(-INCOG(3)+PMASS(J)%ZGG) 
		
	  ELSE IF (PMASS(J)%X2>=SECTPSN(1,I).AND.PMASS(J)%X1<=SECTPSN(1,I)) THEN
	    T1 = PMASS(J)%MI/(PMASS(J)%X2-PMASS(J)%X1)*(SECTPSN(1,I)-PMASS(J)%X1)

		I44 = I44 + T1*PMASS(J)%YR**2 + T1*(INCOG(3)-PMASS(J)%ZGG)**2 
		I55 = I55 + T1*((SECTPSN(1,I)-PMASS(J)%X1)/4)**2 + & 
		      T1*((INCOG(1)-(PMASS(J)%X1+SECTPSN(1,I))*0.5)**2+(INCOG(3)-PMASS(J)%ZGG)**2) 
		     
		I66 = I66 + T1*( (SECTPSN(1,I)-PMASS(J)%X1)/4 )**2 + &
		       T1*(INCOG(1)-(PMASS(J)%X1+SECTPSN(1,I))*0.5)**2 + T1*(INCOG(2)-PMASS(J)%YR)**2
		      
		I46 = I46 + T1*(-INCOG(1)+(PMASS(J)%X1+SECTPSN(1,I))*0.5)*(-INCOG(3)+PMASS(J)%ZGG) 
		
	  END IF
   END DO
   I46 = -I46
   SECTMASSMATRIX(1,1,I) = M       ; SECTMASSMATRIX(2,2,I) = M       ; SECTMASSMATRIX(3,3,I) = M
   SECTMASSMATRIX(4,4,I) = I44     ; SECTMASSMATRIX(5,5,I) = I55     ; SECTMASSMATRIX(6,6,I) = I66

   SECTMASSMATRIX(1,5,I) = M*ZGC   ; SECTMASSMATRIX(1,6,I) = -M*YGC  ; SECTMASSMATRIX(2,4,I) = -M*ZGC
   SECTMASSMATRIX(2,6,I) = M*XGC   ; SECTMASSMATRIX(3,4,I) = M*YGC   ; SECTMASSMATRIX(3,5,I) = -M*XGC
   SECTMASSMATRIX(4,5,I) = I45     ; SECTMASSMATRIX(4,6,I) = I46     ; SECTMASSMATRIX(5,6,I) = I56
   ! GET SYMMETRY
   SECTMASSMATRIX(5,1,I) = SECTMASSMATRIX(1,5,I)  ; SECTMASSMATRIX(6,1,I) = SECTMASSMATRIX(1,6,I)  ; SECTMASSMATRIX(4,2,I) = SECTMASSMATRIX(2,4,I)
   SECTMASSMATRIX(6,2,I) = SECTMASSMATRIX(2,6,I)  ; SECTMASSMATRIX(4,3,I) = SECTMASSMATRIX(3,4,I)  ; SECTMASSMATRIX(5,3,I) = SECTMASSMATRIX(3,5,I)
   SECTMASSMATRIX(5,4,I) = SECTMASSMATRIX(4,5,I)  ; SECTMASSMATRIX(6,4,I) = SECTMASSMATRIX(4,6,I)  ; SECTMASSMATRIX(6,5,I) = SECTMASSMATRIX(5,6,I)
END DO

END SUBROUTINE SECTMASS

!**********************************************************************************************
SUBROUTINE GLBMASS(PMASS)
IMPLICIT NONE

REAL(8)::ttDeflect
REAL(8),dimension(1:nitem)::Deflect
INTEGER :: I , J 
TYPE(MS),INTENT(INOUT),DIMENSION(1:NITEM) :: PMASS  
REAL(8):: M , I44 , I55 , I66 , I45 , I46 , I56
REAL(8) :: XGC , YGC , ZGC
REAL(8),DIMENSION(1:3):: COG ,INCOG           ! CENTRE GRAVITY OF THE BODY

INCOG (1)=XG; INCOG(2)=YG; INCOG(3)=ZG
COG = 0 ; PMASS(:)%MI = PMASS(:)%MI * TTMASS/SUM(PMASS(:)%MI)

Deflect=0.
!-------------重力与浮力差加权平均到各质量段上-----------------------------------
IF(TRIM(ADJUSTL(ModiDeflect))=="YES")THEN
     ttDeflect=rou*vol-ttmass   !by sun
     do i=1,nitem
        Deflect(i)=(PMASS(i)%X2-PMASS(i)%X1)/(PMASS(Nitem)%X2-PMASS(1)%X1)*ttDeflect
        PMASS(I)%MI=PMASS(I)%MI+Deflect(i)
     end do
end if
!---------------------------------------------------------------------------------
M = 0. ; I44 = 0 ;  I55 = 0 ; I66 = 0 ; I45 = 0 ; I46 = 0 ; I56 = 0
XGC = 0 ; YGC = 0 ; ZGC = 0 

DO I = 1 , NITEM
   COG(1) = COG(1) + PMASS(I)%MI*PMASS(I)%XGG
   COG(3) = COG(3) + PMASS(I)%MI*PMASS(I)%ZGG
   COG(2) = COG(2) + PMASS(I)%MI
END DO

M=COG(2)
COG(:) = COG(:) / COG(2) ; COG(2) = 0

PMASS(:)%X1 = PMASS(:)%X1 - ( COG(1)-INCOG(1) )
PMASS(:)%X2 = PMASS(:)%X2 - ( COG(1)-INCOG(1) )
PMASS(:)%XGG = PMASS(:)%XGG - ( COG(1)-INCOG(1) )
PMASS(:)%ZGG = PMASS(:)%ZGG - ( COG(3)-INCOG(3) )

COG(:)=INCOG(:)

XGC = COG(1) - INCOG(1)
YGC = COG(2) - INCOG(2)
ZGC = COG(3) - INCOG(3) 

DO J = 1 , NITEM
      IF(PMASS(J)%X2<=1000) THEN

		I44 = I44 + PMASS(J)%MI*PMASS(J)%YR**2 + PMASS(J)%MI*(INCOG(3)-PMASS(J)%ZGG)**2
		I55 = I55 + PMASS(J)%MI*((PMASS(J)%X2-PMASS(J)%X1)/4)**2 + & 
		       PMASS(J)%MI*((INCOG(1)-PMASS(J)%XGG)**2 + (INCOG(3)-PMASS(J)%ZGG)**2)        
		I66 = I66 + PMASS(J)%MI*((PMASS(J)%X2-PMASS(J)%X1)/4)**2 + &
		       PMASS(J)%MI*(INCOG(1)-PMASS(J)%XGG)**2  + PMASS(J)%MI*(INCOG(2)-PMASS(J)%YR)**2
		I46 = I46 + PMASS(J)%MI*(-INCOG(1)+PMASS(J)%XGG)*(-INCOG(3)+PMASS(J)%ZGG) 
		
	  END IF
END DO
   I46 = -I46
   MASSMATRIX(1,1) = M       ; MASSMATRIX(2,2) = M       ; MASSMATRIX(3,3) = M
   MASSMATRIX(4,4) = I44     ; MASSMATRIX(5,5) = I55     ; MASSMATRIX(6,6) = I66

   MASSMATRIX(1,5) = M*ZGC   ; MASSMATRIX(1,6) = -M*YGC  ; MASSMATRIX(2,4) = -M*ZGC
   MASSMATRIX(2,6) = M*XGC   ; MASSMATRIX(3,4) = M*YGC   ; MASSMATRIX(3,5) = -M*XGC
   MASSMATRIX(4,5) = I45     ; MASSMATRIX(4,6) = I46     ; MASSMATRIX(5,6) = I56
   ! GET SYMMETRY
   MASSMATRIX(5,1) = MASSMATRIX(1,5)  ; MASSMATRIX(6,1) = MASSMATRIX(1,6)  ; MASSMATRIX(4,2) = MASSMATRIX(2,4)
   MASSMATRIX(6,2) = MASSMATRIX(2,6)  ; MASSMATRIX(4,3) = MASSMATRIX(3,4)  ; MASSMATRIX(5,3) = MASSMATRIX(3,5)
   MASSMATRIX(5,4) = MASSMATRIX(4,5)  ; MASSMATRIX(6,4) = MASSMATRIX(4,6)  ; MASSMATRIX(6,5) = MASSMATRIX(5,6)

END SUBROUTINE GLBMASS

!----------------------------------------------------------------------
!求多体各质量段的分段质量矩阵（相对于全船重心）
!----------------------------------------------------------------------
SUBROUTINE MULTIMASS  ()
IMPLICIT   NONE
!REAL(8),   DIMENSION(1:NITEM)      ::  XM,YM,ZM,MM,IXX,IYY,IZZ,IXY,IXZ,IYZ
!REAL(8),   DIMENSION(1:3,1:NBSECT) ::  SECTPSN
!INTEGER,   DIMENSION(1:NBSECT)     ::  SECTTYPE
REAL(8) :: M
INTEGER  ::I,J,K
!--------------------求出质心位置 ----------------------------------------------
M=SUM(MM(:))
XG=SUM(MM(:)*XM(:))/M
YG=SUM(MM(:)*YM(:))/M
ZG=SUM(MM(:)*ZM(:))/M
!--------------------进行坐标转换---------------------------
XM(:)=XM(:)-XG
YM(:)=YM(:)-YG
ZM(:)=ZM(:)-ZG

SECTPSN(1,:)=SECTPSN(1,:)-XG
SECTPSN(2,:)=SECTPSN(2,:)-YG
SECTPSN(3,:)=SECTPSN(3,:)-ST  
!-------------------求分段质量矩阵--------------------------
SECTMASSMATRIX=0.0

DO  I=1,  NBSECT
      DO  J=1, NITEM    
            IF(SECTTYPE(I)==1.AND.XM(J)<=SECTPSN(1,I))THEN              
                    SECTMASSMATRIX(1,1,I) = SECTMASSMATRIX(1,1,I)+MM(J)
                    SECTMASSMATRIX(2,2,I) = SECTMASSMATRIX(1,1,I)
                    SECTMASSMATRIX(3,3,I) = SECTMASSMATRIX(1,1,I)
                    
                    SECTMASSMATRIX(4,4,I) = SECTMASSMATRIX(4,4,I)+Ixx(J)+MM(J)*YM(J)*YM(J)+MM(J)*ZM(J)*ZM(J)
                    SECTMASSMATRIX(5,5,I) = SECTMASSMATRIX(5,5,I)+Iyy(J)+MM(J)*XM(J)*XM(J)+MM(J)*ZM(J)*ZM(J)
                    SECTMASSMATRIX(6,6,I) = SECTMASSMATRIX(6,6,I)+Izz(J)+MM(J)*XM(J)*XM(J)+MM(J)*YM(J)*YM(J)
                    
                    SECTMASSMATRIX(1,5,I) = SECTMASSMATRIX(1,5,I)+MM(J)*ZM(J)
                    SECTMASSMATRIX(5,1,I) = SECTMASSMATRIX(1,5,I)
                    SECTMASSMATRIX(2,4,I) =-SECTMASSMATRIX(1,5,I)
                    SECTMASSMATRIX(4,2,I) =-SECTMASSMATRIX(1,5,I)
                    
                    SECTMASSMATRIX(1,6,I) = SECTMASSMATRIX(1,6,I)-MM(J)*YM(J)
                    SECTMASSMATRIX(6,1,I) = SECTMASSMATRIX(1,6,I)
                    SECTMASSMATRIX(3,4,I) =-SECTMASSMATRIX(1,6,I)
                    SECTMASSMATRIX(4,3,I) =-SECTMASSMATRIX(1,6,I)
                    
                    SECTMASSMATRIX(2,6,I) = SECTMASSMATRIX(2,6,I)+MM(J)*XM(J)
                    SECTMASSMATRIX(6,2,I) = SECTMASSMATRIX(2,6,I)
                    SECTMASSMATRIX(3,5,I) =-SECTMASSMATRIX(2,6,I)
                    SECTMASSMATRIX(5,3,I) = SECTMASSMATRIX(3,5,I)
                    
                    SECTMASSMATRIX(4,5,I) = SECTMASSMATRIX(4,5,I)-MM(J)*XM(J)*YM(J)-Ixy(J)
                    SECTMASSMATRIX(5,4,I) = SECTMASSMATRIX(4,5,I)
                    
                    SECTMASSMATRIX(4,6,I) = SECTMASSMATRIX(4,6,I)-MM(J)*XM(J)*ZM(J)-Ixz(J)
                    SECTMASSMATRIX(6,4,I) = SECTMASSMATRIX(4,6,I)
                    
                    SECTMASSMATRIX(5,6,I) = SECTMASSMATRIX(5,6,I)-Iyz(J)-MM(J)*YM(J)*ZM(J)
                    SECTMASSMATRIX(6,5,I) = SECTMASSMATRIX(5,6,I)
               ELSEIF(SECTTYPE(I)==2.AND.YM(J)<=SECTPSN(2,I))THEN
                
                    SECTMASSMATRIX(1,1,I) = SECTMASSMATRIX(1,1,I)+MM(J)
                    SECTMASSMATRIX(2,2,I) = SECTMASSMATRIX(1,1,I)
                    SECTMASSMATRIX(3,3,I) = SECTMASSMATRIX(1,1,I)
                    
                    SECTMASSMATRIX(4,4,I) = SECTMASSMATRIX(4,4,I)+Ixx(J)+MM(J)*YM(J)*YM(J)+MM(J)*ZM(J)*ZM(J)
                    SECTMASSMATRIX(5,5,I) = SECTMASSMATRIX(5,5,I)+Iyy(J)+MM(J)*XM(J)*XM(J)+MM(J)*ZM(J)*ZM(J)
                    SECTMASSMATRIX(6,6,I) = SECTMASSMATRIX(6,6,I)+Izz(J)+MM(J)*XM(J)*XM(J)+MM(J)*YM(J)*YM(J)
                    
                    SECTMASSMATRIX(1,5,I) = SECTMASSMATRIX(1,5,I)+MM(J)*ZM(J)
                    SECTMASSMATRIX(5,1,I) = SECTMASSMATRIX(1,5,I)
                    SECTMASSMATRIX(2,4,I) =-SECTMASSMATRIX(1,5,I)
                    SECTMASSMATRIX(4,2,I) =-SECTMASSMATRIX(1,5,I)
                    
                    SECTMASSMATRIX(1,6,I) = SECTMASSMATRIX(1,6,I)-MM(J)*YM(J)
                    SECTMASSMATRIX(6,1,I) = SECTMASSMATRIX(1,6,I)
                    SECTMASSMATRIX(3,4,I) =-SECTMASSMATRIX(1,6,I)
                    SECTMASSMATRIX(4,3,I) =-SECTMASSMATRIX(1,6,I)
                    
                    SECTMASSMATRIX(2,6,I) = SECTMASSMATRIX(2,6,I)+MM(J)*XM(J)
                    SECTMASSMATRIX(6,2,I) = SECTMASSMATRIX(2,6,I)
                    SECTMASSMATRIX(3,5,I) =-SECTMASSMATRIX(2,6,I)
                    SECTMASSMATRIX(5,3,I) = SECTMASSMATRIX(3,5,I)
                    
                    SECTMASSMATRIX(4,5,I) = SECTMASSMATRIX(4,5,I)-MM(J)*XM(J)*YM(J)-Ixy(J)
                    SECTMASSMATRIX(5,4,I) = SECTMASSMATRIX(4,5,I)
                    
                    SECTMASSMATRIX(4,6,I) = SECTMASSMATRIX(4,6,I)-MM(J)*XM(J)*ZM(J)-Ixz(J)
                    SECTMASSMATRIX(6,4,I) = SECTMASSMATRIX(4,6,I)
                    
                    SECTMASSMATRIX(5,6,I) = SECTMASSMATRIX(5,6,I)-Iyz(J)-MM(J)*YM(J)*ZM(J)
                    SECTMASSMATRIX(6,5,I) = SECTMASSMATRIX(5,6,I) 
                  ENDIF              
                                                                              
        END   DO
    
END   DO
!-------------------------求整体质量矩阵------------------------------

DO  J=1, NITEM

    MASSMATRIX(1,1) = MASSMATRIX(1,1)+MM(J)
    MASSMATRIX(2,2) = MASSMATRIX(1,1)
    MASSMATRIX(3,3) = MASSMATRIX(1,1)
                
    MASSMATRIX(4,4) = MASSMATRIX(4,4)+Ixx(J)+MM(J)*YM(J)*YM(J)+MM(J)*ZM(J)*ZM(J)
    MASSMATRIX(5,5) = MASSMATRIX(5,5)+Iyy(J)+MM(J)*XM(J)*XM(J)+MM(J)*ZM(J)*ZM(J)
    MASSMATRIX(6,6) = MASSMATRIX(6,6)+Izz(J)+MM(J)*XM(J)*XM(J)+MM(J)*YM(J)*YM(J)
                
    MASSMATRIX(1,5) = MASSMATRIX(1,5)+MM(J)*ZM(J)
    MASSMATRIX(5,1) = MASSMATRIX(1,5)
    MASSMATRIX(2,4) =-MASSMATRIX(1,5)
    MASSMATRIX(4,2) =-MASSMATRIX(1,5)
                
    MASSMATRIX(1,6) = MASSMATRIX(1,6)-MM(J)*YM(J)
    MASSMATRIX(6,1) = MASSMATRIX(1,6)
    MASSMATRIX(3,4) =-MASSMATRIX(1,6)
    MASSMATRIX(4,3) =-MASSMATRIX(1,6)
                
    MASSMATRIX(2,6) = MASSMATRIX(2,6)+MM(J)*XM(J)
    MASSMATRIX(6,2) = MASSMATRIX(2,6)
    MASSMATRIX(3,5) =-MASSMATRIX(2,6)
    MASSMATRIX(5,3) = MASSMATRIX(3,5)
                
    MASSMATRIX(4,5) = MASSMATRIX(4,5)-MM(J)*XM(J)*YM(J)-Ixy(J)
    MASSMATRIX(5,4) = MASSMATRIX(4,5)
                
    MASSMATRIX(4,6) = MASSMATRIX(4,6)-MM(J)*XM(J)*ZM(J)-Ixz(J)
    MASSMATRIX(6,4) = MASSMATRIX(4,6)
                
    MASSMATRIX(5,6) = MASSMATRIX(5,6)-Iyz(J)-MM(J)*YM(J)*ZM(J)
    MASSMATRIX(6,5) = MASSMATRIX(5,6)   
               
ENDDO                            

ENDSUBROUTINE MULTIMASS

!----------------------------------------------------------------------
!hydrostatic restoring matrix!求整体恢复力矩阵HRM(1:6,1:6)
!----------------------------------------------------------------------
SUBROUTINE RESTOMATRIX(Aw,Sy,hx,hy)
implicit none

real(8),INTENT(INOUT)::Aw,Sy,hx,hy
REAL(8)::XF
CHARACTER(LEN=300)::MOLINPUT
real(8)::mlT,mlH 	!单根锚链预张力,锚链水深

real(8)::mlXG,mlZG						
integer::mlN2						!单根锚链分段数
real(8),allocatable,dimension(:)::	REF,OMG,EL1,W1
real(8),dimension(1:3)::tempcoor
real(8),dimension(1:6,1:6)::	CML				!锚链刚度系数
integer::i

HRM=0.0
YB=0.0


IF (trim(adjustl(INPUTSTATIC))=="NO")THEN
        if(Factor==1)then
        do i=1,NB
		        HRM(3,3)=HRM(3,3)+e(3,3,i)*ea(i)
		        HRM(3,4)=HRM(3,4)+xav(2,i)*e(3,3,i)*ea(i)
		        HRM(3,5)=HRM(3,5)-xav(1,i)*e(3,3,i)*ea(i)
		        HRM(4,4)=HRM(4,4)+xav(2,i)*xav(2,i)*e(3,3,i)*ea(i)
		        HRM(4,5)=HRM(4,5)-xav(1,i)*xav(2,i)*e(3,3,i)*ea(i)
		        HRM(5,5) =HRM(5,5)+xav(1,i)*xav(1,i)*e(3,3,i)*ea(i)
        enddo
        else if(Factor==2)then
        do i=1,NB
		        HRM(3,3)=HRM(3,3)+2*e(3,3,i)*ea(i)
		        HRM(3,4)=0
		        HRM(3,5)=HRM(3,5)-2*xav(1,i)*e(3,3,i)*ea(i)
		        HRM(4,4)=HRM(4,4)+2*xav(2,i)*xav(2,i)*e(3,3,i)*ea(i)
		        HRM(4,5)=0
		        HRM(5,5)=HRM(5,5)+2*xav(1,i)*xav(1,i)*e(3,3,i)*ea(i)
        enddo
                else if(Factor==4)then
        do i=1,NB
		        HRM(3,3)=HRM(3,3)+4*e(3,3,i)*ea(i)
		        HRM(3,4)=0 
		        HRM(3,5)=0
		        HRM(4,4)=HRM(4,4)+4*xav(2,i)*xav(2,i)*e(3,3,i)*ea(i)
		        HRM(4,5)=0
		        HRM(5,5)=HRM(5,5)+4*xav(1,i)*xav(1,i)*e(3,3,i)*ea(i)
        enddo
        end if
        HRM(4,4)=HRM(4,4)+vol*(ZB-ZG)
        HRM(5,5)=HRM(5,5)+vol*(ZB-ZG)
!!!!!!!!!! HRM(4,6)=-vol*XB;HRM(5,6)=-vol*YB !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !-新增---------------------
        HRM(5,3)=HRM(3,5)
        HRM(5,4)=HRM(4,5)
        HRM(4,3)=HRM(3,4)
        Aw=HRM(3,3)
	    Sy=-HRM(3,5)
	    XF=Sy/Aw
	    hx=HRM(4,4)/vol
	    hy=HRM(5,5)/vol
	    HRM=HRM*rou*g0
 ELSEIF(trim(adjustl(INPUTSTATIC))=="YES")THEN
 	    HRM(3,3)=rou*g0*Aw
	    HRM(3,5)=-rou*g0*Sy
	    HRM(5,3)=HRM(3,5)
	    HRM(4,4)=rou*g0*vol*hx
	    HRM(5,5)=rou*g0*vol*hy
 ENDIF
!----------------------------check file-------------------------------
write(12,'(a)')"#------------------------------------------------------------ "
write(12,*)
write(12,'(a)')"[静水恢复力矩阵]"
write(12,*)
write(12,'(a,e12.4)')"AW= ",AW
write(12,'(a,e12.4)')"Sy= ",Sy
write(12,'(a,f12.4)')"hx= ",hx
write(12,'(a,f12.4)')"hy= ",hy
write(12,'(a,e12.4)')"C33= ",HRM(3,3)
write(12,'(a,f12.4)')"C34=C43= ",HRM(3,4)
write(12,'(a,e12.4)')"C35=C53= ",HRM(3,5) 
write(12,'(a,e12.4)')"C44= ",HRM(4,4)
write(12,'(a,f12.4)')"C45=C54= ",HRM(4,5) 
!!!!!!write(12,'(a,f12.4)')"C46= ",HRM(4,6) !!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(12,'(a,e12.4)')"C55= ",HRM(5,5)
!!!!!!!!write(12,'(a,f12.4)')"C56= ",HRM(5,6) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(12,*)
!----------------------------------------------------------------------------------------
if((trim(adjustl(MoorLine))=="YES"))then       !读入锚链数据
    MOLINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.MOL'
    OPEN (27,FILE=MOLINPUT,STATUS="OLD")
	call IUTMP(27)
	READ(27,*)mlN1,mlN2,mlT
	mlH=H                                    !工作水深
	allocate(REF(1:mlN1),OMG(1:mlN2),EL1(1:mlN2),W1(1:mlN2))
    mlxg=xg
	mlzg=zg
	read(27,*)(tempcoor(i),i=1,3)
	allocate(xMl(1:mlN1),yMl(1:mlN1),zMl(1:mlN1))
	do	i=1,	mlN1
		read(27,*)xMl(i),yMl(i),zMl(i)
	enddo
	xMl(:)=xMl(:)+tempcoor(1)
	yMl(:)=yMl(:)+tempcoor(2)
	zMl(:)=zMl(:)+tempcoor(3)
	READ(27,*) (REF(I),I=1,mlN1)
	READ(27,*) (OMG(I),I=1,mlN2)
	READ(27,*) (W1(I),I=1,mlN2-1)
	READ(27,*)  (EL1(I),I=1,mlN2)

	CML=0.0
	allocate(SML(1:6,1:6,1:mlN1))
	SML=0.0

	do	i=1,	mlN1
		call ML(1,mlN2,mlT,mlH,mlXG,mlZG,XML(i),ZML(i),REF,OMG,EL1,W1,SML(:,:,i))
		CML(:,:)=CML(:,:)+SML(:,:,i)
	end	do
	write(12,*)
	write(12,'(a)')"The stiffness coefficients of the moor lines:"
	write(12,*)
	write(12,'(a,e12.4)')"CML33= ",CML(3,3)
	write(12,'(a,e12.4)')"CML35= ",CML(3,5)
	write(12,'(a,e12.4)')"CML53= ",CML(5,3)
	write(12,'(a,e12.4)')"CML55= ",CML(5,5)
	write(12,*)

	HRM=HRM+CML
  
    DEALLOCATE(REF,OMG,EL1,W1)
ENDIF

END SUBROUTINE RESTOMATRIX
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!计算局部恢复力矩阵
!-------------------------------------------------------------------------------
SUBROUTINE RESTOMATRIX_SECT(Aw,Sy,hx,hy)  
IMPLICIT NONE

real(8),INTENT(INOUT)::Aw,Sy,hx,hy
real(8),ALLOCATABLE,dimension(:,:,:)::eu !各节点的位移及其导数
real(8)::TEMPHRM(1:6,1:6),hrm1(1:6,1:6),ES(1:nb)

REAL(8)::XF
INTEGER::I,J,K,L,II,count
CHARACTER(LEN=100)::MOLINPUT
real(8)::mlT,mlH 	!单根锚链预张力,锚链水深
real(8)::mlXG,mlZG						
integer::mlN2						!单根锚链分段数
real(8),allocatable,dimension(:)::	REF,OMG,EL1,W1
real(8),allocatable,dimension(:,:)::x00
real(8),dimension(1:3)::tempcoor
real(8),dimension(1:6,1:6)::	CML				!锚链刚度系数
real(8),allocatable,dimension(:,:,:)::	ee
ALLOCATE(EU(1:3,1:6,1:NB),x00(1:3,1:NB))
ALLOCATE(SECTHRM(1:6,1:6,1:NBSECT),SECTEA(1:NBSECT))
Allocate(ee(1:6,1:3,1:NB))
SECTHRM=0.0
SECTEA=0.0
wetea=0.0
HRM=0.0
do count=1,TarNum,Space
select case(count)
       case(1)
              x00(1,:)=xav(1,:)
              x00(2,:)=xav(2,:)
              ee(1,3,:)=e(1,3,:)
              ee(2,3,:)=e(2,3,:)
              ee(3,3,:)=e(3,3,:)
        case(2)
              x00(1,:)=-xav(1,:)
              x00(2,:)=xav(2,:)   
              ee(1,3,:)=-e(1,3,:)
              ee(2,3,:)=e(2,3,:)
              ee(3,3,:)=e(3,3,:)    
       case(3)
             x00(1,:)=-xav(1,:)
             x00(2,:)=-xav(2,:)
              ee(1,3,:)=-e(1,3,:)
              ee(2,3,:)=-e(2,3,:)
              ee(3,3,:)=e(3,3,:)
        case(4)
            x00(1,:)=xav(1,:)
            x00(2,:)=-xav(2,:)
              ee(1,3,:)=e(1,3,:)
              ee(2,3,:)=-e(2,3,:)
              ee(3,3,:)=e(3,3,:)
 end select

eu=0.0
eu(1,1,1:NB)=1.0;eu(2,2,1:NB)=1.0;eu(3,3,1:NB)=1.0
eu(2,4,1:NB)=-(xav(3,1:NB)-zg);eu(3,4,1:NB)=x00(2,1:NB)
eu(1,5,1:NB)=xav(3,1:NB)-zg;eu(3,5,1:NB)=-x00(1,1:NB)
eu(1,6,1:NB)=-x00(2,1:NB);eu(2,6,1:NB)=x00(1,1:NB)
ES=EA
do k=1,NB	
    do i=1,6
	    do j=1,6
			TEMPHRM(i,j)=rou*g0*dot_product(ee(1:3,3,k),eu(1:3,i,k))*eu(3,j,k)*ea(k)
		enddo
	enddo
	DO II=1,NBSECT
        IF(SECTTYPE(II)==1)THEN
            L=1
        ELSEIF(SECTTYPE(II)==2)THEN
            L=2
        ENDIF
        IF(X00(L,K)<=SECTPSN(L,II))THEN
            SECTHRM(:,:,II)=SECTHRM(:,:,II)+TEMPHRM(:,:)
            SECTEA(II)=SECTEA(II)+EA(K)  !计算局部湿面积
        ENDIF
    ENDDO
    wetea=wetea+ea(k)
    HRM(:,:)=HRM(:,:)+TEMPHRM(:,:)
    HRM1=HRM
enddo
end do

DO II=1,NBSECT
    SECTHRM(2,4,II)=-(-SECTMASSMATRIX(1,1,II))*G0+SECTHRM(2,4,II)
    SECTHRM(4,4,II)=-(SECTMASSMATRIX(1,5,II))*G0+SECTHRM(4,4,II)
    SECTHRM(6,4,II)=-(-SECTMASSMATRIX(2,6,II))*G0+SECTHRM(6,4,II) 
    SECTHRM(1,5,II)=-(SECTMASSMATRIX(1,1,II))*G0+SECTHRM(1,5,II)
    SECTHRM(5,5,II)=-(-SECTMASSMATRIX(2,4,II))*G0+SECTHRM(5,5,II)
    SECTHRM(6,5,II)=-(SECTMASSMATRIX(1,6,II))*G0+SECTHRM(6,5,II)          
ENDDO
HRM(2,4)=-(-MASSMATRIX(1,1))*G0+HRM(2,4)
HRM(4,4)=-(MASSMATRIX(1,5))*G0+HRM(4,4)
HRM(6,4)=-(-MASSMATRIX(2,6))*G0+HRM(6,4) 
HRM(1,5)=-(MASSMATRIX(1,1))*G0+HRM(1,5)
HRM(5,5)=-(-MASSMATRIX(2,4))*G0+HRM(5,5)
HRM(6,5)=-(MASSMATRIX(1,6))*G0+HRM(6,5)
HRM1=HRM
IF(trim(adjustl(INPUTSTATIC))=="NO")THEN
        Aw=HRM(3,3)/rou/g0
        Sy=-HRM(3,5)/rou/g0
        XF=Sy/Aw/rou/g0
        hx=HRM(4,4)/vol/rou/g0
        hy=HRM(5,5)/vol/rou/g0
 ELSEIF(trim(adjustl(INPUTSTATIC))=="YES")THEN
 	    HRM(3,3)=rou*g0*Aw
	    HRM(3,5)=-rou*g0*Sy
	    HRM(5,3)=HRM(3,5)
	    HRM(4,4)=rou*g0*vol*hx
	    HRM(5,5)=rou*g0*vol*hy
 ENDIF
 
 DEALLOCATE(EU)
!----------------------------check file-------------------------------
write(12,'(a)')"[The hydrostatic restoring matrix]"
write(12,*)
write(12,'(a,e12.4)')"AW= ",AW
write(12,'(a,e12.4)')"Sy= ",Sy
write(12,'(a,f12.4)')"hx= ",hx
write(12,'(a,f12.4)')"hy= ",hy
write(12,'(a,e12.4)')"C33= ",HRM(3,3)
write(12,'(a,f12.4)')"C34=C43= ",HRM(3,4)
write(12,'(a,e12.4)')"C35=C53= ",HRM(3,5) 
write(12,'(a,e12.4)')"C44= ",HRM(4,4)
write(12,'(a,f12.4)')"C45=C54= ",HRM(4,5) 
write(12,'(a,f12.4)')"C46= ",HRM(4,6) 
write(12,'(a,e12.4)')"C55= ",HRM(5,5)
write(12,'(a,f12.4)')"C56= ",HRM(5,6) 
write(12,*)
!----------------------------------------------------------------------------------------
if((trim(adjustl(MoorLine))=="YES"))then       !读入锚链数据
    MOLINPUT=trim(adjustl(InAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.MOL'
    OPEN (27,FILE=MOLINPUT,STATUS="OLD")
	call IUTMP(27)
	READ(27,*)mlN1,mlN2,mlT
	mlH=H                                    !工作水深
	allocate(REF(1:mlN1),OMG(1:mlN2),EL1(1:mlN2),W1(1:mlN2))
    mlxg=xg
	mlzg=zg
	read(27,*)(tempcoor(i),i=1,3)
	allocate(xMl(1:mlN1),yMl(1:mlN1),zMl(1:mlN1))
	do	i=1,	mlN1
		read(27,*)xMl(i),yMl(i),zMl(i)
	enddo
	xMl(:)=xMl(:)+tempcoor(1)
	yMl(:)=yMl(:)+tempcoor(2)
	zMl(:)=zMl(:)+tempcoor(3)
	READ(27,*) (REF(I),I=1,mlN1)
	READ(27,*) (OMG(I),I=1,mlN2)
	READ(27,*) (W1(I),I=1,mlN2-1)
	READ(27,*)  (EL1(I),I=1,mlN2)

	CML=0.0
	allocate(SML(1:6,1:6,1:mlN1))
	SML=0.0

	do	i=1,	mlN1
		call ML(1,mlN2,mlT,mlH,mlXG,mlZG,XML(i),ZML(i),REF,OMG,EL1,W1,SML(:,:,i))
		CML(:,:)=CML(:,:)+SML(:,:,i)
	end	do
	write(12,*)
	write(12,'(a)')"The stiffness coefficients of the moor lines:"
	write(12,*)
	write(12,'(a,e12.4)')"CML33= ",CML(3,3)
	write(12,'(a,e12.4)')"CML35= ",CML(3,5)
	write(12,'(a,e12.4)')"CML53= ",CML(5,3)
	write(12,'(a,e12.4)')"CML55= ",CML(5,5)
	write(12,*)

	HRM=HRM+CML
  
    DEALLOCATE(REF,OMG,EL1,W1)
ENDIF
 deallocate(x00,ee)    
END SUBROUTINE RESTOMATRIX_SECT

!---------------------------------------------------------------------
SUBROUTINE MASSMATRIXPUTOUT()
IMPLICIT NONE
INTEGER(4)::I,J,k
write(12,"(a)")"[ Centre of gravity of the body in the mesh reference : ]"
write(12,"(a)")"[	XG	YG	 ZG  ]"
write(12,"(3(2x,f12.3))")XG,YG,ZG
write(12,*)
write(12,"(a)")"[ The Inertia Matrix of total floating body : ]"
do i=1,6
        write(12,"(6(e12.5,4x))")(massmatrix(i,j),j=1,6)
enddo
write(12,*)
write(12,"(a)")"[ The Inertia Matrix of each section : ]"
do k=1,nbsect
    write(12,"(a,4x,i8)")"Sect numeber",K
    do i=1,6
            write(12,"(6(e16.8,4x))")(SECTmassmatrix(i,j,K),j=1,6)
    enddo
enddo
end subroutine massmatrixputout


!************************************************************************************
!程序功能：获得浮体各计算剖面处的静水剪力及静水弯矩
!输入参数：SECTMASSMATRIX(1,1,J)为每段的质量（单位：吨）；
!                SECTCOG(1,J)为每段的重心纵向坐标（即重力作用点）   
!输出结果：每段的静水剪力SF(1:NBSECT)及静水弯矩SM(1:NBSECT)    
!输出文件：.hst文件
!BY SUN 2012.9.12   
!************************************************************************************
subroutine getHYDROSTATIC(X0,XALL)
implicit none

REAL(8),INTENT(IN)::X0,XALL
REAL(8),ALLOCATABLE,DIMENSION(:)::VOL_SECT
REAL(8),ALLOCATABLE,DIMENSION(:)::XB_SECT   !每段的浮心纵向坐标（即浮力作用点）
REAL(8),ALLOCATABLE,DIMENSION(:)::FB_SECT   !每段的浮力（单位KN）
REAL(8),ALLOCATABLE,DIMENSION(:)::SF,SM      !每段的静水剪力及静水弯矩
CHARACTER(LEN=300)::OUTPUTHST
integer::i,J,count
REAL(8),ALLOCATABLE,DIMENSION(:,:)::x00
REAL(8)::VOL_ALL,XB_ALL,SM_ALL,TEMP,TEMP1,TEMP2

ALLOCATE(VOL_SECT(1:NBSECT),XB_SECT(1:NBSECT),FB_SECT(1:NBSECT))
ALLOCATE(SF(1:NBSECT),SM(1:NBSECT))
ALLOCATE(x00(1,1:NB))
 VOL_SECT=0.0;XB_SECT=0.0

!DO I=1,NB
!     VOL_ALL=VOL_ALL-e(3,3,i)*xav(3,i)*ea(i)
!     XB_ALL=XB_ALL-e(3,3,i)*xav(3,i)*ea(i)*xav(1,i)
!ENDDO
!XB_ALL=XB_ALL/VOL_ALL
!SM_ALL=MASSMATRIX(1,1)*G0*(XALL-XG)-ROU*G0*VOL_ALL*(XALL-XG-XB_ALL)

do count=1,TarNum,Space
select case(count)
       case(1)
              x00(1,:)=xav(1,:)
        case(2)
              x00(1,:)=-xav(1,:)      
       case(3)
              x00(1,:)=-xav(1,:)
        case(4)
              x00(1,:)=xav(1,:)
 end select
DO J=1,NBSECT
    do i=1,NB
       IF(X00(1,I)<=SECTPSN(1,J))THEN
            VOL_SECT(J)= VOL_SECT(J)-e(3,3,i)*xav(3,i)*ea(i)
            XB_SECT(J)=XB_SECT(J)-e(3,3,i)*xav(3,i)*ea(i)*x00(1,i)
       ENDIF
    enddo
    FB_SECT(J)=ROU*G0*vol_SECT(J)
    SF(J)=SECTMASSMATRIX(1,1,J)*G0-FB_SECT(J)
  SM(J)=SECTMASSMATRIX(1,1,J)*G0*(SECTPSN(1,J)-SECTCOG(1,J))+ROU*G0*(XB_SECT(J)-VOL_SECT(J)*SECTPSN(1,J))
!    SM(J)=SM(J)-SM_ALL*(SECTPSN(1,J)+XG-X0)/(XALL-X0)
ENDDO
end do
!要求SECTPSN(:,NBSECT)为浮体的艏部终剖面
DO J=1,NBSECT
       SF(J)=SF(J)-SF(NBSECT)*(SECTPSN(1,J)+XG)/(SECTPSN(1,NBSECT)+XG)
       SM(J)=SM(J)-SM(NBSECT)*(SECTPSN(1,J)+XG)/(SECTPSN(1,NBSECT)+XG)
ENDDO
!-------------------输出.hst文件----------------------------------------------------
OUTPUTHST=trim(adjustl(OutAccess))//'\'//TRIM(ADJUSTL(PROJNAME)) //'.hst'
OPEN (262,FILE=OUTPUTHST)
write(262,'(a)')"#              [ WALCS Version 2.0 ]"
write(262,*)
write(262,'(a)')"#              [ 静水载荷计算结果 ]"
write(262,*)
write(262,'(a,5x,a,3x,a)')"#","静水剪力(KN)","静水弯矩(KN*m)"
do j=1,Nbsect
!	write(262,'(i3,E13.4,E13.4)')j,SF(J),SM(J)
	write(262,'(i3,E13.4,E13.4,E13.4)')j,SF(J),SM(J)
enddo


DEALLOCATE(VOL_SECT,XB_SECT,FB_SECT,SF,SM)

end subroutine getHYDROSTATIC


END MODULE GETM_MOD
!**************************************************************************************************************************************