MODULE MAINCAL_MOD  
  
USE MESH_MOD,ONLY:SL,xn,xav,e,EA,xq,NPANEL,NB,Factor,TarNum,Space  !Factor,TarNum,Space考虑对称性所新添加的变量。
USE ENVIRONMENT_MOD,ONLY:NBU,NBHEAD,NBOME,PHASE,HEAD,OME0,UKN,OME,AMP,G0,PI,H,NBSECT,OME1
USE PVPOINT_MOD,ONLY:NPPOINT,NFEPOINT
USE CAL_MOD,ONLY:MASSSOLVE,WATERDEPTH,HULLTYPE
USE PRINT_MOD,ONLY:PRTMOTION,PRTSECTLOAD,PRTPPRESS,PRTFEPRESS,PRTva,BILGEFORCE
USE WAVEKJ_MOD
USE BASICSOURCE_MOD
USE GREEN_MOD
USE SOLVEPHI_MOD
USE getHydroCoef_MOD
USE SHIPMOTION_MOD
USE SHIPPRESS_MOD
USE SECTLOAD_MOD
USE FEMPRESS_MOD
USE OUTPUTFILE_MOD
USE POINTVA_MOD
USE FEMPHIPRESS_MOD
USE POINTPRESS_MOD
!USE BILGEFORCE_MOD   !舭龙骨节点力计算模块
USE CPUTime_MOD

IMPLICIT NONE
PRIVATE
PUBLIC::MAINCAL

CONTAINS
!----------------------------------------------------------
SUBROUTINE MAINCAL()
IMPLICIT NONE

INTEGER::IU,IB,IP  !（用于计数用）航速、浪向、频率
REAL(8),DIMENSION(0:NBOME)::k0
REAL(8)::PNU0,PNU
REAL(8)::OMEIN                                      !插值频率  byshiyuyun 2013.12.12                        
integer(4),	parameter::	wkjNum=40
real(8),dimension(1:wkjNum)	::	wkj				!波数
real(8),ALLOCATABLE,dimension(:,:,:,:)::SRPHI	!诱导速度及单层势
real(8),ALLOCATABLE,dimension(:,:,:)::DRPHI	  !诱导速度法向投影
real(8),ALLOCATABLE,DIMENSION(:,:)::GcD,GsD     !格林函数法向导数
real(8),ALLOCATABLE,DIMENSION(:,:,:)::GcS,GsS   !格林函数梯度及本身
real(8),ALLOCATABLE,DIMENSION(:,:,:,:)::G0Dc ,G0Ds !格林函数本身以及导数
real(8),allocatable,dimension(:,:,:,:,:)::Phi   !速度势
real(8),allocatable,dimension(:,:,:)::F       !水动力系数及波浪干扰力
real(8),allocatable,dimension(:,:,:,:)::sectf !局部水动力系数及波浪干扰力
real(8),allocatable,dimension(:,:,:)::MOTION    !船舶六自由度运动
real(8),allocatable,dimension(:,:)::MOTIONPOLY    !船舶POLY
real(8),allocatable,dimension(:,:)::M4x    !附加横摇阻尼力矩
real(8),allocatable,dimension(:,:,:)::ps,pw,pd,pr,pall	!船舶表面压力
real(8),allocatable,dimension(:,:,:,:)::Phi_press   !速度势
real(8),ALLOCATABLE,DIMENSION(:,:,:,:)::SSTRENGTH     !格林函数法向导数
real(8),ALLOCATABLE,DIMENSION(:,:,:,:)::S_PRESS
integer(4)::i,j
!real(8)::M4x(1:2)  !附加横摇阻尼力矩
real(8)::bv44,PGL
REAL::RATIO,RATIOPTS
real(8)::ometh

write(*,'(3(a,i3))')" 航速计算数目=",NBU,"    浪向角计算数目=",NBHEAD,"    自然频率计算数目=",NBOME
write(*,*)
write(*,'(2(a,f5.1),a)')" 浪向角计算从 ",HEAD(1)," 度 到 ",HEAD(NBHEAD)," 度"
write(*,*)
write(*,'(2(a,f7.3),a)')" 自然频率计算从 ",OME0(1)," rad/s 到",OME0(NBOME)," rad/s"

!****************************SHIYUYUN***********************************
if(TRIM(ADJUSTL(PRTMOTION))=="YES")THEN
ALLOCATE(MOTION(1:2,1:6,0:NBOME),MOTIONPOLY(1:2,1:6),M4x(1:2,0:NBOME))
MOTION=0.0;MOTIONPOLY=0.0;M4x=0.0
END IF
!***********************************************************************
if(TRIM(ADJUSTL(PRTFEPRESS))/="NO")THEN
        CALL FEPPOINT(RATIO)    !读入有限元网格文件并处理得到有限元网格信息
        write(*,*)
        write(*,'(a,i7)')" 有限元网格数目:",NFEPOINT
ENDIF

if(TRIM(ADJUSTL(PRTVA))=="YES")CALL GETVAPOINT()   !读入速度加速度计算点坐标文件.pva
if(TRIM(ADJUSTL(PRTPPRESS))=="YES") then !读入压力计算点坐标文件.pts
     CALL READPPOINT(RATIOPTS)
     ALLOCATE(S_PRESS(1:NPPoint,1:NPanel,1:4,1:Factor))
	 call GR_P(s_press)     !计算基本源
end if
!if(TRIM(ADJUSTL(BILGEFORCE))=="YES")CALL GETBILGENODE(PGL) !读入舭龙骨附近定义的节点坐标文件.pvf

ALLOCATE(SRPHI(1:NPANEL,1:NPANEL,1:4,1:Factor),DRPHI(1:NPANEL,1:NPANEL,1:Factor))
CALL OUTPUTFILE()
CALL GR(SRPHI,DRPHI)   !计算基本源



ometh=sqrt(2.0*PI*g0/SL) 	

do IU=1,NBU	! -----------------------对航速循环---------------------------

	UKN(IU)=UKN(IU)*0.5144  !航速单位转换为m/s

	do IB=1,NBHEAD	 !---------------------对浪向角循环---------------------------

		   HEAD(IB)=HEAD(IB)/180.0*PI
           write(*,*)      		
           Call GetStartTime()
           do IP=0, NBOME	!----------------------对波浪频率循环------------------------- 
                if(IP>=1)     then      
			        write(*,'(1x,I4,2x,a,f7.3,a,f7.1,a,f7.3,a)')(IB-1)*NBOME+IP," 航速= ",UKN(IU)/0.5144," 节，     浪向角= ",HEAD(IB)*180.0/PI," 度，     自然频率= ",ome0(IP)," rad/s"
			    endif
                !-----------------------------------------------------------------------------------------------------
                !***********************SHIYUYUN2014***********************
                !				EVERY UKN AND HEAD CHANGES
                
                IF(IP==0)THEN
				    IF((HEAD(IB)/PI*180<90.0).AND.(UKN(IU)/=0.0))THEN !如果浪向角小于90度且航速不为0,可能存在遭遇频率为负数的情况
				         ome0(IP)=(SQRT(g0**2+4*g0*UKN(IU)*cos(HEAD(IB))*ome(IP))-g0)/(2.0*UKN(IU)*cos(HEAD(IB))) !g0=9.81
				    ELSEIF((HEAD(IB)/PI*180>90.0).AND.(UKN(IU)/=0.0))THEN
				        ome0(IP)=(-SQRT(g0**2+4*g0*UKN(IU)*cos(HEAD(IB))*ome(IP))-g0)/(2.0*UKN(IU)*cos(HEAD(IB)))
				        IF(ome0(IP)<0)THEN
				            ome0(IP)=(-SQRT(g0**2+4*g0*UKN(IU)*cos(HEAD(IB))*ome(IP))-g0)/(2.0*UKN(IU)*cos(HEAD(IB)))
				        END IF
				    ELSE
				        ome0(IP)=ome(IP)
				    END IF
				END IF
            
!                IF(IP==0.and.UKN(IU)>0)THEN
!                    ome0(IP)=(-SQRT(g0**2+4*g0*UKN(IU)*cos(HEAD(IB))*ome(IP))-g0)/(2.0*UKN(IU)*cos(HEAD(IB)))
!                else
!                    ome0(IP)=ome(IP)
!                ENDIF
                
                !**********************************************************	 
				Pnu0=ome0(IP)**2/g0 !无限水深色散关系 Pnuo=波数
				if((trim(adjustl(WATERDEPTH))=="INFINITE"))then  !无限水深
					wkj(1)=Pnu0 !波数
				elseif((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
					call wavekj(wkjNum,wkj,H,Pnu0) !求自然频率波数
				endif
				k0(IP)=wkj(1)
				ome(IP)=ome0(IP)+wkj(1)*UKN(IU)*cos(HEAD(IB)) !转化成遭遇频率
				Pnu=ome(IP)**2/g0 !转化成遭遇频率后的波数Pnu
				OME1(IP)=ome(IP)
				ome(IP)=ABS(ome(IP))
				if((trim(adjustl(WATERDEPTH))=="FINITE"))then !有限水深
					call wavekj(wkjNum,wkj,H,Pnu) !求遭遇频率波数
				endif
                !------------------------------------------------------	
                
                allocate(GcS(1:NPANEL,1:NPANEL,1:4))
                allocate(GsS(1:NPANEL,1:NPANEL,1:4))	
                
                allocate(GcD(1:NPANEL,1:NPANEL),GsD(1:NPANEL,1:NPANEL))		
			    allocate(G0Dc(1:NPANEL,1:NPANEL,1:3,1:Factor),G0Ds(1:NPANEL,1:NPANEL,1:3,1:Factor))					
				call coef(wkjNum,wkj,Pnu,SRPHI,DRPHI,GcS,GsS,GcD,GsD,G0Dc,G0Ds) !格林函数计算
				DEALLOCATE(GcS,GsS,GcD,GsD)
			    ALLOCATE(PHI(1:2,1:2,1:7,1:Factor,1:NB))            
                ALLOCATE(SSTRENGTH(1:2,1:7,1:Factor,1:NPANEL))
				call solve_gmres(G0Dc,G0Ds,ome0(IP),HEAD(IB),k0(IP),Pnu0,Phi,SSTRENGTH) !GMRES法求解速度势
				DEALLOCATE(G0Dc,G0Ds)
				IF(trim(adjustl(PRTPPRESS))=="YES") THEN
				    ALLOCATE(PHI_PRESS(1:2,1:2,1:7,1:NPPoint))
                    CALL GetPressPhi(wkjNum,wkj,Pnu,S_PRESS,SSTRENGTH,Phi_press)
				END IF
				
                DEALLOCATE(SSTRENGTH)
													                                       
              !-----------------------------求解水动力系数及波浪干扰力---------------------------------
			   ALLOCATE(F(1:2,1:6,1:9))
			   if(trim(adjustl(HULLTYPE))=="SINGLE")then
		   		   call getHydroCoef(UKN(IU),ea(1:NB),e(:,:,1:NB),xav(:,1:NB),ome(IP),phi(:,:,:,:,1:NB),F(:,:,1:7)) !计算水动力系数
                   call getWaveExistForce(OME(IP),OME0(IP),AMP(IP),k0(IP),HEAD(IB),ea(1:NB),e(:,:,1:NB),xav(:,1:NB),F(:,:,1:9),IP) !计算波浪干扰力
               ELSEif(TRIM(ADJUSTL(HULLTYPE))=="MULTI")	then	
                    allocate(sectf(1:2,1:6,1:nbsect,1:9))
                    call getHydroCoef_sect(UKN(IU),ea(1:NB),e(:,:,1:NB),xav(:,1:NB),ome(IP),phi(:,:,:,:,1:NB),F(:,:,1:7),sectf(:,:,:,1:7))
                    call getWaveExistForce_sect(OME(IP),OME0(IP),AMP(IP),k0(IP),HEAD(IB),ea(1:NB),e(:,:,1:NB),xav(:,1:NB),F(:,:,1:9),sectf(:,:,:,1:9),IP)
               ENDIF	   	   			       

		      !------求解运动--------------频率<0.1进行插值求解运动byshiyuyun2013.12.12----------------------------------------------------------------		       
		        !***************************SHIYUYUN********************************
		       if(trim(adjustl(PRTMOTION))=="YES")	then
                !	motion(:,:,IP)=0.0
					IF(IP==0)THEN
					    call shipmotion	(HEAD(IB),OME(IP),F(:,:,:),motionPOLY(:,:),UKN(IU),M4x(:,IP),bv44)  
				    ELSE 
				        call shipmotion	(HEAD(IB),OME(IP),F(:,:,:),motion(:,:,IP),UKN(IU),M4x(:,IP),bv44)
				    END IF

                !*******************************************************************
		      !------求解船舶表面压力--------------------------------------------------------------------------
			           ALLOCATE(ps(1:2,1:Factor,1:NB),pw(1:2,1:Factor,1:NB),pd(1:2,1:Factor,1:NB),pr(1:2,1:Factor,1:NB),pall(1:2,1:Factor,1:NB))			
			           call shippress(ome(IP),ome1(IP),UKN(IU),HEAD(IB),amp(IP),k0(IP),xav(:,1:NB),phi(:,:,:,:,1:NB),motion(:,:,IP),ps,pw,pd,pr,pall,IP)
				   endif

				   DEALLOCATE(F,PHI)
		      !-------------------------------求解船舶各站剖面载荷-----------------------------------------------------
			       if(trim(adjustl(PRTSECTLOAD))=="YES")then
					        if(TRIM(ADJUSTL(MASSSOLVE))=="SECT")	then				
					        	call ShipSectLoad (HEAD(IB),e(:,:,1:NB),ea(1:NB),pall,xav(:,1:NB),OME0(IP),ome(IP),ome1(IP),PHASE(IP),AMP(IP),motion(:,:,IP),M4x(:,IP),IP)
							ELSEIF(TRIM(ADJUSTL(HULLTYPE))=="MULTI") then
							    call shipsectload_lt (HEAD(IB),OME0(IP),ome(IP),ome1(IP),PHASE(IP),AMP(IP),motion(:,:,IP),bv44,sectf,IP)
					        endif											
			       endif
		           if(TRIM(ADJUSTL(HULLTYPE))=="MULTI")DEALLOCATE(SECTF)
		      !-------------------------------求解有限元网格上的波浪压力----------------------------    
			       if(TRIM(ADJUSTL(PRTFEPRESS))/="NO")	then			
                        CALL FEPRESS(IU,IB,IP,HEAD(IB),OME0(IP),PHASE(IP),PS,PW,PD,PR,PALL,RATIO)
			       endif
                   if(trim(adjustl(PRTMOTION))=="YES")DEALLOCATE(ps,pw,pd,pr,pall)

              !------求解给定点压力-----------------------------------------------------
			       IF(trim(adjustl(PRTPPRESS))=="YES") THEN
						call shippointpress(ome0(IP),ome(IP),ome1(IP),UKN(IU),HEAD(IB),amp(IP),k0(IP),phi_press(:,:,:,:),&
											motion(:,:,IP),RATIOPTS,PHASE(IP),IP)
						DEALLOCATE(PHI_PRESS)
				   end if 
		      !------求解给定点速度加速度-----------------------------------------------------
			       if(trim(adjustl(PRTva))=="YES")then
				          call shippointva(HEAD(IB),OME0(IP),ome(IP),ome1(IP),K0(IP),PHASE(IP),motion(:,:,IP),IP)
			       endif
			  !--------------------求舭龙骨附近位置的等效节点力----------------------------------
!			      if(TRIM(ADJUSTL(BILGEFORCE))=="YES")then
!			              CALL  GETBILGEFORCE (HEAD(IB),OME0(IP),OME(IP),PHASE(IP),MOTION(:,4,IP),bv44,PGL)
!			      endif		

				 IF((abs(Ome(IP))<0.15).AND.Ome0(IP)>ometh.and.UKN(IU)>0) THEN
				    
				        motion(:,:,IP)=1.0/0.15*abs(ome(IP))*motionPOLY(:,:)
				        
				 END IF

		   END DO  !结束频率循环
		   
		   
           Call GetEndTime()
           if(trim(adjustl(PRTMOTION))=="YES")call m_putouts(IB,motion,M4x) 

	END DO  !结束浪向循环

END DO  !结束航速循环


deallocate(xn,xav,e,EA,xq)
DEALLOCATE(UKN,HEAD,OME0,AMP,PHASE,OME)
if(TRIM(ADJUSTL(PRTMOTION))=="YES")DEALLOCATE(MOTION)
if(TRIM(ADJUSTL(PRTPPRESS))=="YES")DEALLOCATE(S_PRESS)
if(TRIM(ADJUSTL(PRTFEPRESS))=="PATRAN")	then
	open(301,file='tempp.dat',status='old',dispose='delete')
endif

write(*,*)
write(*,*)"规则波中的所有计算已经完成！"
write(*,*)
write(*,*)"感谢您使用三维波浪载荷计算系统(COMPASS-WALCS-BASIC V1.0)."
write(*,*)
close(12)

END SUBROUTINE MAINCAL
!------------------------------------------------------------


END MODULE MAINCAL_MOD
