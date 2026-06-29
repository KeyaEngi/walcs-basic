!--------------------------------------------------------------
!模块功能：计算指数积分
!程序编制：孙葳 2011.12.26
!--------------------------------------------------------------
MODULE DEI_MOD   
                 
PRIVATE          
PUBLIC::DEI      
                 
CONTAINS         
!------------------------------------
FUNCTION DEI(X)                            
DOUBLE PRECISION DEI,X,Y,R,G               
! 	EXTERNAL EIX                             
! 	DOUBLE PRECISION EIX                     
IF (X<0.0) THEN                         
  DEI=EX(X)                             
  RETURN                                   
END IF                                     
Y=-X                                       
R=0.57721566490153286060651                
CALL FLRGS(0.0D0,Y,EIX,1.0E-07,G)          
DEI=G                                      
DEI=R+LOG(X)+DEI                           
RETURN                                     
END FUNCTION DEI                           
!---------------------------                                           
FUNCTION EIX(X)                            
DOUBLE PRECISION EIX,X                     
EIX=(EXP(-X)-1.0)/X                        
RETURN                                     
END FUNCTION EIX   
!---------------------------

	SUBROUTINE FLRGS(A,B,F,EPS,G)                                
	DIMENSION T(5),C(5)                                          
	DOUBLE PRECISION A,B,F,G,T,C,S,P,H,AA,BB,W,X,Q               
	DATA T/-0.9061798459,-0.5384693101,0.0,0.5384693101,0.9061798459/                     
	DATA C/0.2369268851,0.4786286705,0.5688888889,0.4786286705,0.2369268851/                      
	M=1                                                          
	S=(B-A)*0.001                                                
	P=0.0                                                        
10	H=(B-A)/M                                                  
	G=0.0                                                        
	DO 30 I=1,M                                                  
	  AA=A+(I-1)*H                                               
	  BB=A+I*H                                                   
	  W=0.0                                                      
	  DO 20 J=1,5                                                
	    X=((BB-AA)*T(J)+(BB+AA))/2.0                             
	    W=W+F(X)*C(J)                                            
20	  CONTINUE                                                 
	  G=G+W                                                      
30	CONTINUE                                                   
	G=G*H/2.0                                                    
	Q=ABS(G-P)/(1.0+ABS(G))                                      
	IF ((Q.GE.EPS).AND.(ABS(H).GT.ABS(S))) THEN                  
	  P=G                                                        
	  M=M+1                                                      
	  GOTO 10                                                    
	END IF                                                       
	RETURN                                                       
	END SUBROUTINE FLRGS                                         
!------------------------------------------------                                                               
FUNCTION EX(X) RESULT(X_RESULT)
IMPLICIT NONE
REAL(8)::X,X_RESULT
REAL(8)::R1,R2,R3,R4,R5,R6
REAL(8)::T1,T2,T3,T4
REAL(8)::S1,S2,S3,S4
REAL(8)::TEMP1,TEMP2,y

DATA R1,R2,R3,R4,R5,R6/0.107857D-2,-0.976004D-2,0.5519968D-1,-0.24991055D0,0.99999193D0,-0.57721566D0/
DATA T1,T2,T3,T4/8.5733287401D0,1.8059016973D1,8.6347608925D0,0.2677737343D0/
DATA S1,S2,S3,S4/9.5733223454D0,2.56329561486D1,2.10996530827D1,3.9584969228D0/

X_RESULT=0.
y=ABS(x)

IF(0<Y<1) THEN
   X_RESULT=-(R6+R5*y+R4*(y**2)+R3*(y**3)+R2*(y**4)+R1*(y**5)-LOG(Y))
ELSE IF(Y>=1) THEN
   TEMP1=T4+T3*Y+T2*(Y**2)+T1*(Y**3)+Y**4
   TEMP2=S4+S3*Y+S2*(Y**2)+S1*(Y**3)+Y**4
   X_RESULT=-((TEMP1/TEMP2)/(EXP(Y)*Y))
END IF

END FUNCTION EX
!-----------------------------------------------------
                             

END MODULE DEI_MOD     
                        

