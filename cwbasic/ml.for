!************************************************************************************
!程序功能：计算锚链刚度系数，仅考虑垂荡，纵摇耦合情况
!输入参数：T-单根锚链预张力，N1-锚链根数，N2-单根锚链分段数,H-锚链工作水深
!          XG-系泊油船重心与其艉间的距离
!          ZG-系泊油船重心与其基线间的距离
!          XL-转塔与系泊油船艉间的距离
!          ZL-着链点与系泊油船基线间的距离，若在基线下，则为负
!          (REF(I),I=1,N1)-每根锚链与系泊油船对称轴线间的夹角，按逆时针方向选取
!          (OMG(I),I=1,N2)-单根锚链中每段的水中重量(顺序由上至下)
!          (W1(I),I=1,N2-1)- 单根锚链中每个配重块的水中重量(顺序由上至下)，共N2-1块
!          (EL1(I),I=1,N2)- 单根锚链中每段的长度(顺序由上至下)
!输出参数：CML-锚链刚度系数
!备注：
!程序编制：张海彬    时间：2002年12月        
!************************************************************************************
	 SUBROUTINE  ML(N1,N2,T,h,XG,ZG,XL,ZL,REF,OMG,EL1,W1,CML)
C      计算刚度系数，仅考虑垂荡，纵摇耦合情况(song)

        IMPLICIT REAL*8 (A-H,O-Z)
!------Add by zhb------------------------
      integer(4),intent(in)::N1,N2
	real(8),intent(in)::T,h,XG,ZG,XL,ZL
	real(8),intent(inout),dimension(1:N1)::REF
	real(8),intent(in),dimension(1:N2)::OMG,EL1
	real(8),intent(inout),dimension(1:N2)::W1
	real(8),intent(out),dimension(1:6,1:6)::CML
!----------------------------------------     
  	DIMENSION    TW1(20),TW2(20),          
     *     s(20),S0(20) !,ref(20),el1(20),W1(20),OMG(20)
      common/COE1/RXX,RXY,RYX,RYY
 	COMMON/FUN/IC,NT,TH,TF,AWI(10),AC(10),AC0(10)
  	COMMON/PA/AW(10),ASL(10),TLM0
!	COMMON/AMT/XL,LL
	COMMON/AMM/TM(6,6) 
C      XG MEANS THE DISTANCE BETWEEN SHIP COG AND BASELINE	  
C      ZG MEANS THE DISTANCE  BETWEEN  SHIP  COG AND STERN
C      XL MEANS THE DISTANCE  BETWEEN  TURRET AND STERN
C      ZL MEANS THE DISTANCE  BETWEEN  MOORLINE  AND BASELINE      
C      IF MOORLINE UNDER BASELINE ,ZL<0
      
C      应设W1(N2)=0
C          OPEN(101,FILE='RLINE2.DAT')

!		READ(101,*)T,n1,n2,h
!		READ(101,*)xg,zg
!		READ(101,*)xl,zl 
		
!	    READ(101,*) (REF(I),I=1,N1)
!	    READ(101,*) (OMG(I),I=1,N2)
!	    READ(101,*) (W1(I),I=1,N2-1)
!	    READ(101,*)  (EL1(I),I=1,N2)
C	    CLOSE(101)
C	N1:MOOR LINE NUMBERS
C     N2: HOW MANY PARTS OF A SINGLE MOOR LINE

	w1(n2)=0.0
	PI=3.1415926
 	DEG=PI/180
      X0=XL-XG
 	Z0=ZL-ZG


C 	OPEN(5,FILE='WLINE.DAT')
! 	WRITE(5,*)'MOORLINE PARAMETERS'
! 	WRITE(5,13) t ,n1, n2, h
! 13	FORMAT('T n1 n2 h ',1F10.3,2X,I2,2X,I2,2X,1F10.3)

!	WRITE(5,23)xg,zg
 !23	FORMAT('xg zg ',1F10.3,2X,1F10.3)
!
!	WRITE(5,33)xl,zl
! 33	FORMAT('xl zl ',1F10.3,2X,1F10.3) 
		
!	WRITE(5,*) '(REF(I),I=1,N1)'
!	WRITE(5,5) (REF(I),I=1,N1)
 	
!	WRITE(5,*)' (OMG(I),I=1,N2) '
!	WRITE(5,5) (OMG(I),I=1,N2)

!	WRITE(5,*) '(W1(I),I=1,N2)'
!	WRITE(5,5) (W1(I),I=1,N2)
!	WRITE(5,*)  '(EL1(I),I=1,N2)'
 !     WRITE(5,5)  (EL1(I),I=1,N2)
	w1(n2)=0
 	DO 225 I=1,N1
	REF(I)=REF(I)*DEG
225   CONTINUE	  
 


	  

C1      预张力列表
C11      情况1：锚链触海底       
      	DO 25 J=N2,1,-1
          S(J)=OMG(J)*EL1(J)
      	S0(J)=0.0
	IF(J.GT.1) THEN
	   DO 20  I=J-1,1,-1
	S(I)=S(I+1)+OMG(I)*EL1(I)+W1(I)
20	S0(I)=S0(I+1)+OMG(I+1)*EL1(I+1)+W1(I)
      END IF
		                
25	TW1(J)= ACO(J,H,S,S0,OMG)
!	WRITE(5,*) '(TW1(I),I=1,N2)'
!	WRITE(5,5) (TW1(I),I=1,N2)
!5	FORMAT(4F10.4)
 


	
C12     情况2 ：块触海底     
	  DO 35 J=N2,1,-1
      S(J)=OMG(J)*EL1(J)+W1(J)
      S0(J)=W1(J)
 	IF(J.GT.1) THEN
	DO 30  I=J-1,1,-1
    	S(I)=S(I+1)+OMG(I)*EL1(I)+W1(I)
30	S0(I)=S0(I+1)+OMG(I+1)*EL1(I+1)+W1(I)
      END IF

35	TW2(J)=ACO(J,H,S,S0,OMG)
! 	WRITE(5,*) '(TW2(I),I=1,N2)'
! 	WRITE(5,5) (TW2(I),I=1,N2)
		 
C2     求解与给定预张力对应的锚链状态及刚度系数
      I=0
40	I=I+1
	N=I
	if(n.gt.n2) stop 'error data'
	IF (T.LE.0.9*TW1(I)) GOTO 50
	IF (T.LE.TW2(I))  GOTO 80
	GOTO 40
 50   TLM0=0
      TH=H
	TF=T
	IC=1
	NT=N
	IF(NT.GT.1) THEN
	DO 55  I=NT-1,1,-1
 55   AW(NT-I)=W1(I)
      DO  60  I=NT-1,1,-1
60	ASL(NT+1-I)=EL1(I)
      ENDIF
	DO 65 I=NT,1,-1
65	AWI(NT+1-I)=OMG(i)
	CALL M_LINE
	
	GOTO 140
	 
 80   TLM0=0
      TH=H
 	TF=T
	IC=2
	NT=N
	IF(NT.GT.1) THEN
	DO 85  I=NT-1,1,-1
 85   AW(NT-I)=W1(I)
      DO  90  I=NT-1,1,-1
90	ASL(NT+1-I)=EL1(I)
      ENDIF
	DO 95 I=NT,1,-1
95	AWI(NT+1-I)=OMG(i)
      DO 100 I=1,NT
100	TLM0=TLM0+EL1(I)
	CALL M_LINE
	
140   CONTINUE
	
	A=RXX
	B=RXY
	C=RYX
	D=RYY

	DO 145  I=1,4
	DO 145  J=1,4
	TM(I,J)=0.0
145   CONTINUE
	 
     	 
	DO  150  I=1,N1
 	TM1=D
	TM2=-D*X0-C*Z0*COS(REF(I))
	TM(1,1)=-TM1+TM(1,1)
	
	TM(1,3)=-TM2+TM(1,3)
	

	TM(3,1)=TM1*X0+TM(3,1)

	TM(3,3)=TM2*X0+TM(3,3)

150   CONTINUE
      
      TM(2,2)=TM(1,1)
	TM(2,4)=TM(1,3)
	TM(4,2)=TM(3,1)
	TM(4,4)=TM(3,3)
!-----Add by zhb-----------------------
      CML(3,3)=TM(1,1)
	CML(3,5)=TM(1,3)
	CML(5,3)=TM(3,1)
	CML(5,5)=TM(3,3)
!-----------------------------
!	write(5,*) '[(TM(I,J),J=1,4),I=1,4]'
!	DO 89 I=1,4
!	WRITE(5,195)(TM(I,J),J=1,4)
!89    CONTINUE
!195    FORMAT(4F18.4)
C      CLOSE(5)
!      WRITE(*,*)'TM'
      RETURN
	END

	FUNCTION   ACO(N2,H,S,S0,OMG)
	IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION  S(N2),S0(N2),OMG(N2)
	T=1.5*S(1)
10    G=0.0
	DO 15  I=1,N2
	G=G+(SQRT(T**2-S(1)**2+S(I)**2)-
     *     SQRT(T**2-S(1)**2+S0(I)**2))/OMG(I)
15    CONTINUE

      F=G-H
	 
	IF(ABS(F).LE.1.0-6) GOTO 25
	SW=T
	G=0.0
	DO 20 I=1,N2
	G=G+T*(1/SQRT(T**2-S(1)**2+S(I)**2)-
     *   1/SQRT(T**2-S(1)**2+S0(I)**2))/OMG(I)
 20   CONTINUE
	T=T-F/G
	IF(ABS(T-SW).GT.1.0E-5) GOTO  10
25    ACO=T
      RETURN
	END
   	
	subroutine  M_LINE
C       Mooring Line
C         shape parameter {l,s,h,m,a}
C         stiffness coefficient {Rxx,Rxy=Ryx,Ryy}
C       IC=1 :: {mo,h,T}       IC=2 :: {l,h,T}
C       IC=1 :: dm0=0          IC=2 :: dl=0
C.      TS,TH,TL,TM,TM0,TN,TN0,TF,TFH,TFV, NT
C.      AW(1:NT-1),AWI(1:NT), ASS(1:NT),ASH(1:NT),ASL(1:NT),
C.      ASM0(1:NT),ASM(1:NT),ASN0(1:NT),ASN(1:NT)
C.      RXX,RXY,RYX,RYY
C       DATE : 2000.11.18
	IMPLICIT REAL*8 (A-H,O-Z)
        common/COE1/RXX,RXY,RYX,RYY
	COMMON/FUN/IC,NT,TH,TF,AWI(10),AC(10),AC0(10)
	COMMON/PA/AW(10),ASL(10),TLM0
C       CHARACTER RF*10,WF*10,FI*7,F*2
C        WRITE(*,1)
C1       FORMAT(/5X,'[OPERAND FILE NO.] =___')
C        READ(*,'(A)') F
C        WRITE(*,2) F
C2       FORMAT(/10X,'[FILE NO.] =',A3,10X,'OK!')
C        L=LEN_TRIM(F)
C        FI='ML'//F(1:L)//'.'
C        LI=LEN_TRIM(FI)
C        RF='R'//FI(1:LI)//'DAT'
C        WF='W'//FI(1:LI)//'DAT'
C        OPEN(1,FILE=RF,STATUS='OLD')
C        OPEN(2,FILE=WF)
C        WRITE(5,3) RF
C3       FORMAT(/20X,'[Input-File]: ',A11)
        CALL MOORINGLINE
C        CLOSE(1)
C        CLOSE(2)
        RETURN
        END
C-----------------------------------------------------------------------
C***********************************************************************
        SUBROUTINE MOORINGLINE
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER REM*80

        INTEGER IC,IT,NT
        REAL*8 TS,TH,TL,TM0,TM,TN0,TN,TF,TFH,TFV
        DIMENSION ASS(10),ASH(10),
     *  ASM0(10),ASM(10),ASN0(10),ASN(10)
        REAL*8 RXX,RXY,RYX,RYY
        
        COMMON /FUN/IC,NT,TH,TF,AWI(10),AC(10),AC0(10)
	COMMON/PA/AW(10),ASL(10),TLM0

	  common/COE1/RXX,RXY,RYX,RYY
        CALL INPUT0(REM,IC,IT,NT,TH,TL,TM0,TN0,TF,AW,AWI,TLM0,ASL)
100     CALL ACAC0(IC,NT,TF,TN0,ASL,AW,AWI,AC,AC0)
!        WRITE(5,1) IC,NT,(I,AC0(I),AC(I),I=1,NT)
!1       FORMAT(//2X,'IC=',I2,5X,'NT=',I2,9X,'AC0(),AC():',/(I5,2F20.5))
        IF(IT.EQ.-1.AND.NT.EQ.1) THEN
          IF(IC.EQ.1) THEN
            TN=TN0/(1-AWI(1)*TH/TF)
            TM=SQRT(TN**2-1)
            TFH=TF/TN
            TL=TFH/AWI(1)*(TM-TM0)
            ASL(1)=TL
          ELSE
            W=(TH+AWI(1)*(TL**2-TH**2)/TF/2)/TL
            IF(W.GE.1) STOP 'Error! { w*(l+h)/2 >= T }'
            TN=1/SQRT(1-W**2)
            TM=SQRT(TN**2-1)
            TN0=TN*(1-AWI(1)*TH/TF)
            TM0=SQRT(TN0**2-1)
            TFH=TF/TN
            ASL(1)=TL
          END IF
        ELSE
          X0=1
          DX=.1
          CALL SEARCH(X0,DX,X,F,NL)
          TN=X
          TM=SQRT(TN**2-1)
          TFH=TF/TN
          IF(IC.EQ.1) THEN
            SM1=TM-TN*AC(NT)
            SL1=TFH/AWI(1)*(SM1-TM0)
            ASL(1)=SL1
            TL=0
            DO 10 I=1,NT
10          TL=TL+ASL(1)
          ELSE
            TM0=TM-TN*AC(NT)
            TN0=SQRT(1+TM0**2)
          END IF
        END IF
        CALL PARAMETERS(NT,TFH,TM0,ASL,AW,AWI,ASH,ASS,ASM0,ASM,ASN0,
     *  ASN,TS)
        TFV=TFH*TM
!        WRITE(5,2) (I,ASL(I),ASH(I),ASS(I),ASM0(I),ASM(I),ASN0(I),
!     *  ASN(I),I=1,NT)
!2       FORMAT(//2X,'PARAMETERS:{l,h,s,m0,m,n0,n}',/(I3,3F10.3,4F11.5))
        CALL STIFFNESS(NT,TFH,TS,ASH,ASL,ASM0,ASN0,ASN,RXX,RXY,RYX,RYY)
!        WRITE(5,3) TF,TFH,TFV,TH,TS,TL,TM0,TM,TN0,TN,RXX,RXY,RYX,RYY
!3       FORMAT(//2X,'TF,TFH,TFV=',3F12.5,/4X,'TH,TS,TL=',3F12.3,
!     *  /4X,'TM0,TM,TN0,TN=',4F15.5,/2X,'RXX,RXY,RYX,RYY=',4F15.5)
C200     READ(1,*) TF
         TF=0
        IF(TF.GT.0) THEN
!          WRITE(5,4) TF
!4         FORMAT(//,80('='),/25X,'TF=',F12.5)
          GOTO 100
        END IF
        RETURN
        END
C***********************************************************************
        SUBROUTINE INPUT0(REM,IC,IT,NT,TH,TL,TM0,TN0,TF,AW,AWI,TLM0,ASL)
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER REM*80,TIT(2)*3
        INTEGER IC,IT,NT
        REAL*8 TH,TL,TM0,TF
        DIMENSION ASL(10),AWI(10),AW(10)
	

        DATA TIT/'TM0','TL'/
C        READ(1,'(A)') REM
C        READ(1,*) IC,NT,TH,TLM0,TF
        REM='WENCHANG'

        IT=NT
        NT=ABS(NT)
C       WRITE(*,*) REM
        IF(IC.EQ.1) THEN
          TM0=TLM0
          TN0=SQRT(1+TM0**2)
        ELSE IF(IC.EQ.2) THEN
          TL=TLM0
        ELSE
          STOP 'Error! [IC]'
        END IF
C        IF(NT.GT.1) READ(1,*) (AW(I),I=1,NT-1)
C       READ(1,*) (AWI(I),I=1,NT)
C        IF(NT.GT.1) READ(1,*) (ASL(I),I=2,NT)
!        WRITE(5,1) REM,TIT(IC),IC,NT,TH,TLM0,TF
!1       FORMAT(//20X,'INPUT  DATA',//A80,//2X,'[IC-NT-TH-',A3,
!     *  '-TF]',/2I5,3F15.5)
!        IF(NT.GT.1) WRITE(5,2) (AW(I),I=1,NT-1)
!2       FORMAT(/2X,'[AW(1:NT-1)]',/(6F13.5))
!        WRITE(5,3) (AWI(I),I=1,NT)
!3       FORMAT(/2X,'[AWI(1:NT)]',/,(6F13.5))
!        IF(NT.GT.1) WRITE(5,4) (ASL(I),I=2,NT)
!4       FORMAT(/2X,'[(ASL(2:NT)]',/,(6F13.5))
!        WRITE(5,5)
!5       FORMAT(//,80('='))
        IF(IC.EQ.2) THEN
          SL1=TL
          DO 10 I=2,NT
10        SL1=SL1-ASL(I)
          ASL(1)=SL1
        END IF
        RETURN
        END
C***********************************************************************
        SUBROUTINE ACAC0(IC,NT,TF,TN0,ASL,AW,AWI,AC,AC0)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER IC,NT
        REAL*8 TF,TN0
        DIMENSION ASL(10),AW(10),AWI(10),AC(10),AC0(10)
        IF(IC.EQ.1) THEN
          AC0(1)=TN0
          AC(1)=0
        ELSE
          AC0(1)=0
          AC(1)=AWI(1)*ASL(1)/TF
        END IF
        DO 10 I=2,NT
          AC0(I)=AC(I-1)+AW(I-1)/TF
          AC(I)=AC0(I)+AWI(I)*ASL(I)/TF
10      CONTINUE
        RETURN
        END
C***********************************************************************
        SUBROUTINE SEARCH(X0,DX,X,F,NL)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NL
        REAL*8 X0,DX,X,F,X1,F1,X2,F2,EPX,EPF
        PARAMETER (NMAX=200,EPX=1.E-8)
!        WRITE(5,1)
!1       FORMAT(//25X,'***** CHECK *****')
        X=X0
        NL=1
        F=FUNC_ALP(X)
        F1=F
10      X1=X
        F0=F1
        F1=F
20      X=X1+DX
        NL=NL+1
        IF(NL.GT.NMAX) STOP 'Error (data) !'
        F=FUNC_ALP(X)
        FF=F1*F
        IF(FF.GT.0.) GOTO 10
        IF(FF.LT.0..AND.(F1-F0)*(F1-F).GE.0.) GOTO 10
30      EPF=SQRT(ABS(F1*F))*EPX
        X2=X
        F2=F
        CALL SECANT(X1,F1,X2,F2,EPX,EPF,X,NL)
        RETURN
        END
C***********************************************************************
        SUBROUTINE SECANT(XA,FA,XB,FB,EPX,EPF,X,NL)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NL
        REAL*8 XA,FA,XB,FB,EPX,EPF,X
        PARAMETER (NMAX=200)
        X=XA
        IF(ABS(FA).LE.EPF) RETURN
        X=XB
        IF(ABS(FB).LE.EPF) RETURN
10      X=XB-FB*(XB-XA)/(FB-FA)
        NL=NL+1
        IF(NL.GT.NMAX) STOP '[Unsuitable step !]'
        F=FUNC_ALP(X)
        IF(ABS(F).LE.EPF) RETURN
        IF(FB*F.LT.0) THEN
          XA=XB
          FA=FB
        ELSE
          FA=FA*FB/(FB+F)
        END IF
        XB=X
        FB=F
        IF(ABS(XB-XA).GT.EPX) GOTO 10
        IF(ABS(FA).LT.ABS(FB)) X=XA
        RETURN
        END
C***********************************************************************
        FUNCTION FUNC_ALP(ALP)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NT
        REAL*8 TH,TF,ALP
        DIMENSION AWI(10),AC(10),AC0(10)
        REAL*8 MO,N0,N
        COMMON /FUN/IC,NT,TH,TF,AWI,AC,AC0
        MO=SQRT(ALP**2-1)-AC(NT)*ALP
        F=-AWI(1)*TH/TF*ALP
        DO 10 I=1,NT
        N=SQRT(1+(MO+AC(I)*ALP)**2)
        IF(IC.EQ.1.AND.I.EQ.1) THEN
          N0=AC0(1)
        ELSE
          N0=SQRT(1+(MO+AC0(I)*ALP)**2)
        END IF
        F=F+AWI(1)/AWI(I)*(N-N0)
10      CONTINUE
!        WRITE(5,1) ALP,F
!1       FORMAT(2X,'alp=',E20.8,5X,'F(alp)=',E20.8)
        FUNC_ALP=F
        RETURN
        END
C***********************************************************************
        SUBROUTINE PARAMETERS(NT,TFH,TM0,ASL,AW,AWI,ASH,ASS,
     *  ASM0,ASM,ASN0,ASN,TS)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NT
        REAL*8 TS,TFH,TM0
        DIMENSION AW(10),AWI(10),ASL(10),ASH(10),ASS(10),
     *  ASM0(10),ASM(10),ASN0(10),ASN(10)
        TS=0
        DO 10 I=1,NT
          IF(I.EQ.1) THEN
            SM0=TM0
          ELSE
            SM0=ASM(I-1)+AW(I-1)/TFH
          END IF
          SA=TFH/AWI(I)
          SM=SM0+ASL(I)/SA
          SN0=SQRT(1+SM0**2)
          SN=SQRT(1+SM**2)
          ASINH0=LOG(SM0+SQRT(SM0**2+1))
          ASINH=LOG(SM+SQRT(SM**2+1))
          ASM0(I)=SM0
          ASM(I)=SM
          ASN0(I)=SN0
          ASN(I)=SN
          ASH(I)=SA*(SN-SN0)
          ASS(I)=SA*(ASINH-ASINH0)
          TS=TS+ASS(I)
10      CONTINUE
        RETURN
        END
C***********************************************************************
        SUBROUTINE STIFFNESS(NT,TFH,TS,ASH,ASL,ASM0,ASN0,ASN,
     *  RXX,RXY,RYX,RYY)
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NT
        REAL*8 TS,TFH,RXX,RXY,RYX,RYY
        DIMENSION ASH(10),ASL(10),ASM0(10),ASN0(10),ASN(10)
        REAL*8 H,L,M0,N0,N
        A=0
        B=0
         DO 10 I=1,NT
           H=ASH(I)
          L=ASL(I)
          M0=ASM0(I) 
		N=ASN(I)
          N0=ASN0(I)
          A=A+(N0*L-M0*H)/N/N0
          B=B+H/N/N0
10      CONTINUE
        C=B
        D=TS-A
        E=A*D-B*C
        RXX=TFH*A/E
        RXY=TFH*B/E
        RYX=TFH*C/E
        RYY=TFH*D/E
        RETURN
	END
C******************************************************************C
C      THE END OF THIS  PROGRAMME                                  C
C******************************************************************C
