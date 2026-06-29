      MODULE GMRES_MOD

	PRIVATE
	PUBLIC::BGMRES

	CONTAINS
!    N - 方程组阶数，M-重启动数（取50至100间），NS-方程组右端项个数，（同一系数矩阵，对辐射问题取6）
!    A-系数矩阵，Q-解，B-右端项
      SUBROUTINE BGMRES(A,BB,Q,N,NS,M,ITR)
!      USE BGMRES_ARRAYS
	PARAMETER( ITERMAX=200,TOL=1.0E-5)
	PARAMETER (ZERO = (0.0, 0.0) )
C	COMPLEX   BB(N,NS)
C	COMPLEX BV(N,NS,M)
	COMPLEX BH(M+1,M)
	
	COMPLEX  G(M+1),Y(M+1)
	COMPLEX S(M),   HIJ,HIJ1 ,A(N,N), RNORM1,BNORM1
	COMPLEX  BB(N,*),BV(N,NS,M),Q(N,*),P(N,NS),R(N,NS),AP(N,NS) 
	
       REAL C(M)
!           A(N,N)- COEFFICIENT MATRIX,BH(M+1,M+1)-HESSENBERG MATRIX
!           BB(N,NS)-RIGHT  HAND SIDE, P(N,NS),AP(N,NS)-INTERMEDIATE WORKSPACE
!            Q(N,NS) -LINERA SYSTEM SOLUTION, R(N)-RESIDUAL
!             Y(M,NS)-LEAST SQUARE  SOLUTION 
!           C(M,NS),S(M,NS)-COSINE ,SINE OF THE GIVENS ROTATION
!           ITR-ACTUAL ITERAION NUMBER
        

	TOL1=1
	ITEROUT=0
	KSTART=0
!	OPEN(151,FILE='IHISTORY.DAT')
!	OPEN(66,FILE='RESULT.DAT')
!	WRITE(6,*) 'BGMRES'
C      DO J = 1,N
C        DO I = 1,N
C          A(I,J) = ZERO
C        ENDDO
C        WORK(J,1) = ONEC
C	WORK(J,2)=TWO
C      ENDDO

   



C      DO I = 1,N
C        A(I,I) = (4.0, 0.0)
C	  BB(I,1)=ZERO
C	  BB(I,2)=ZERO
C      ENDDO
C      DO I = 1,N-1
C        A(I,I+1) = (-2.0, 1.0)
C        A(I+1,I) = (-1.0, 1.0)
C      ENDDO
*      
C      DO J=1,NS
C       DO I=1,N
	  
C         DO K=1,N
C         BB(I,J)=BB(I,J)+A(I,K)*WORK(K,J)
C         ENDDO
C	  ENDDO
C	ENDDO
 
      CALL FNORM(BB,BB,N,NS,BNORM1)
	BNORM=REAL(BNORM1)
	BNORM=SQRT(BNORM)

	DO  I=1,N
	   DO J=1,NS
          Q(I,J)=ZERO
	   ENDDO
	ENDDO
   
       CALL EQUAL(R,BB,N,NS)
  88   CONTINUE	   
!   SET UP V^1 AND G^0
	CALL FNORM(R,R,N,NS,RNORM1)
C	WRITE(6,*) RNORM1
	RNORM=REAL(RNORM1)
	RNORM=SQRT(RNORM)
	IF(RNoRM<=1.0E-9)then
	P=0
	else
	  DO I=1,N
	     DO J=1,NS
           P(I,J)=R(I,J)/RNORM
	     ENDDO
	  ENDDO
	  end if

       DO  I=1,M+1
	  G(I)=ZERO
	 ENDDO
       G(1)=CMPLX(RNORM,0.0)
!      BEGIN ITERATION
       ITER=1
!	WRITE(66,*) 'BETA',G(1)
       DO WHILE((ITER.LE.M).AND.(TOL1.GT.TOL).AND.(ITEROUT*M+ITER.LE
     $	 .ITERMAX))

! SAVE  P AS THE V(ITER)
        CALL EQUAL(BV(1,1,ITER),P,N,NS)
 !       WRITE(66,*)'ITER',ITER
!	  WRITE(66,*) ((BV(I1,J1,ITER),I1=1,N),J1=1,NS)

!     FORM  AV(ITER)	   
	  CALL  MATAV(A,P,AP,N,NS)
!	IF(ITER.EQ.1) WRITE(*,*) 'AP', AP(3,1)
!  INITIALIZE V^(ITER+1) WITH AV(ITER)
	  CALL EQUAL(P,AP,N,NS)

! MAKE V^(ITER+1) ORTHOGONAL TO V^(I),I<=ITER
       DO  J=1,ITER
	   CALL FNORM(BV(1,1,J),AP,N,NS,HIJ)
 
	   DO  I=1,N
	     DO K=1,NS
		 P(I,K)=P(I,K)-HIJ*BV(I,K,J)
		 ENDDO
	   ENDDO 
	   BH(J,ITER)=HIJ
	 ENDDO 
!NORMALIZE V^(ITER+1)
	  CALL FNORM(P,P,N,NS,RNORM1)
	RNORM =REAL(RNORM1)
 
	  RNORM =SQRT(RNORM)
	IF(RNoRM<=1.0E-9)then
	DO  I=1,N
	   DO J=1,NS
          Q(I,J)=ZERO
	   ENDDO
	ENDDO
       return
      end if
	  DO  I=1,N
	     DO J=1,NS
		 P(I,J)=P(I,J)/RNORM
		 ENDDO
	  ENDDO
	BH(ITER+1,ITER)=CMPLX(RNORM,0.0) 
!	 WRITE(66,*)'H OLD BEFORE ', (BH(I,ITER),I=1,ITER+1)    

! APPLY ROTATIONS TO NEW H COLUMN
       DO I=1,ITER-1
	   HIJ=BH(I,ITER)
	   HIJ1=BH(I+1,ITER)
	   BH(I,ITER)=HIJ*C(I) +HIJ1*CONJG(S(I))
	   BH(I+1,ITER)=-HIJ*S(I)+HIJ1*C(I)
	 ENDDO
!	WRITE(66,*)'H OLD AFTER ', (BH(I,ITER),I=1,ITER)  
!  COMPUTE NEW ROTATIONS
        HIJ=BH(ITER,ITER)
	  HIJ1=BH(ITER+1,ITER) 	   
        CALL GIVENS(HIJ,HIJ1,C(ITER),S(ITER))
!        WRITE(66,*) 'GIVENS',C(ITER),S(ITER)
!  APPLY NEW ROTATIONS
       BH(ITER,ITER)= HIJ*C(ITER)+HIJ1*CONJG(S(ITER))
       BH(ITER+1,ITER)=ZERO
!       WRITE(66,*)'H NEW AFTER ', (BH(I,ITER),I=1,ITER+1) 

        HIJ=G(ITER)
       G(ITER)=C(ITER)*HIJ
	G(ITER+1)=-S(ITER)*HIJ
       TOL1=CABS(G(ITER+1))/BNORM
!       WRITE(66,*)'G  AFTER ', (G(I ),I=1,ITER+1) 
C      TOL1=CABS(G(ITER+1))
!      WRITE(6,*) 'ITER',ITER
!	WRITE(6,*) 'TOL1',TOL1
!	WRITE(66,*) 'TOLERANCE',TOL1
!     WRITE(151,*) ITEROUT*M+ITER,LOG10(TOL1)
	ITER=ITER+1
  
	ENDDO
 
       ITER=ITER-1
C       WRITE(6,*) 'TOL1',TOL1
 
!     COMPUTE  LEAST SQUARE  SOLUTION
       DO I=1,ITER
	  Y(I)=G(I)
	 ENDDO
!	WRITE(66,*)'G'
!	WRITE(66,*) (G(I),I=1,ITER)
	  DO I=ITER ,1,-1 
 	     Y(I)=Y(I)/BH(I,I)   
 	     DO J=I-1,1,-1
            Y(J)=Y(J)-BH(J,I)*Y(I)
           ENDDO
       ENDDO
!	WRITE(66,*) 'Y'
!	WRITE(66,*)(Y(J),J=1,ITER)
!	 WRITE(66,*)'BV'
	 DO  K=1,ITER
!	 WRITE(66,*)'ITER',K
	
!       WRITE(66,*)(BV(I,1,K),I=1,N)
	ENDDO

!      ACQUIRE ITERATION SOLUTION
	  DO I=1,N
	   DO J=1,NS
	    DO K=1,ITER
         Q(I,J)=Q(I,J)+Y(K)*BV(I,J,K)
	   ENDDO
	  ENDDO
	 ENDDO
!	WRITE(66,*) 'Q'
!      WRITE(66,*) (Q(I,1),I=1,N)

      IF((TOL1.LT.TOL).AND.(ITEROUT*M+ITER.LE.ITERMAX)) THEN
!      ITERATIONS END
	 GOTO 98
	ENDIF
      KSTART=0
	IF((TOL1.GT.TOL).AND.(ITEROUT*M+ITER.LT.ITERMAX)) KSTART=1

       IF( (ITEROUT*M+ITER.GT.ITERMAX).AND.(TOL1.GT.TOL) ) THEN
!	   WRITE(6,*) 'WARNING NO CONVERGENCE  AFTER ' 
!	   WRITE(6,*) ITEROUT*M+ITER,'ITERATIONS'
!	   WRITE(6,*) 'TOLRENCE  IS',TOL1
C	   WRITE(6,*) 'ITEROUT',ITEROUT
          GOTO  109
	 ENDIF

	IF(KSTART.EQ.1) THEN

	ITEROUT=ITEROUT+1
C  　    
!       WRITE(6,*) 'ITEROUT',ITEROUT                                                                                                                                                                                                                                                       
       CALL  MATAV(A,Q,AP,N,NS)
       DO  I=1,N
	    DO J=1,NS
	     R(I,J)=BB(I,J)-AP(I,J)
	    ENDDO
	ENDDO 
!     ITERATIONS BEGIN AGAIN
      
	GOTO 88
      ENDIF
   98  CONTINUE
!         WRITE(6,*) 'CONGRATULATIONS CONVERGENCE ACQUIRE ' 
!	   WRITE(6,*)    'AFTER', ITEROUT*M+ITER,' ITERATIONS'
!	   WRITE(6,*) 'TOLRENCE  IS',TOL1  

  99    CONTINUE 
       CALL EQUAL(BB,Q,N,NS)
	ITR=ITEROUT*M+ITER
	CLOSE(66)
      



C	DO J=1,NS  
C	DO I=1,N
 
C  	WRITE(*,*)BB(I,J)
C	ENDDO
C	ENDDO
!      DEALLOCATE(BV,Q,P,R,AP)
  109  CONTINUE
C       WRITE(*,*)"USED"
	RETURN
	END  SUBROUTINE BGMRES
	  

	 
	 	  
 


	SUBROUTINE MATAV(A,V,AP,N,NS)
	COMPLEX A(N,N),V(N,NS),AP(N,NS)
	PARAMETER (ZERO = (0.0, 0.0), ONE = (1.0, 0.0))
	DO  J=1,NS
	  DO I=1,N
	   AP(I,J)=ZERO
	   ENDDO
	ENDDO
	 DO I=1,N   
	   DO  J=1,NS
	     DO  K=1,N
	     AP(I,J)=AP(I,J)+A(I,K)*V(K,J)
	     ENDDO
	   ENDDO
	ENDDO
	RETURN
	END SUBROUTINE MATAV


	SUBROUTINE FNORM(V,AP,N,NS,HIJ)
!  CALCULATE FROBENIUS INNER PRODUCT 
	COMPLEX V(N,NS),AP(N,NS),HIJ
	PARAMETER (ZERO = (0.0, 0.0), ONE = (1.0, 0.0))

	 HIJ=ZERO
	DO I=1,NS
	    DO K=1,N
	        HIJ=HIJ+CONJG(V(K,I))*AP(K,I)
	    ENDDO
	ENDDO
	 
	RETURN
	END  SUBROUTINE FNORM
	
	SUBROUTINE EQUAL(A,B,N,NS)
	COMPLEX A(N,NS),B(N,NS)
	DO I=1,N
	  DO J=1,NS
	  A(I,J)=B(I,J)
	  ENDDO
	ENDDO
      RETURN
	END  SUBROUTINE EQUAL
      
		     
       SUBROUTINE GIVENS(Z1,Z2,C,S)
*
       COMPLEX  Z1,Z2,S
       REAL      C
*
       REAL  NORMZ
       INTRINSIC CABS , SQRT 
*
       COMPLEX  ZERO, ONE
       PARAMETER (ZERO = (0.0E0, 0.0E0) , ONE = (1.0E0, 0.0E0) )
       REAL  DZRO,DONE
       PARAMETER (DZRO = 0.0E0 , DONE = 1.0E0 ) 
*
       NORMZ = SQRT ( CABS (Z1)**2 + CABS (Z2)**2)
       IF (CABS (Z1).NE.DZRO) THEN
         S = Z2* CABS (Z1)/(Z1*NORMZ)
         C = CABS (Z1)/NORMZ
       ELSE IF (CABS (Z2).NE.DZRO) THEN
         C = DZRO
         S = Z2/ CABS (Z2)
       ELSE
         C = DONE
         S = ZERO
       ENDIF
*
       RETURN
       END   SUBROUTINE GIVENS     
      END MODULE GMRES_MOD
