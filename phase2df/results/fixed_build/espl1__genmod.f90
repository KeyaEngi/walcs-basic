        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ESPL1__genmod
          INTERFACE 
            SUBROUTINE ESPL1(X,Y,N,DY1,DYN,XX,M,S,DS,DDS)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: X(N)
              REAL(KIND=8) :: Y(N)
              REAL(KIND=8) :: DY1
              REAL(KIND=8) :: DYN
              REAL(KIND=8) :: XX(M)
              REAL(KIND=8) :: S(M)
              REAL(KIND=8) :: DS(M)
              REAL(KIND=8) :: DDS(M)
            END SUBROUTINE ESPL1
          END INTERFACE 
        END MODULE ESPL1__genmod
