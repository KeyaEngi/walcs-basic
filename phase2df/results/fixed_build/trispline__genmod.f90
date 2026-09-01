        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:34 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TRISPLINE__genmod
          INTERFACE 
            SUBROUTINE TRISPLINE(N,POINT1,POINT2,DS1,DSN,NUMN,INTERX,NS,&
     &NDS,NDDS)
              INTEGER(KIND=4) :: NUMN
              INTEGER(KIND=4) :: N
              REAL(KIND=8) :: POINT1(N,1)
              REAL(KIND=8) :: POINT2(N,1)
              REAL(KIND=8) :: DS1
              REAL(KIND=8) :: DSN
              REAL(KIND=8) :: INTERX(NUMN)
              REAL(KIND=8) :: NS(NUMN)
              REAL(KIND=8) :: NDS(NUMN)
              REAL(KIND=8) :: NDDS(NUMN)
            END SUBROUTINE TRISPLINE
          END INTERFACE 
        END MODULE TRISPLINE__genmod
