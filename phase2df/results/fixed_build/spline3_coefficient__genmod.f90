        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:31 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE3_COEFFICIENT__genmod
          INTERFACE 
            SUBROUTINE SPLINE3_COEFFICIENT(NUM,POINT,DS1,DSN,AA,D)
              INTEGER(KIND=4) :: NUM
              REAL(KIND=8) :: POINT(NUM,2)
              REAL(KIND=8) :: DS1
              REAL(KIND=8) :: DSN
              REAL(KIND=8) :: AA(NUM,NUM)
              REAL(KIND=8) :: D(NUM)
            END SUBROUTINE SPLINE3_COEFFICIENT
          END INTERFACE 
        END MODULE SPLINE3_COEFFICIENT__genmod
