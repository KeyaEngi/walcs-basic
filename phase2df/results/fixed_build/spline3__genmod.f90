        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:31 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE3__genmod
          INTERFACE 
            SUBROUTINE SPLINE3(NUM,POINT,DS1,DSN,SSS)
              INTEGER(KIND=4) :: NUM
              REAL(KIND=8) :: POINT(NUM,2)
              REAL(KIND=8) :: DS1
              REAL(KIND=8) :: DSN
              REAL(KIND=8) :: SSS(NUM-1,4)
            END SUBROUTINE SPLINE3
          END INTERFACE 
        END MODULE SPLINE3__genmod
