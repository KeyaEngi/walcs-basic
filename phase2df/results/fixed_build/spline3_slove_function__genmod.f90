        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:31 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE3_SLOVE_FUNCTION__genmod
          INTERFACE 
            SUBROUTINE SPLINE3_SLOVE_FUNCTION(NUM,POINT,MM,SS)
              INTEGER(KIND=4) :: NUM
              REAL(KIND=8) :: POINT(NUM,2)
              REAL(KIND=8) :: MM(NUM)
              REAL(KIND=8) :: SS(NUM-1,4)
            END SUBROUTINE SPLINE3_SLOVE_FUNCTION
          END INTERFACE 
        END MODULE SPLINE3_SLOVE_FUNCTION__genmod
