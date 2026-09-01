        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GAUSS_JORDAN__genmod
          INTERFACE 
            SUBROUTINE GAUSS_JORDAN(ROW,COL,MATRIXINI,MATRIX1)
              INTEGER(KIND=4) :: COL
              INTEGER(KIND=4) :: ROW
              REAL(KIND=8) :: MATRIXINI(ROW,COL)
              REAL(KIND=8) :: MATRIX1(ROW,COL)
            END SUBROUTINE GAUSS_JORDAN
          END INTERFACE 
        END MODULE GAUSS_JORDAN__genmod
