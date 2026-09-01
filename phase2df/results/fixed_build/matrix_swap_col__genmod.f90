        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE MATRIX_SWAP_COL__genmod
          INTERFACE 
            SUBROUTINE MATRIX_SWAP_COL(ROW,COL,MATRIX,K1,K2)
              INTEGER(KIND=4) :: COL
              INTEGER(KIND=4) :: ROW
              REAL(KIND=8) :: MATRIX(ROW,COL)
              INTEGER(KIND=4) :: K1
              INTEGER(KIND=4) :: K2
            END SUBROUTINE MATRIX_SWAP_COL
          END INTERFACE 
        END MODULE MATRIX_SWAP_COL__genmod
