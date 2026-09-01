        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:48:57 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BDFEMESH__genmod
          INTERFACE 
            SUBROUTINE BDFEMESH(TEMNUM_ELE,XAV0,BDFELE_UR)
              USE CONSTANT, ONLY :                                      &
     &          NR
              INTEGER(KIND=4) :: TEMNUM_ELE
              REAL(KIND=8) :: XAV0(TEMNUM_ELE,3)
              REAL(KIND=8) :: BDFELE_UR(TEMNUM_ELE,NR,3)
            END SUBROUTINE BDFEMESH
          END INTERFACE 
        END MODULE BDFEMESH__genmod
