        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:12 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SECTMOTION__genmod
          INTERFACE 
            SUBROUTINE SECTMOTION(Y,T)
              USE CONSTANT, ONLY :                                      &
     &          NR,                                                     &
     &          COEFT
              REAL(KIND=8) :: Y(2*NR)
              REAL(KIND=8) :: T
            END SUBROUTINE SECTMOTION
          END INTERFACE 
        END MODULE SECTMOTION__genmod
