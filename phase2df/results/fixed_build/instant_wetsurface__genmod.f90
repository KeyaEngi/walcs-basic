        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:03 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INSTANT_WETSURFACE__genmod
          INTERFACE 
            SUBROUTINE INSTANT_WETSURFACE(Y,IT,T,SMTF,FORCEIS)
              USE CONSTANT
              USE PANELGEOMETRY
              REAL(KIND=8) :: Y(1:2*NR)
              INTEGER(KIND=4) :: IT
              REAL(KIND=8) :: T
              REAL(KIND=8) :: SMTF
              REAL(KIND=8) :: FORCEIS(NR)
            END SUBROUTINE INSTANT_WETSURFACE
          END INTERFACE 
        END MODULE INSTANT_WETSURFACE__genmod
