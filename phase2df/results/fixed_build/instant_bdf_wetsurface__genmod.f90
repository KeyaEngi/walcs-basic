        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:02 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INSTANT_BDF_WETSURFACE__genmod
          INTERFACE 
            SUBROUTINE INSTANT_BDF_WETSURFACE(Y,IT,T,SMTF,INSTFORCE)
              USE CONSTANT
              USE PANELGEOMETRY
              REAL(KIND=8) :: Y(1:2*NR)
              INTEGER(KIND=4) :: IT
              REAL(KIND=8) :: T
              REAL(KIND=8) :: SMTF
              REAL(KIND=8) :: INSTFORCE(NR)
            END SUBROUTINE INSTANT_BDF_WETSURFACE
          END INTERFACE 
        END MODULE INSTANT_BDF_WETSURFACE__genmod
