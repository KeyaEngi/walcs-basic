        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:48:57 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELASTIC_SURFACEINTEGRAL__genmod
          INTERFACE 
            SUBROUTINE ELASTIC_SURFACEINTEGRAL(DTBPHI,FORCE)
              USE PANELGEOMETRY
              USE CONSTANT, ONLY :                                      &
     &          ROU,                                                    &
     &          NR
              REAL(KIND=8) :: DTBPHI(1:NWH,1:NL)
              REAL(KIND=8) :: FORCE(1:NR)
            END SUBROUTINE ELASTIC_SURFACEINTEGRAL
          END INTERFACE 
        END MODULE ELASTIC_SURFACEINTEGRAL__genmod
