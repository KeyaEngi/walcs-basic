        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:31 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SURFACEINTEGRAL1__genmod
          INTERFACE 
            SUBROUTINE SURFACEINTEGRAL1(NWH,NL,NODE0,DTBPHI,FORCE,INSTUR&
     &)
              USE CONSTANT, ONLY :                                      &
     &          ROU,                                                    &
     &          NR
              INTEGER(KIND=4) :: NL
              INTEGER(KIND=4) :: NWH
              REAL(KIND=8) :: NODE0(1:NWH,1:NL,1:3)
              REAL(KIND=8) :: DTBPHI(1:NWH,1:NL)
              REAL(KIND=8) :: FORCE(1:NR)
              REAL(KIND=8) :: INSTUR(NWH,NL,NR,3)
            END SUBROUTINE SURFACEINTEGRAL1
          END INTERFACE 
        END MODULE SURFACEINTEGRAL1__genmod
