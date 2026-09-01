        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:12 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SECTIONLOAD__genmod
          INTERFACE 
            SUBROUTINE SECTIONLOAD(Y,T,SMTF,IT)
              USE SHIPHULLVAR, ONLY :                                   &
     &          ELOADR,                                                 &
     &          NBSECT
              USE CONSTANT, ONLY :                                      &
     &          NR,                                                     &
     &          NRAMP
              REAL(KIND=8) :: Y(2*NR)
              REAL(KIND=8) :: T
              REAL(KIND=8) :: SMTF
              INTEGER(KIND=4) :: IT
            END SUBROUTINE SECTIONLOAD
          END INTERFACE 
        END MODULE SECTIONLOAD__genmod
