        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAMCASE3__genmod
          INTERFACE 
            SUBROUTINE SLAMCASE3(SMTF,NUMT,TIME0,TIME1,DIS0,DIS1,       &
     &MESECLOAD0,MESECLOAD1,T0SLAMFORCE,MAINFORCE)
              USE CONSTANT, ONLY :                                      &
     &          PI,                                                     &
     &          U0,                                                     &
     &          ROU,                                                    &
     &          NR
              USE SLAMMING
              INTEGER(KIND=4) :: NUMT
              REAL(KIND=8) :: SMTF
              REAL(KIND=8) :: TIME0
              REAL(KIND=8) :: TIME1
              REAL(KIND=8) :: DIS0(6)
              REAL(KIND=8) :: DIS1(6)
              REAL(KIND=8) :: MESECLOAD0(1:6,1:SLAMNUMLINE)
              REAL(KIND=8) :: MESECLOAD1(1:6,1:SLAMNUMLINE)
              REAL(KIND=8) :: T0SLAMFORCE(NR)
              REAL(KIND=8) :: MAINFORCE(NR)
            END SUBROUTINE SLAMCASE3
          END INTERFACE 
        END MODULE SLAMCASE3__genmod
