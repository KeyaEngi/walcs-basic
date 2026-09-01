        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:28 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_WAGCONDITION__genmod
          INTERFACE 
            SUBROUTINE SLAM_WAGCONDITION(NUM0,POINT0,SDX,NUMINTC,INTC,  &
     &INTNC,DERNC,TEMRISE)
              INTEGER(KIND=4) :: NUMINTC
              INTEGER(KIND=4) :: NUM0
              REAL(KIND=8) :: POINT0(1:NUM0,1:2)
              REAL(KIND=8) :: SDX(1:NUM0)
              REAL(KIND=8) :: INTC(NUMINTC)
              REAL(KIND=8) :: INTNC(NUMINTC)
              REAL(KIND=8) :: DERNC(NUMINTC+1)
              REAL(KIND=8) :: TEMRISE(NUMINTC)
            END SUBROUTINE SLAM_WAGCONDITION
          END INTERFACE 
        END MODULE SLAM_WAGCONDITION__genmod
