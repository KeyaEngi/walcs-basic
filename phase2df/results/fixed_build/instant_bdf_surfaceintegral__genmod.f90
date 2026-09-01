        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:02 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INSTANT_BDF_SURFACEINTEGRAL__genmod
          INTERFACE 
            SUBROUTINE INSTANT_BDF_SURFACEINTEGRAL(NUMN,ELEKIND,ELE_PIS,&
     &ELE_NODE,ELE_UR,FORCE)
              USE CONSTANT, ONLY :                                      &
     &          ROU,                                                    &
     &          NR
              INTEGER(KIND=4) :: NUMN
              INTEGER(KIND=4) :: ELEKIND(NUMN)
              REAL(KIND=8) :: ELE_PIS(NUMN,4)
              REAL(KIND=8) :: ELE_NODE(NUMN,4,3)
              REAL(KIND=8) :: ELE_UR(NUMN,4,NR,3)
              REAL(KIND=8) :: FORCE(1:NR)
            END SUBROUTINE INSTANT_BDF_SURFACEINTEGRAL
          END INTERFACE 
        END MODULE INSTANT_BDF_SURFACEINTEGRAL__genmod
