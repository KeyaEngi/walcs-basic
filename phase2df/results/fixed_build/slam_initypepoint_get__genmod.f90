        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:25 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_INITYPEPOINT_GET__genmod
          INTERFACE 
            SUBROUTINE SLAM_INITYPEPOINT_GET(NUM0,POINTX0,NUMLP,LPOINT, &
     &NUMANGAREA,ANGBOUNDARY,ANG_VALUE,NUMSWCOF,SW_BOUNDARYX,SW_COEF)
              INTEGER(KIND=4) :: NUMSWCOF
              INTEGER(KIND=4) :: NUMANGAREA
              INTEGER(KIND=4) :: NUMLP
              INTEGER(KIND=4) :: NUM0
              REAL(KIND=8) :: POINTX0(NUM0)
              REAL(KIND=8) :: LPOINT(NUMLP,3)
              REAL(KIND=8) :: ANGBOUNDARY(NUMANGAREA)
              REAL(KIND=8) :: ANG_VALUE(NUMANGAREA)
              REAL(KIND=8) :: SW_BOUNDARYX(NUMSWCOF,2)
              REAL(KIND=8) :: SW_COEF(NUMSWCOF)
            END SUBROUTINE SLAM_INITYPEPOINT_GET
          END INTERFACE 
        END MODULE SLAM_INITYPEPOINT_GET__genmod
