        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:23 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_CURVEOPTIMIZE__genmod
          INTERFACE 
            SUBROUTINE SLAM_CURVEOPTIMIZE(NUM0,POINT0,M,NUMBOUNPOINT,   &
     &BPOINT2,BDX,NUMINCLINE,INCLINEB)
              INTEGER(KIND=4) :: M
              INTEGER(KIND=4) :: NUM0
              REAL(KIND=8) :: POINT0(NUM0,2)
              INTEGER(KIND=4) :: NUMBOUNPOINT
              REAL(KIND=8) :: BPOINT2(M,2)
              REAL(KIND=8) :: BDX(M)
              INTEGER(KIND=4) :: NUMINCLINE
              REAL(KIND=8) :: INCLINEB(100,2)
            END SUBROUTINE SLAM_CURVEOPTIMIZE
          END INTERFACE 
        END MODULE SLAM_CURVEOPTIMIZE__genmod
