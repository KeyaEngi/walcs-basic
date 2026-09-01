        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:23 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_COARSECURVE__genmod
          INTERFACE 
            SUBROUTINE SLAM_COARSECURVE(SECID,TEMSLAMCASE,TYPEANGLE,    &
     &TEMNUMSECTP2,TEMSECTP2)
              USE SLAMMING, ONLY :                                      &
     &          SLAMLIBNUMZ,                                            &
     &          SLAMLIBNUMPORT,                                         &
     &          SLAMLIBPORTNODE,                                        &
     &          SLAMLIBTYPE
              INTEGER(KIND=4) :: TEMNUMSECTP2
              INTEGER(KIND=4) :: SECID
              INTEGER(KIND=4) :: TEMSLAMCASE
              REAL(KIND=8) :: TYPEANGLE
              REAL(KIND=8) :: TEMSECTP2(TEMNUMSECTP2,3)
            END SUBROUTINE SLAM_COARSECURVE
          END INTERFACE 
        END MODULE SLAM_COARSECURVE__genmod
