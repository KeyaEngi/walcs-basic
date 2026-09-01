        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:26 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_PRESPACIALINTEGRAL__genmod
          INTERFACE 
            SUBROUTINE SLAM_PRESPACIALINTEGRAL(DIS,NUMLINE,NUMPRENODE,  &
     &PRENODE,SLAMPRE,SFORCE,SLINEFORCE)
              USE CONSTANT, ONLY :                                      &
     &          ROU,                                                    &
     &          NR
              USE SLAMMING, ONLY :                                      &
     &          SLAMINTNUMP,                                            &
     &          SLAMNUMLP,                                              &
     &          SLAMPOINT,                                              &
     &          SLAMNUMINTEGRATION,                                     &
     &          SLAMUR,                                                 &
     &          SLAMRUR,                                                &
     &          SLAMANGLE,                                              &
     &          SLAMWIDTH,                                              &
     &          SLAMNODE
              INTEGER(KIND=4) :: NUMPRENODE
              REAL(KIND=8) :: DIS(6)
              INTEGER(KIND=4) :: NUMLINE
              REAL(KIND=8) :: PRENODE(0:NUMPRENODE,1:2)
              REAL(KIND=8) :: SLAMPRE(0:NUMPRENODE)
              REAL(KIND=8) :: SFORCE(NR)
              REAL(KIND=8) :: SLINEFORCE(3)
            END SUBROUTINE SLAM_PRESPACIALINTEGRAL
          END INTERFACE 
        END MODULE SLAM_PRESPACIALINTEGRAL__genmod
