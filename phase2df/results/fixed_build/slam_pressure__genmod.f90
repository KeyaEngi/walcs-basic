        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:27 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_PRESSURE__genmod
          INTERFACE 
            SUBROUTINE SLAM_PRESSURE(NUMLINE,INIPENETRATION,INIRISE,HT, &
     &VEL,ACE,NUMPRENODE,PRENODE,SLAMPRE)
              USE SLAMMING, ONLY :                                      &
     &          SLAMINTNUMP,                                            &
     &          SLAMNUMLP,                                              &
     &          SLAMPOINT,                                              &
     &          SLAMDX,                                                 &
     &          SLAMNUMINTC,                                            &
     &          SLAMNUMINTEGRATION,                                     &
     &          SLAMINTC,                                               &
     &          SLAMNC,                                                 &
     &          SLAMRISEC,                                              &
     &          SLAMDERNC,                                              &
     &          SLAMNUMABDENT,                                          &
     &          SLAMABDENTBZ,                                           &
     &          NUMEXPSLAMCP,                                           &
     &          EXPSLAMCP
              INTEGER(KIND=4) :: NUMLINE
              REAL(KIND=8) :: INIPENETRATION
              REAL(KIND=8) :: INIRISE
              REAL(KIND=8) :: HT
              REAL(KIND=8) :: VEL
              REAL(KIND=8) :: ACE
              INTEGER(KIND=4) :: NUMPRENODE
              REAL(KIND=8) :: PRENODE(0:SLAMNUMINTEGRATION,1:2)
              REAL(KIND=8) :: SLAMPRE(0:SLAMNUMINTEGRATION)
            END SUBROUTINE SLAM_PRESSURE
          END INTERFACE 
        END MODULE SLAM_PRESSURE__genmod
