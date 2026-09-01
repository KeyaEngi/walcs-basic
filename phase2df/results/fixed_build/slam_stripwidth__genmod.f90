        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:27 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_STRIPWIDTH__genmod
          INTERFACE 
            SUBROUTINE SLAM_STRIPWIDTH(NUMLINE,TEMSTRIPWIDTH,           &
     &TEMSTRIPBNODE)
              USE SLAMMING, ONLY :                                      &
     &          SLAMNUMLINE,                                            &
     &          SLAMANGLE,                                              &
     &          SLAMSTRIPBEXIST,                                        &
     &          SLAMSTRIPBANGLE,                                        &
     &          SLAMNUMLP,                                              &
     &          SLAMNODE,                                               &
     &          SLAMINTNUMP,                                            &
     &          SLAMSTRIPTYP,                                           &
     &          SLAMWIDTHMCOEF
              INTEGER(KIND=4) :: NUMLINE
              REAL(KIND=8) :: TEMSTRIPWIDTH(SLAMINTNUMP)
              REAL(KIND=8) :: TEMSTRIPBNODE(SLAMINTNUMP,2,3)
            END SUBROUTINE SLAM_STRIPWIDTH
          END INTERFACE 
        END MODULE SLAM_STRIPWIDTH__genmod
