        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:48:58 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EMESH__genmod
          INTERFACE 
            SUBROUTINE EMESH(NWH,NL,NODE0)
              USE CONSTANT, ONLY :                                      &
     &          NR
              USE SHIPHULLVAR, ONLY :                                   &
     &          STRUSECT,                                               &
     &          STRUXN,                                                 &
     &          DRM,                                                    &
     &          STRUZN,                                                 &
     &          STRUYN,                                                 &
     &          STRUZSC,                                                &
     &          UR,                                                     &
     &          DUR,                                                    &
     &          POINTCOR,                                               &
     &          NITEM,                                                  &
     &          POINTCORE,                                              &
     &          POINTCORYN,                                             &
     &          POINTCORZN,                                             &
     &          POINTCORZSC
              INTEGER(KIND=4) :: NL
              INTEGER(KIND=4) :: NWH
              REAL(KIND=8) :: NODE0(NWH,NL,3)
            END SUBROUTINE EMESH
          END INTERFACE 
        END MODULE EMESH__genmod
