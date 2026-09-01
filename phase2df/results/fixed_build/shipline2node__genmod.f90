        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:17 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SHIPLINE2NODE__genmod
          INTERFACE 
            SUBROUTINE SHIPLINE2NODE(NUMLINE,NUMPOINT,POINT,TF,TA)
              USE PANELGEOMETRY, ONLY :                                 &
     &          NWH,                                                    &
     &          NL,                                                     &
     &          NODE,                                                   &
     &          NODEB,                                                  &
     &          INSTBREAKNODE,                                          &
     &          BREAKKEXI
              INTEGER(KIND=4) :: NUMLINE
              INTEGER(KIND=4) :: NUMPOINT(NUMLINE)
              REAL(KIND=8) :: POINT(NUMLINE,100,3)
              REAL(KIND=8) :: TF
              REAL(KIND=8) :: TA
            END SUBROUTINE SHIPLINE2NODE
          END INTERFACE 
        END MODULE SHIPLINE2NODE__genmod
