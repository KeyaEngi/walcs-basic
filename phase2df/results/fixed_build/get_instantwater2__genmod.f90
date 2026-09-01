        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_INSTANTWATER2__genmod
          INTERFACE 
            SUBROUTINE GET_INSTANTWATER2(NWH,NL,NODE0,INST_FETA,        &
     &BREAKNODE0,BREAKINST_FETA,NODE1,SUBMERGE)
              INTEGER(KIND=4) :: NL
              INTEGER(KIND=4) :: NWH
              REAL(KIND=8) :: NODE0(NWH,NL,3)
              REAL(KIND=8) :: INST_FETA(NWH,NL)
              REAL(KIND=8) :: BREAKNODE0(NL,3)
              REAL(KIND=8) :: BREAKINST_FETA(NL)
              REAL(KIND=8) :: NODE1(NWH,NL,3)
              INTEGER(KIND=4) :: SUBMERGE
            END SUBROUTINE GET_INSTANTWATER2
          END INTERFACE 
        END MODULE GET_INSTANTWATER2__genmod
