        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:01 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INST_EMESH__genmod
          INTERFACE 
            SUBROUTINE INST_EMESH(NWH,NL,NODE0,INSTUR,INSTDUR)
              USE CONSTANT, ONLY :                                      &
     &          NR
              INTEGER(KIND=4) :: NL
              INTEGER(KIND=4) :: NWH
              REAL(KIND=8) :: NODE0(NWH,NL,3)
              REAL(KIND=8) :: INSTUR(NWH,NL,NR,3)
              REAL(KIND=8) :: INSTDUR(NWH,NL,NR,3,3)
            END SUBROUTINE INST_EMESH
          END INTERFACE 
        END MODULE INST_EMESH__genmod
