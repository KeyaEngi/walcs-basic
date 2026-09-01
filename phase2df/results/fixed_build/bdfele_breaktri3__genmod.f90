        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE BDFELE_BREAKTRI3__genmod
          INTERFACE 
            SUBROUTINE BDFELE_BREAKTRI3(NODE0,ZZ,ELEEXIST,ELEKIND,      &
     &NEWNODE,NEWELEAREA,NEWXAV)
              REAL(KIND=8) :: NODE0(3,3)
              REAL(KIND=8) :: ZZ(3)
              INTEGER(KIND=4) :: ELEEXIST
              INTEGER(KIND=4) :: ELEKIND
              REAL(KIND=8) :: NEWNODE(4,3)
              REAL(KIND=8) :: NEWELEAREA
              REAL(KIND=8) :: NEWXAV(3)
            END SUBROUTINE BDFELE_BREAKTRI3
          END INTERFACE 
        END MODULE BDFELE_BREAKTRI3__genmod
