        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_BDF_INSTANTWATER__genmod
          INTERFACE 
            SUBROUTINE GET_BDF_INSTANTWATER(NUM0,NODE0,INSFATE0,        &
     &NUMNEW_ELE,NEWKIND,NEWNODE)
              INTEGER(KIND=4) :: NUM0
              REAL(KIND=8) :: NODE0(4,3)
              REAL(KIND=8) :: INSFATE0(4)
              INTEGER(KIND=4) :: NUMNEW_ELE
              INTEGER(KIND=4) :: NEWKIND(2)
              REAL(KIND=8) :: NEWNODE(2,4,3)
            END SUBROUTINE GET_BDF_INSTANTWATER
          END INTERFACE 
        END MODULE GET_BDF_INSTANTWATER__genmod
