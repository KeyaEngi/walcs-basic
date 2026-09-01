        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:00 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GET_DYDX__genmod
          INTERFACE 
            SUBROUTINE GET_DYDX(IT,CA,AIR_NUM,AIR_POINT,DYDX)
              USE LIFT, ONLY :                                          &
     &          DYDX
              INTEGER(KIND=4) :: AIR_NUM
              INTEGER(KIND=4) :: IT
              REAL(KIND=8) :: CA
              REAL(KIND=8) :: AIR_POINT(AIR_NUM,2)
              REAL(KIND=8) :: DYDX(IT)
            END SUBROUTINE GET_DYDX
          END INTERFACE 
        END MODULE GET_DYDX__genmod
