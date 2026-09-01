        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:27 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SLAM_ROMBERG__genmod
          INTERFACE 
            SUBROUTINE SLAM_ROMBERG(ROMK,A,B,CT,NUM0,POINT1,INTAB)
              INTEGER(KIND=4) :: NUM0
              INTEGER(KIND=4) :: ROMK
              REAL(KIND=8) :: A
              REAL(KIND=8) :: B
              REAL(KIND=8) :: CT
              REAL(KIND=8) :: POINT1(NUM0,2)
              REAL(KIND=8) :: INTAB
            END SUBROUTINE SLAM_ROMBERG
          END INTERFACE 
        END MODULE SLAM_ROMBERG__genmod
