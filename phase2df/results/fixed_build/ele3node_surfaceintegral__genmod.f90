        !COMPILER-GENERATED INTERFACE MODULE: Wed Aug 12 14:49:02 2026
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE ELE3NODE_SURFACEINTEGRAL__genmod
          INTERFACE 
            SUBROUTINE ELE3NODE_SURFACEINTEGRAL(NR,NODE,PS,UR,TEMFORCE)
              INTEGER(KIND=4) :: NR
              REAL(KIND=8) :: NODE(3,3)
              REAL(KIND=8) :: PS(3)
              REAL(KIND=8) :: UR(3,NR,3)
              REAL(KIND=8) :: TEMFORCE(NR)
            END SUBROUTINE ELE3NODE_SURFACEINTEGRAL
          END INTERFACE 
        END MODULE ELE3NODE_SURFACEINTEGRAL__genmod
