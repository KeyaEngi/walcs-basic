        !COMPILER-GENERATED INTERFACE MODULE: Wed May 02 20:59:02 2018
        MODULE STIFFNESS__genmod
          INTERFACE 
            SUBROUTINE STIFFNESS(NT,TFH,TS,ASH,ASL,ASM0,ASN0,ASN,RXX,RXY&
     &,RYX,RYY)
              INTEGER(KIND=4) :: NT
              REAL(KIND=8) :: TFH
              REAL(KIND=8) :: TS
              REAL(KIND=8) :: ASH(10)
              REAL(KIND=8) :: ASL(10)
              REAL(KIND=8) :: ASM0(10)
              REAL(KIND=8) :: ASN0(10)
              REAL(KIND=8) :: ASN(10)
              REAL(KIND=8) :: RXX
              REAL(KIND=8) :: RXY
              REAL(KIND=8) :: RYX
              REAL(KIND=8) :: RYY
            END SUBROUTINE STIFFNESS
          END INTERFACE 
        END MODULE STIFFNESS__genmod
