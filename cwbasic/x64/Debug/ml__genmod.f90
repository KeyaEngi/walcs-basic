        !COMPILER-GENERATED INTERFACE MODULE: Wed May 02 20:59:27 2018
        MODULE ML__genmod
          INTERFACE 
            SUBROUTINE ML(N1,N2,T,H,XG,ZG,XL,ZL,REF,OMG,EL1,W1,CML)
              INTEGER(KIND=4), INTENT(IN) :: N2
              INTEGER(KIND=4), INTENT(IN) :: N1
              REAL(KIND=8), INTENT(IN) :: T
              REAL(KIND=8), INTENT(IN) :: H
              REAL(KIND=8), INTENT(IN) :: XG
              REAL(KIND=8), INTENT(IN) :: ZG
              REAL(KIND=8), INTENT(IN) :: XL
              REAL(KIND=8), INTENT(IN) :: ZL
              REAL(KIND=8), INTENT(INOUT) :: REF(1:N1)
              REAL(KIND=8), INTENT(IN) :: OMG(1:N2)
              REAL(KIND=8), INTENT(IN) :: EL1(1:N2)
              REAL(KIND=8), INTENT(INOUT) :: W1(1:N2)
              REAL(KIND=8), INTENT(OUT) :: CML(1:,1:)
            END SUBROUTINE ML
          END INTERFACE 
        END MODULE ML__genmod
