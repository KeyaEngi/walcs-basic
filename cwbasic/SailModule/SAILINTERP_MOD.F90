MODULE SAILINTERP_MOD
!------------------------------------------------------------------------------
! Linear interpolation for sail aerodynamic coefficients.
!------------------------------------------------------------------------------
USE SAILPARAM_MOD, ONLY: SAIL_DP
USE SAILDATABASE_MOD, ONLY: SailDbN, SailDbAlphaDeg, SailDbCL, SailDbCD

IMPLICIT NONE
SAVE
PRIVATE

PUBLIC :: GetSailCoeff

CONTAINS

SUBROUTINE GetSailCoeff(alphaDeg, clValue, cdValue, ierr)
IMPLICIT NONE
REAL(SAIL_DP), INTENT(IN) :: alphaDeg
REAL(SAIL_DP), INTENT(OUT) :: clValue
REAL(SAIL_DP), INTENT(OUT) :: cdValue
INTEGER, INTENT(OUT) :: ierr

INTEGER :: i
REAL(SAIL_DP) :: ratio

ierr = 0
clValue = 0.0_SAIL_DP
cdValue = 0.0_SAIL_DP

IF (SailDbN <= 0) THEN
    ierr = 1
    RETURN
END IF

IF (SailDbN == 1) THEN
    clValue = SailDbCL(1)
    cdValue = SailDbCD(1)
    RETURN
END IF

IF (alphaDeg <= SailDbAlphaDeg(1)) THEN
    clValue = SailDbCL(1)
    cdValue = SailDbCD(1)
    RETURN
END IF

IF (alphaDeg >= SailDbAlphaDeg(SailDbN)) THEN
    clValue = SailDbCL(SailDbN)
    cdValue = SailDbCD(SailDbN)
    RETURN
END IF

DO i = 1, SailDbN - 1
    IF (alphaDeg >= SailDbAlphaDeg(i) .AND. alphaDeg <= SailDbAlphaDeg(i + 1)) THEN
        ratio = (alphaDeg - SailDbAlphaDeg(i)) / &
                (SailDbAlphaDeg(i + 1) - SailDbAlphaDeg(i))
        clValue = SailDbCL(i) + ratio * (SailDbCL(i + 1) - SailDbCL(i))
        cdValue = SailDbCD(i) + ratio * (SailDbCD(i + 1) - SailDbCD(i))
        RETURN
    END IF
END DO

ierr = 2

END SUBROUTINE GetSailCoeff

END MODULE SAILINTERP_MOD
