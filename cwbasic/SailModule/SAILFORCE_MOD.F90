MODULE SAILFORCE_MOD
!------------------------------------------------------------------------------
! Convert sail aerodynamic coefficients to a six-degree-of-freedom external load.
!
! Coordinate convention:
!   X forward, Y to port, Z upward.
!   apparentWindAngleDeg is measured in the XY plane from +X toward +Y.
!   Drag acts along the apparent wind velocity direction.
!   Lift is taken 90 degrees counter-clockwise from the drag direction in XY.
!
! Output:
!   Fsail(1,1:6) = real part of force/moment vector.
!   Fsail(2,1:6) = imaginary part, currently zero for quasi-steady loading.
!------------------------------------------------------------------------------
USE SAILPARAM_MOD, ONLY: SAIL_DP, SailAirDensity, SailDegToRad

IMPLICIT NONE
SAVE
PRIVATE

PUBLIC :: GetSailForce

CONTAINS

SUBROUTINE GetSailForce(clValue, cdValue, apparentWindSpeed, apparentWindAngleDeg, &
                        sailArea, xCE, yCE, zCE, Fsail, ierr)
IMPLICIT NONE
REAL(SAIL_DP), INTENT(IN) :: clValue
REAL(SAIL_DP), INTENT(IN) :: cdValue
REAL(SAIL_DP), INTENT(IN) :: apparentWindSpeed
REAL(SAIL_DP), INTENT(IN) :: apparentWindAngleDeg
REAL(SAIL_DP), INTENT(IN) :: sailArea
REAL(SAIL_DP), INTENT(IN) :: xCE
REAL(SAIL_DP), INTENT(IN) :: yCE
REAL(SAIL_DP), INTENT(IN) :: zCE
REAL(SAIL_DP), DIMENSION(1:2,1:6), INTENT(OUT) :: Fsail
INTEGER, INTENT(OUT) :: ierr

REAL(SAIL_DP) :: qDyn
REAL(SAIL_DP) :: angleRad
REAL(SAIL_DP) :: dragForce
REAL(SAIL_DP) :: liftForce
REAL(SAIL_DP) :: fx
REAL(SAIL_DP) :: fy
REAL(SAIL_DP) :: fz

ierr = 0
Fsail = 0.0_SAIL_DP

IF (apparentWindSpeed < 0.0_SAIL_DP) THEN
    ierr = 1
    RETURN
END IF

IF (sailArea < 0.0_SAIL_DP) THEN
    ierr = 2
    RETURN
END IF

angleRad = SailDegToRad(apparentWindAngleDeg)
qDyn = 0.5_SAIL_DP * SailAirDensity * apparentWindSpeed * apparentWindSpeed * sailArea
dragForce = qDyn * cdValue
liftForce = qDyn * clValue

fx = dragForce * COS(angleRad) - liftForce * SIN(angleRad)
fy = dragForce * SIN(angleRad) + liftForce * COS(angleRad)
fz = 0.0_SAIL_DP

Fsail(1,1) = fx
Fsail(1,2) = fy
Fsail(1,3) = fz
Fsail(1,4) = yCE * fz - zCE * fy
Fsail(1,5) = zCE * fx - xCE * fz
Fsail(1,6) = xCE * fy - yCE * fx

END SUBROUTINE GetSailForce

END MODULE SAILFORCE_MOD
