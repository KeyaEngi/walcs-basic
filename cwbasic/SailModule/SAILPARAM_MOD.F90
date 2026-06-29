MODULE SAILPARAM_MOD
!------------------------------------------------------------------------------
! Shared parameters and small utilities for the independent sail module.
! This module is intentionally not connected to the WALCS main program yet.
!------------------------------------------------------------------------------
IMPLICIT NONE
SAVE

INTEGER, PARAMETER :: SAIL_DP = SELECTED_REAL_KIND(12, 60)
INTEGER, PARAMETER :: SAIL_MAX_DB_POINTS = 512
REAL(SAIL_DP), PARAMETER :: SAIL_PI = 3.1415926535897932384626433832795_SAIL_DP

REAL(SAIL_DP) :: SailAirDensity = 1.225_SAIL_DP

PUBLIC :: SAIL_DP
PUBLIC :: SAIL_MAX_DB_POINTS
PUBLIC :: SAIL_PI
PUBLIC :: SailAirDensity
PUBLIC :: SailDegToRad
PUBLIC :: SailRadToDeg

CONTAINS

FUNCTION SailDegToRad(angleDeg) RESULT(angleRad)
IMPLICIT NONE
REAL(SAIL_DP), INTENT(IN) :: angleDeg
REAL(SAIL_DP) :: angleRad

angleRad = angleDeg * SAIL_PI / 180.0_SAIL_DP

END FUNCTION SailDegToRad

FUNCTION SailRadToDeg(angleRad) RESULT(angleDeg)
IMPLICIT NONE
REAL(SAIL_DP), INTENT(IN) :: angleRad
REAL(SAIL_DP) :: angleDeg

angleDeg = angleRad * 180.0_SAIL_DP / SAIL_PI

END FUNCTION SailRadToDeg

END MODULE SAILPARAM_MOD
