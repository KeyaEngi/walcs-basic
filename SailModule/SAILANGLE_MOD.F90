!===============================================================================
! SAILANGLE_MOD
!
! Compute the representative sail chord direction, relative-wind upstream
! direction, signed angle of attack, and its 180-degree database mapping.
!===============================================================================
MODULE SAILANGLE_MOD
  USE SAILPARAM_MOD, ONLY: &
    DP, IDX_X, IDX_Y, IDX_Z, DEG_TO_RAD, RAD_TO_DEG, ANGLE_TOL_DEG, &
    SAIL_ANGLE_PERIOD_DEG, VREL_MIN, SAIL_OK, SAIL_ERR_INVALID_INPUT, &
    SAIL_ERR_LOW_WIND_SPEED
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeSailAngle

CONTAINS

  SUBROUTINE ComputeSailAngle( &
      V_REL_H_BODY, V_REL_H_MAG, DELTA_S_DEG, DELTA_NORMALIZED_DEG, &
      C_CHORD_BODY, E_UPSTREAM_BODY, ALPHA_RAW_DEG, ALPHA_DB_DEG, &
      IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: V_REL_H_BODY(3)
    REAL(DP), INTENT(IN) :: V_REL_H_MAG
    REAL(DP), INTENT(IN) :: DELTA_S_DEG
    REAL(DP), INTENT(OUT) :: DELTA_NORMALIZED_DEG
    REAL(DP), INTENT(OUT) :: C_CHORD_BODY(3)
    REAL(DP), INTENT(OUT) :: E_UPSTREAM_BODY(3)
    REAL(DP), INTENT(OUT) :: ALPHA_RAW_DEG
    REAL(DP), INTENT(OUT) :: ALPHA_DB_DEG
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    REAL(DP) :: V_MAG_CHECK
    REAL(DP) :: MAG_TOL
    REAL(DP) :: DELTA_RAD
    REAL(DP) :: CHORD_NORM
    REAL(DP) :: UPSTREAM_NORM
    REAL(DP) :: CROSS_Z
    REAL(DP) :: DOT_VALUE
    REAL(DP) :: CHORD_CANDIDATE(3)
    REAL(DP) :: UPSTREAM_CANDIDATE(3)
    REAL(DP) :: DELTA_CANDIDATE
    REAL(DP) :: ALPHA_RAW_CANDIDATE
    REAL(DP) :: ALPHA_DB_CANDIDATE
    CHARACTER(LEN=512) :: DIAGNOSTIC

    DELTA_NORMALIZED_DEG = 0.0_DP
    C_CHORD_BODY = 0.0_DP
    E_UPSTREAM_BODY = 0.0_DP
    ALPHA_RAW_DEG = 0.0_DP
    ALPHA_DB_DEG = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''

    IF (.NOT. ALL(IEEE_IS_FINITE(V_REL_H_BODY))) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind vector contains NaN or infinity.'
      RETURN
    END IF
    IF (.NOT. IEEE_IS_FINITE(V_REL_H_MAG)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind speed is NaN or infinity.'
      RETURN
    END IF
    IF (.NOT. IEEE_IS_FINITE(DELTA_S_DEG)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: sail angle is NaN or infinity.'
      RETURN
    END IF
    IF (V_REL_H_MAG < 0.0_DP) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind speed must be nonnegative.'
      RETURN
    END IF

    V_MAG_CHECK = SQRT(V_REL_H_BODY(IDX_X)**2 + &
                       V_REL_H_BODY(IDX_Y)**2)
    IF (.NOT. IEEE_IS_FINITE(V_MAG_CHECK)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind vector magnitude overflowed.'
      RETURN
    END IF
    MAG_TOL = MAX(1.0E-12_DP, &
                  1.0E-10_DP * MAX(1.0_DP, V_MAG_CHECK))

    IF (ABS(V_REL_H_BODY(IDX_Z)) > MAG_TOL) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: horizontal relative wind has a nonzero z component.'
      RETURN
    END IF
    IF (ABS(V_REL_H_MAG - V_MAG_CHECK) > MAG_TOL) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind scalar and vector magnitudes disagree.'
      RETURN
    END IF
    IF (V_REL_H_MAG < VREL_MIN) THEN
      IERR = SAIL_ERR_LOW_WIND_SPEED
      MESSAGE = 'Horizontal relative wind is too small; angle of attack is undefined.'
      RETURN
    END IF

    DELTA_CANDIDATE = NormalizePeriodicAngle(DELTA_S_DEG)
    DELTA_RAD = DELTA_CANDIDATE * DEG_TO_RAD

    CHORD_CANDIDATE = 0.0_DP
    CHORD_CANDIDATE(IDX_X) = -SIN(DELTA_RAD)
    CHORD_CANDIDATE(IDX_Y) = -COS(DELTA_RAD)

    UPSTREAM_CANDIDATE = 0.0_DP
    UPSTREAM_CANDIDATE(IDX_X) = &
      -V_REL_H_BODY(IDX_X) / V_REL_H_MAG
    UPSTREAM_CANDIDATE(IDX_Y) = &
      -V_REL_H_BODY(IDX_Y) / V_REL_H_MAG

    CHORD_NORM = VectorNorm2D(CHORD_CANDIDATE)
    UPSTREAM_NORM = VectorNorm2D(UPSTREAM_CANDIDATE)
    IF (ABS(CHORD_NORM - 1.0_DP) > 1.0E-10_DP .OR. &
        ABS(UPSTREAM_NORM - 1.0_DP) > 1.0E-10_DP) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid result: chord or upstream direction is not a unit vector.'
      RETURN
    END IF

    CROSS_Z = CHORD_CANDIDATE(IDX_X) * UPSTREAM_CANDIDATE(IDX_Y) - &
              CHORD_CANDIDATE(IDX_Y) * UPSTREAM_CANDIDATE(IDX_X)
    DOT_VALUE = CHORD_CANDIDATE(IDX_X) * UPSTREAM_CANDIDATE(IDX_X) + &
                CHORD_CANDIDATE(IDX_Y) * UPSTREAM_CANDIDATE(IDX_Y)
    ALPHA_RAW_CANDIDATE = ATAN2(CROSS_Z, DOT_VALUE) * RAD_TO_DEG
    ALPHA_DB_CANDIDATE = NormalizePeriodicAngle(ALPHA_RAW_CANDIDATE)

    IF (.NOT. IEEE_IS_FINITE(DELTA_CANDIDATE) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(CHORD_CANDIDATE)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(UPSTREAM_CANDIDATE)) .OR. &
        .NOT. IEEE_IS_FINITE(ALPHA_RAW_CANDIDATE) .OR. &
        .NOT. IEEE_IS_FINITE(ALPHA_DB_CANDIDATE)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid result: sail-angle computation produced NaN or infinity.'
      RETURN
    END IF

    DELTA_NORMALIZED_DEG = DELTA_CANDIDATE
    C_CHORD_BODY = CHORD_CANDIDATE
    E_UPSTREAM_BODY = UPSTREAM_CANDIDATE
    ALPHA_RAW_DEG = ALPHA_RAW_CANDIDATE
    ALPHA_DB_DEG = ALPHA_DB_CANDIDATE

    WRITE(DIAGNOSTIC, '(A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,A)') &
      'Sail angle computed successfully; delta=', DELTA_S_DEG, &
      ' deg, normalized delta=', DELTA_NORMALIZED_DEG, &
      ' deg, alpha_raw=', ALPHA_RAW_DEG, ', alpha_db=', ALPHA_DB_DEG, &
      ' deg, wind speed=', V_REL_H_MAG, ' m/s.'
    MESSAGE = TRIM(DIAGNOSTIC)
  END SUBROUTINE ComputeSailAngle


  PURE FUNCTION NormalizePeriodicAngle(ANGLE_DEG) RESULT(NORMALIZED_DEG)
    REAL(DP), INTENT(IN) :: ANGLE_DEG
    REAL(DP) :: NORMALIZED_DEG

    NORMALIZED_DEG = MODULO(ANGLE_DEG, SAIL_ANGLE_PERIOD_DEG)
    IF (NORMALIZED_DEG <= ANGLE_TOL_DEG .OR. &
        SAIL_ANGLE_PERIOD_DEG - NORMALIZED_DEG <= ANGLE_TOL_DEG) THEN
      NORMALIZED_DEG = 0.0_DP
    END IF
  END FUNCTION NormalizePeriodicAngle


  PURE FUNCTION VectorNorm2D(VECTOR) RESULT(NORM_VALUE)
    REAL(DP), INTENT(IN) :: VECTOR(3)
    REAL(DP) :: NORM_VALUE

    NORM_VALUE = SQRT(VECTOR(IDX_X)**2 + VECTOR(IDX_Y)**2)
  END FUNCTION VectorNorm2D

END MODULE SAILANGLE_MOD
