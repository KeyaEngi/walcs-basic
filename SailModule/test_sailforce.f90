PROGRAM TEST_SAILFORCE
  USE SAILPARAM_MOD
  USE SAILFORCE_MOD
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE

  REAL(DP), PARAMETER :: FORCE_TOL = 1.0E-10_DP
  REAL(DP), PARAMETER :: VECTOR_TOL = 1.0E-12_DP
  INTEGER :: N_PASS
  INTEGER :: N_FAIL

  N_PASS = 0
  N_FAIL = 0

  CALL TestStarBaseline
  CALL TestPureDrag
  CALL TestPositiveLift
  CALL TestNegativeLift
  CALL TestPositiveXWind
  CALL TestNegativeXWind
  CALL TestNegativeYWind
  CALL TestObliqueWind
  CALL TestCenterOfGravity
  CALL TestZeroCoefficients
  CALL TestSpeedSquaredScaling
  CALL TestLowWind
  CALL TestMinimumWind
  CALL TestMagnitudeMismatch
  CALL TestNonzeroZ
  CALL TestNegativeDrag
  CALL TestTinyNegativeDrag
  CALL TestNonFiniteInputs
  CALL TestCrossProductOrder
  CALL TestLoadConsistency
  CALL TestDirectionScan
  CALL TestRepeatedCalls
  CALL TestInputsUnchanged

  WRITE(*, '(A)') '========================================'
  WRITE(*, '(A)') 'SAILFORCE TEST SUMMARY'
  WRITE(*, '(A,I0)') 'PASS: ', N_PASS
  WRITE(*, '(A,I0)') 'FAIL: ', N_FAIL
  WRITE(*, '(A)') '========================================'

  IF (N_FAIL > 0) STOP 1

CONTAINS

  SUBROUTINE RecordResult(NAME, PASSED)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    LOGICAL, INTENT(IN) :: PASSED

    IF (PASSED) THEN
      N_PASS = N_PASS + 1
      WRITE(*, '(A,A)') '[PASS] ', TRIM(NAME)
    ELSE
      N_FAIL = N_FAIL + 1
      WRITE(*, '(A,A)') '[FAIL] ', TRIM(NAME)
    END IF
  END SUBROUTINE RecordResult


  LOGICAL FUNCTION NearScalar(A, B, TOLERANCE)
    REAL(DP), INTENT(IN) :: A
    REAL(DP), INTENT(IN) :: B
    REAL(DP), INTENT(IN) :: TOLERANCE

    NearScalar = ABS(A - B) <= TOLERANCE
  END FUNCTION NearScalar


  LOGICAL FUNCTION NearVector(A, B, TOLERANCE)
    REAL(DP), INTENT(IN) :: A(:)
    REAL(DP), INTENT(IN) :: B(:)
    REAL(DP), INTENT(IN) :: TOLERANCE

    NearVector = SIZE(A) == SIZE(B)
    IF (NearVector) NearVector = MAXVAL(ABS(A - B)) <= TOLERANCE
  END FUNCTION NearVector


  LOGICAL FUNCTION OutputsAreZero(E_DRAG, E_LIFT, Q, FORCE, MOMENT, LOAD)
    REAL(DP), INTENT(IN) :: E_DRAG(3), E_LIFT(3), Q
    REAL(DP), INTENT(IN) :: FORCE(3), MOMENT(3), LOAD(N_DOF)

    OutputsAreZero = MAXVAL(ABS(E_DRAG)) <= FORCE_TOL .AND. &
      MAXVAL(ABS(E_LIFT)) <= FORCE_TOL .AND. ABS(Q) <= FORCE_TOL .AND. &
      MAXVAL(ABS(FORCE)) <= FORCE_TOL .AND. &
      MAXVAL(ABS(MOMENT)) <= FORCE_TOL .AND. &
      MAXVAL(ABS(LOAD)) <= FORCE_TOL
  END FUNCTION OutputsAreZero


  SUBROUTINE CrossProductReference(R, F, M)
    REAL(DP), INTENT(IN) :: R(3), F(3)
    REAL(DP), INTENT(OUT) :: M(3)

    M(1) = R(2) * F(3) - R(3) * F(2)
    M(2) = R(3) * F(1) - R(1) * F(3)
    M(3) = R(1) * F(2) - R(2) * F(1)
  END SUBROUTINE CrossProductReference


  SUBROUTINE ReportFailureDetails(ACTUAL_IERR, EXPECTED_IERR, &
      E_DRAG, EXPECTED_DRAG, E_LIFT, EXPECTED_LIFT, Q, EXPECTED_Q, &
      FORCE, EXPECTED_FORCE, MOMENT, EXPECTED_MOMENT, LOAD, &
      EXPECTED_LOAD, MESSAGE)
    INTEGER, INTENT(IN) :: ACTUAL_IERR, EXPECTED_IERR
    REAL(DP), INTENT(IN), CONTIGUOUS :: E_DRAG(:), EXPECTED_DRAG(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: E_LIFT(:), EXPECTED_LIFT(:)
    REAL(DP), INTENT(IN) :: Q, EXPECTED_Q
    REAL(DP), INTENT(IN), CONTIGUOUS :: FORCE(:), EXPECTED_FORCE(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: MOMENT(:), EXPECTED_MOMENT(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: LOAD(:), EXPECTED_LOAD(:)
    CHARACTER(LEN=*), INTENT(IN) :: MESSAGE

    WRITE(*, '(A,I0,A,I0)') '  actual/expected IERR: ', ACTUAL_IERR, '/', &
      EXPECTED_IERR
    WRITE(*, '(A,3ES16.7)') '  actual E_DRAG:   ', E_DRAG
    WRITE(*, '(A,3ES16.7)') '  expected E_DRAG: ', EXPECTED_DRAG
    WRITE(*, '(A,3ES16.7)') '  actual E_LIFT:   ', E_LIFT
    WRITE(*, '(A,3ES16.7)') '  expected E_LIFT: ', EXPECTED_LIFT
    WRITE(*, '(A,ES16.7,A,ES16.7)') '  actual/expected q: ', Q, '/', EXPECTED_Q
    WRITE(*, '(A,3ES16.7)') '  actual FORCE:   ', FORCE
    WRITE(*, '(A,3ES16.7)') '  expected FORCE: ', EXPECTED_FORCE
    WRITE(*, '(A,3ES16.7)') '  actual MOMENT:   ', MOMENT
    WRITE(*, '(A,3ES16.7)') '  expected MOMENT: ', EXPECTED_MOMENT
    WRITE(*, '(A,6ES16.7)') '  actual LOAD:   ', LOAD
    WRITE(*, '(A,6ES16.7)') '  expected LOAD: ', EXPECTED_LOAD
    WRITE(*, '(A,A)') '  MESSAGE: ', TRIM(MESSAGE)
  END SUBROUTINE ReportFailureDetails


  SUBROUTINE RunForceCase(NAME, WIND, SPEED, CL_VALUE, CD_VALUE, R, &
      EXPECTED_IERR, EXPECTED_DRAG, EXPECTED_LIFT, EXPECTED_Q, &
      EXPECTED_FORCE, EXPECTED_MOMENT)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    REAL(DP), INTENT(IN), CONTIGUOUS :: WIND(:), R(:)
    REAL(DP), INTENT(IN) :: SPEED, CL_VALUE, CD_VALUE
    INTEGER, INTENT(IN) :: EXPECTED_IERR
    REAL(DP), INTENT(IN), CONTIGUOUS :: EXPECTED_DRAG(:), EXPECTED_LIFT(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: EXPECTED_FORCE(:), EXPECTED_MOMENT(:)
    REAL(DP), INTENT(IN) :: EXPECTED_Q
    REAL(DP) :: E_DRAG(3), E_LIFT(3), Q, FORCE(3), MOMENT(3), LOAD(N_DOF)
    REAL(DP) :: EXPECTED_LOAD(N_DOF)
    INTEGER :: IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    EXPECTED_LOAD(1:3) = EXPECTED_FORCE
    EXPECTED_LOAD(4:6) = EXPECTED_MOMENT
    CALL ComputeSailForce(WIND, SPEED, CL_VALUE, CD_VALUE, R, E_DRAG, &
      E_LIFT, Q, FORCE, MOMENT, LOAD, IERR, MESSAGE)
    OK = IERR == EXPECTED_IERR .AND. &
      NearVector(E_DRAG, EXPECTED_DRAG, VECTOR_TOL) .AND. &
      NearVector(E_LIFT, EXPECTED_LIFT, VECTOR_TOL) .AND. &
      NearScalar(Q, EXPECTED_Q, FORCE_TOL) .AND. &
      NearVector(FORCE, EXPECTED_FORCE, FORCE_TOL) .AND. &
      NearVector(MOMENT, EXPECTED_MOMENT, FORCE_TOL) .AND. &
      NearVector(LOAD, EXPECTED_LOAD, FORCE_TOL)
    IF (.NOT. OK) THEN
      CALL ReportFailureDetails(IERR, EXPECTED_IERR, E_DRAG, EXPECTED_DRAG, &
        E_LIFT, EXPECTED_LIFT, Q, EXPECTED_Q, FORCE, EXPECTED_FORCE, &
        MOMENT, EXPECTED_MOMENT, LOAD, EXPECTED_LOAD, MESSAGE)
    END IF
    CALL RecordResult(NAME, OK)
  END SUBROUTINE RunForceCase


  SUBROUTINE RunInvalidCase(NAME, WIND, SPEED, CL_VALUE, CD_VALUE, R)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    REAL(DP), INTENT(IN), CONTIGUOUS :: WIND(:), R(:)
    REAL(DP), INTENT(IN) :: SPEED, CL_VALUE, CD_VALUE
    REAL(DP) :: E_DRAG(3), E_LIFT(3), Q, FORCE(3), MOMENT(3), LOAD(N_DOF)
    INTEGER :: IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK
    REAL(DP) :: ZERO3(3), ZERO6(N_DOF)

    ZERO3 = 0.0_DP
    ZERO6 = 0.0_DP
    CALL ComputeSailForce(WIND, SPEED, CL_VALUE, CD_VALUE, R, E_DRAG, &
      E_LIFT, Q, FORCE, MOMENT, LOAD, IERR, MESSAGE)
    OK = IERR == SAIL_ERR_INVALID_INPUT .AND. &
      OutputsAreZero(E_DRAG, E_LIFT, Q, FORCE, MOMENT, LOAD)
    IF (.NOT. OK) THEN
      CALL ReportFailureDetails(IERR, SAIL_ERR_INVALID_INPUT, E_DRAG, &
        ZERO3, E_LIFT, ZERO3, Q, 0.0_DP, FORCE, ZERO3, MOMENT, ZERO3, &
        LOAD, ZERO6, MESSAGE)
    END IF
    CALL RecordResult(NAME, OK)
  END SUBROUTINE RunInvalidCase


  SUBROUTINE TestStarBaseline
    REAL(DP) :: WIND(3), R(3), E_D(3), E_L(3), F(3), M(3)

    WIND = (/ 0.0_DP, 10.0_DP, 0.0_DP /)
    R = (/ 2.0_DP, 3.0_DP, 4.0_DP /)
    E_D = (/ 0.0_DP, 1.0_DP, 0.0_DP /)
    E_L = (/ 1.0_DP, 0.0_DP, 0.0_DP /)
    F = (/ 183.75_DP, 91.875_DP, 0.0_DP /)
    M = (/ -367.5_DP, 735.0_DP, -367.5_DP /)
    CALL RunForceCase('STAR baseline direction and signs', WIND, 10.0_DP, &
      0.2_DP, 0.1_DP, R, SAIL_OK, E_D, E_L, 61.25_DP, F, M)
    WRITE(*, '(A,3F12.3)') '  STAR force  = ', F
    WRITE(*, '(A,3F12.3)') '  STAR moment = ', M
  END SUBROUTINE TestStarBaseline


  SUBROUTINE TestPureDrag
    CALL RunForceCase('pure drag', (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
      10.0_DP, 0.0_DP, 0.1_DP, (/ 2.0_DP, 3.0_DP, 4.0_DP /), SAIL_OK, &
      (/ 0.0_DP, 1.0_DP, 0.0_DP /), (/ 1.0_DP, 0.0_DP, 0.0_DP /), &
      61.25_DP, (/ 0.0_DP, 91.875_DP, 0.0_DP /), &
      (/ -367.5_DP, 0.0_DP, 183.75_DP /))
  END SUBROUTINE TestPureDrag


  SUBROUTINE TestPositiveLift
    CALL RunForceCase('positive signed lift', (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
      10.0_DP, 0.2_DP, 0.0_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 0.0_DP, 1.0_DP, 0.0_DP /), (/ 1.0_DP, 0.0_DP, 0.0_DP /), &
      61.25_DP, (/ 183.75_DP, 0.0_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestPositiveLift


  SUBROUTINE TestNegativeLift
    CALL RunForceCase('negative signed lift', (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
      10.0_DP, -0.2_DP, 0.0_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 0.0_DP, 1.0_DP, 0.0_DP /), (/ 1.0_DP, 0.0_DP, 0.0_DP /), &
      61.25_DP, (/ -183.75_DP, 0.0_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestNegativeLift


  SUBROUTINE TestPositiveXWind
    CALL RunForceCase('relative wind along +x_b', (/ 10.0_DP, 0.0_DP, 0.0_DP /), &
      10.0_DP, 0.2_DP, 0.1_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 1.0_DP, 0.0_DP, 0.0_DP /), (/ 0.0_DP, -1.0_DP, 0.0_DP /), &
      61.25_DP, (/ 91.875_DP, -183.75_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestPositiveXWind


  SUBROUTINE TestNegativeXWind
    CALL RunForceCase('relative wind along -x_b', (/ -10.0_DP, 0.0_DP, 0.0_DP /), &
      10.0_DP, 0.2_DP, 0.1_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ -1.0_DP, 0.0_DP, 0.0_DP /), (/ 0.0_DP, 1.0_DP, 0.0_DP /), &
      61.25_DP, (/ -91.875_DP, 183.75_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestNegativeXWind


  SUBROUTINE TestNegativeYWind
    CALL RunForceCase('relative wind along -y_b', (/ 0.0_DP, -10.0_DP, 0.0_DP /), &
      10.0_DP, 0.2_DP, 0.1_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 0.0_DP, -1.0_DP, 0.0_DP /), (/ -1.0_DP, 0.0_DP, 0.0_DP /), &
      61.25_DP, (/ -183.75_DP, -91.875_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestNegativeYWind


  SUBROUTINE TestObliqueWind
    CALL RunForceCase('oblique relative wind', (/ 6.0_DP, 8.0_DP, 0.0_DP /), &
      10.0_DP, 0.3_DP, 0.2_DP, (/ 1.0_DP, 2.0_DP, 3.0_DP /), SAIL_OK, &
      (/ 0.6_DP, 0.8_DP, 0.0_DP /), (/ 0.8_DP, -0.6_DP, 0.0_DP /), &
      61.25_DP, (/ 330.75_DP, -18.375_DP, 0.0_DP /), &
      (/ 55.125_DP, 992.25_DP, -679.875_DP /))
  END SUBROUTINE TestObliqueWind


  SUBROUTINE TestCenterOfGravity
    CALL RunForceCase('sail force applied at CG', (/ 6.0_DP, 8.0_DP, 0.0_DP /), &
      10.0_DP, 0.3_DP, 0.2_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 0.6_DP, 0.8_DP, 0.0_DP /), (/ 0.8_DP, -0.6_DP, 0.0_DP /), &
      61.25_DP, (/ 330.75_DP, -18.375_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestCenterOfGravity


  SUBROUTINE TestZeroCoefficients
    CALL RunForceCase('zero aerodynamic coefficients', &
      (/ 6.0_DP, 8.0_DP, 0.0_DP /), 10.0_DP, 0.0_DP, 0.0_DP, &
      (/ 1.0_DP, 2.0_DP, 3.0_DP /), SAIL_OK, &
      (/ 0.6_DP, 0.8_DP, 0.0_DP /), (/ 0.8_DP, -0.6_DP, 0.0_DP /), &
      61.25_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestZeroCoefficients


  SUBROUTINE TestSpeedSquaredScaling
    REAL(DP) :: E_D5(3), E_L5(3), F5(3), M5(3), L5(N_DOF), Q5
    REAL(DP) :: E_D10(3), E_L10(3), F10(3), M10(3), L10(N_DOF), Q10
    REAL(DP) :: WIND5(3), WIND10(3), R(3)
    REAL(DP) :: EXPECTED_F(3), EXPECTED_M(3), EXPECTED_LOAD(N_DOF)
    INTEGER :: IERR5, IERR10
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    WIND5 = (/ 3.0_DP, 4.0_DP, 0.0_DP /)
    WIND10 = (/ 6.0_DP, 8.0_DP, 0.0_DP /)
    R = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    CALL ComputeSailForce(WIND5, 5.0_DP, &
      0.3_DP, 0.2_DP, R, E_D5, E_L5, &
      Q5, F5, M5, L5, IERR5, MESSAGE)
    CALL ComputeSailForce(WIND10, 10.0_DP, &
      0.3_DP, 0.2_DP, R, E_D10, E_L10, &
      Q10, F10, M10, L10, IERR10, MESSAGE)
    EXPECTED_F = 4.0_DP * F5
    EXPECTED_M = 4.0_DP * M5
    EXPECTED_LOAD = 4.0_DP * L5
    OK = IERR5 == SAIL_OK .AND. IERR10 == SAIL_OK .AND. &
      NearScalar(Q10, 4.0_DP * Q5, FORCE_TOL) .AND. &
      NearVector(F10, EXPECTED_F, FORCE_TOL) .AND. &
      NearVector(M10, EXPECTED_M, FORCE_TOL) .AND. &
      NearVector(L10, EXPECTED_LOAD, FORCE_TOL)
    WRITE(*, '(A,F8.3)') '  q(10)/q(5) expected = 4; actual = ', Q10 / Q5
    IF (.NOT. OK) CALL ReportFailureDetails(IERR10, SAIL_OK, E_D10, &
      E_D5, E_L10, E_L5, Q10, 4.0_DP * Q5, F10, EXPECTED_F, &
      M10, EXPECTED_M, L10, EXPECTED_LOAD, MESSAGE)
    CALL RecordResult('wind-speed squared scaling', OK)
  END SUBROUTINE TestSpeedSquaredScaling


  SUBROUTINE TestLowWind
    REAL(DP) :: E_D(3), E_L(3), Q, F(3), M(3), LOAD(N_DOF)
    REAL(DP) :: WIND(3), R(3)
    REAL(DP) :: ZERO3(3), ZERO6(N_DOF)
    INTEGER :: IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    WIND = (/ 0.5_DP * VREL_MIN, 0.0_DP, 0.0_DP /)
    R = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    ZERO3 = 0.0_DP
    ZERO6 = 0.0_DP
    CALL ComputeSailForce(WIND, &
      0.5_DP * VREL_MIN, 0.2_DP, 0.1_DP, R, &
      E_D, E_L, Q, F, M, LOAD, IERR, MESSAGE)
    OK = IERR == SAIL_ERR_LOW_WIND_SPEED .AND. &
      OutputsAreZero(E_D, E_L, Q, F, M, LOAD)
    IF (.NOT. OK) CALL ReportFailureDetails(IERR, &
      SAIL_ERR_LOW_WIND_SPEED, E_D, ZERO3, E_L, ZERO3, Q, 0.0_DP, &
      F, ZERO3, M, ZERO3, LOAD, ZERO6, MESSAGE)
    CALL RecordResult('below VREL_MIN gives zero loads', OK)
  END SUBROUTINE TestLowWind


  SUBROUTINE TestMinimumWind
    REAL(DP) :: SCALE

    SCALE = 0.5_DP * RHO_AIR * VREL_MIN**2 * SAIL_AREA
    CALL RunForceCase('wind speed equal to VREL_MIN', &
      (/ VREL_MIN, 0.0_DP, 0.0_DP /), VREL_MIN, 0.2_DP, 0.1_DP, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 1.0_DP, 0.0_DP, 0.0_DP /), (/ 0.0_DP, -1.0_DP, 0.0_DP /), &
      0.5_DP * RHO_AIR * VREL_MIN**2, &
      (/ 0.1_DP * SCALE, -0.2_DP * SCALE, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestMinimumWind


  SUBROUTINE TestMagnitudeMismatch
    CALL RunInvalidCase('scalar/vector wind magnitude mismatch', &
      (/ 3.0_DP, 4.0_DP, 0.0_DP /), 10.0_DP, 0.2_DP, 0.1_DP, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestMagnitudeMismatch


  SUBROUTINE TestNonzeroZ
    CALL RunInvalidCase('nonzero horizontal-wind z component', &
      (/ 3.0_DP, 4.0_DP, 1.0_DP /), 5.0_DP, 0.2_DP, 0.1_DP, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestNonzeroZ


  SUBROUTINE TestNegativeDrag
    CALL RunInvalidCase('negative CD beyond tolerance', &
      (/ 0.0_DP, 10.0_DP, 0.0_DP /), 10.0_DP, 0.2_DP, -0.1_DP, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestNegativeDrag


  SUBROUTINE TestTinyNegativeDrag
    REAL(DP) :: SMALL_CD, SCALE

    SMALL_CD = -0.5_DP * DATABASE_TOL
    SCALE = 61.25_DP * SAIL_AREA
    CALL RunForceCase('tiny negative CD preserved', &
      (/ 0.0_DP, 10.0_DP, 0.0_DP /), 10.0_DP, 0.0_DP, SMALL_CD, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /), SAIL_OK, &
      (/ 0.0_DP, 1.0_DP, 0.0_DP /), (/ 1.0_DP, 0.0_DP, 0.0_DP /), &
      61.25_DP, (/ 0.0_DP, SMALL_CD * SCALE, 0.0_DP /), &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestTinyNegativeDrag


  SUBROUTINE TestNonFiniteInputs
    REAL(DP) :: NAN_VALUE, INF_VALUE

    NAN_VALUE = IEEE_VALUE(0.0_DP, IEEE_QUIET_NAN)
    INF_VALUE = IEEE_VALUE(0.0_DP, IEEE_POSITIVE_INF)
    CALL RunInvalidCase('NaN CL rejection', (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
      10.0_DP, NAN_VALUE, 0.1_DP, (/ 0.0_DP, 0.0_DP, 0.0_DP /))
    CALL RunInvalidCase('infinite CD rejection', (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
      10.0_DP, 0.2_DP, INF_VALUE, (/ 0.0_DP, 0.0_DP, 0.0_DP /))
    CALL RunInvalidCase('NaN wind-vector rejection', &
      (/ NAN_VALUE, 10.0_DP, 0.0_DP /), 10.0_DP, 0.2_DP, 0.1_DP, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
    CALL RunInvalidCase('infinite position-vector rejection', &
      (/ 0.0_DP, 10.0_DP, 0.0_DP /), 10.0_DP, 0.2_DP, 0.1_DP, &
      (/ 0.0_DP, INF_VALUE, 0.0_DP /))
    CALL RunInvalidCase('NaN wind-speed rejection', &
      (/ 0.0_DP, 10.0_DP, 0.0_DP /), NAN_VALUE, 0.2_DP, 0.1_DP, &
      (/ 0.0_DP, 0.0_DP, 0.0_DP /))
  END SUBROUTINE TestNonFiniteInputs


  SUBROUTINE TestCrossProductOrder
    REAL(DP) :: E_D(3), E_L(3), Q, F(3), M(3), LOAD(N_DOF), M_REF(3)
    REAL(DP) :: R(3), WIND(3)
    REAL(DP) :: EXPECTED_LOAD(N_DOF)
    INTEGER :: IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    R = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    WIND = (/ 6.0_DP, 8.0_DP, 0.0_DP /)
    CALL ComputeSailForce(WIND, 10.0_DP, &
      0.3_DP, 0.2_DP, R, E_D, E_L, Q, F, M, LOAD, IERR, MESSAGE)
    CALL CrossProductReference(R, F, M_REF)
    EXPECTED_LOAD(1:3) = F
    EXPECTED_LOAD(4:6) = M_REF
    OK = IERR == SAIL_OK .AND. NearVector(M, M_REF, FORCE_TOL) .AND. &
      NearVector(M, (/ 55.125_DP, 992.25_DP, -679.875_DP /), FORCE_TOL)
    WRITE(*, '(A,3F12.3)') '  independent R cross F = ', M_REF
    IF (.NOT. OK) CALL ReportFailureDetails(IERR, SAIL_OK, E_D, E_D, &
      E_L, E_L, Q, Q, F, F, M, M_REF, LOAD, EXPECTED_LOAD, MESSAGE)
    CALL RecordResult('moment cross-product order R cross F', OK)
  END SUBROUTINE TestCrossProductOrder


  SUBROUTINE TestLoadConsistency
    REAL(DP) :: E_D(3), E_L(3), Q, F(3), M(3), LOAD(N_DOF)
    REAL(DP) :: WIND(3), R(3)
    REAL(DP) :: EXPECTED_LOAD(N_DOF)
    INTEGER :: IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    WIND = (/ 6.0_DP, 8.0_DP, 0.0_DP /)
    R = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    CALL ComputeSailForce(WIND, 10.0_DP, &
      0.3_DP, 0.2_DP, R, E_D, E_L, &
      Q, F, M, LOAD, IERR, MESSAGE)
    EXPECTED_LOAD(1:3) = F
    EXPECTED_LOAD(4:6) = M
    OK = IERR == SAIL_OK .AND. NearVector(LOAD(1:3), F, FORCE_TOL) .AND. &
      NearVector(LOAD(4:6), M, FORCE_TOL)
    IF (.NOT. OK) CALL ReportFailureDetails(IERR, SAIL_OK, E_D, E_D, &
      E_L, E_L, Q, Q, F, F, M, M, LOAD, EXPECTED_LOAD, MESSAGE)
    CALL RecordResult('six-DOF force and moment consistency', OK)
  END SUBROUTINE TestLoadConsistency


  SUBROUTINE TestDirectionScan
    REAL(DP) :: ANGLE, WIND(3), E_D(3), E_L(3), Q
    REAL(DP) :: F(3), M(3), LOAD(N_DOF), EXPECTED_LIFT(3), R_ZERO(3)
    REAL(DP) :: DRAG_NORM, LIFT_NORM, DIRECTION_DOT
    REAL(DP) :: EXPECTED_DRAG(3), EXPECTED_FORCE(3), EXPECTED_MOMENT(3)
    REAL(DP) :: EXPECTED_LOAD(N_DOF), EXPECTED_Q
    INTEGER :: INDEX, IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    OK = .TRUE.
    R_ZERO = 0.0_DP
    DO INDEX = 0, 11
      ANGLE = REAL(30 * INDEX, DP) * DEG_TO_RAD
      WIND = 10.0_DP * (/ COS(ANGLE), SIN(ANGLE), 0.0_DP /)
      CALL ComputeSailForce(WIND, 10.0_DP, 0.2_DP, 0.1_DP, &
        R_ZERO, E_D, E_L, Q, F, M, LOAD, &
        IERR, MESSAGE)
      EXPECTED_LIFT = (/ E_D(2), -E_D(1), 0.0_DP /)
      EXPECTED_DRAG = WIND / 10.0_DP
      EXPECTED_Q = 0.5_DP * RHO_AIR * 10.0_DP**2
      EXPECTED_FORCE = EXPECTED_Q * SAIL_AREA * &
        (0.1_DP * EXPECTED_DRAG + 0.2_DP * EXPECTED_LIFT)
      EXPECTED_MOMENT = 0.0_DP
      EXPECTED_LOAD(1:3) = EXPECTED_FORCE
      EXPECTED_LOAD(4:6) = EXPECTED_MOMENT
      DRAG_NORM = SQRT(SUM(E_D**2))
      LIFT_NORM = SQRT(SUM(E_L**2))
      DIRECTION_DOT = DOT_PRODUCT(E_D, E_L)
      OK = OK .AND. IERR == SAIL_OK .AND. &
        NearScalar(DRAG_NORM, 1.0_DP, VECTOR_TOL) .AND. &
        NearScalar(LIFT_NORM, 1.0_DP, VECTOR_TOL) .AND. &
        ABS(DIRECTION_DOT) <= VECTOR_TOL .AND. &
        NearVector(E_L, EXPECTED_LIFT, VECTOR_TOL) .AND. &
        ALL(IEEE_IS_FINITE(E_D)) .AND. ALL(IEEE_IS_FINITE(E_L))
      IF (.NOT. OK) THEN
        CALL ReportFailureDetails(IERR, SAIL_OK, E_D, EXPECTED_DRAG, &
          E_L, EXPECTED_LIFT, Q, EXPECTED_Q, F, EXPECTED_FORCE, M, &
          EXPECTED_MOMENT, LOAD, EXPECTED_LOAD, MESSAGE)
        EXIT
      END IF
    END DO
    CALL RecordResult('12-direction unit and orthogonality scan', OK)
  END SUBROUTINE TestDirectionScan


  SUBROUTINE TestRepeatedCalls
    REAL(DP) :: E_D(3), E_L(3), Q, F(3), M(3), LOAD(N_DOF)
    REAL(DP) :: REF_E_D(3), REF_E_L(3), REF_Q, REF_F(3), REF_M(3)
    REAL(DP) :: REF_LOAD(N_DOF)
    REAL(DP) :: WIND(3), R(3)
    INTEGER :: ITERATION, IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    WIND = (/ 6.0_DP, 8.0_DP, 0.0_DP /)
    R = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    CALL ComputeSailForce(WIND, 10.0_DP, &
      0.3_DP, 0.2_DP, R, REF_E_D, REF_E_L, &
      REF_Q, REF_F, REF_M, REF_LOAD, IERR, MESSAGE)
    OK = IERR == SAIL_OK
    DO ITERATION = 1, 1000
      CALL ComputeSailForce(WIND, 10.0_DP, &
        0.3_DP, 0.2_DP, R, E_D, E_L, &
        Q, F, M, LOAD, IERR, MESSAGE)
      OK = OK .AND. IERR == SAIL_OK .AND. &
        NearVector(E_D, REF_E_D, VECTOR_TOL) .AND. &
        NearVector(E_L, REF_E_L, VECTOR_TOL) .AND. &
        NearScalar(Q, REF_Q, FORCE_TOL) .AND. &
        NearVector(F, REF_F, FORCE_TOL) .AND. &
        NearVector(M, REF_M, FORCE_TOL) .AND. &
        NearVector(LOAD, REF_LOAD, FORCE_TOL)
      IF (.NOT. OK) THEN
        CALL ReportFailureDetails(IERR, SAIL_OK, E_D, REF_E_D, E_L, &
          REF_E_L, Q, REF_Q, F, REF_F, M, REF_M, LOAD, REF_LOAD, MESSAGE)
        EXIT
      END IF
    END DO
    CALL RecordResult('1000 repeated calls without state pollution', OK)
  END SUBROUTINE TestRepeatedCalls


  SUBROUTINE TestInputsUnchanged
    REAL(DP) :: WIND(3), SAVED_WIND(3), SPEED, SAVED_SPEED
    REAL(DP) :: CL_VALUE, SAVED_CL, CD_VALUE, SAVED_CD, R(3), SAVED_R(3)
    REAL(DP) :: E_D(3), E_L(3), Q, F(3), M(3), LOAD(N_DOF)
    INTEGER :: IERR
    CHARACTER(LEN=512) :: MESSAGE
    LOGICAL :: OK

    WIND = (/ 6.0_DP, 8.0_DP, 0.0_DP /)
    SPEED = 10.0_DP
    CL_VALUE = 0.3_DP
    CD_VALUE = 0.2_DP
    R = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    SAVED_WIND = WIND
    SAVED_SPEED = SPEED
    SAVED_CL = CL_VALUE
    SAVED_CD = CD_VALUE
    SAVED_R = R
    CALL ComputeSailForce(WIND, SPEED, CL_VALUE, CD_VALUE, R, E_D, E_L, &
      Q, F, M, LOAD, IERR, MESSAGE)
    OK = IERR == SAIL_OK .AND. NearVector(WIND, SAVED_WIND, 0.0_DP) .AND. &
      NearScalar(SPEED, SAVED_SPEED, 0.0_DP) .AND. &
      NearScalar(CL_VALUE, SAVED_CL, 0.0_DP) .AND. &
      NearScalar(CD_VALUE, SAVED_CD, 0.0_DP) .AND. &
      NearVector(R, SAVED_R, 0.0_DP)
    IF (.NOT. OK) THEN
      CALL ReportFailureDetails(IERR, SAIL_OK, E_D, E_D, E_L, E_L, &
        Q, Q, F, F, M, M, LOAD, LOAD, MESSAGE)
      WRITE(*, '(A,3ES16.7)') '  actual/saved wind (actual): ', WIND
      WRITE(*, '(A,3ES16.7)') '  actual/saved wind (saved):  ', SAVED_WIND
      WRITE(*, '(A,4ES16.7)') '  speed, saved, CL, saved CL: ', &
        SPEED, SAVED_SPEED, CL_VALUE, SAVED_CL
      WRITE(*, '(A,2ES16.7)') '  CD, saved CD: ', CD_VALUE, SAVED_CD
      WRITE(*, '(A,3ES16.7)') '  actual/saved R (actual): ', R
      WRITE(*, '(A,3ES16.7)') '  actual/saved R (saved):  ', SAVED_R
    END IF
    CALL RecordResult('input arguments remain unchanged', OK)
  END SUBROUTINE TestInputsUnchanged

END PROGRAM TEST_SAILFORCE
