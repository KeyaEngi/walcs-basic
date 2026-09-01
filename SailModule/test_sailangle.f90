PROGRAM TEST_SAILANGLE
  USE SAILPARAM_MOD
  USE SAILANGLE_MOD
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE

  REAL(DP), PARAMETER :: TEST_TOL = 1.0E-10_DP
  REAL(DP), PARAMETER :: VECTOR_TOL = 1.0E-12_DP
  INTEGER :: N_PASS
  INTEGER :: N_FAIL

  N_PASS = 0
  N_FAIL = 0

  CALL TestZeroAngle
  CALL TestPositiveTen
  CALL TestNegativeTen
  CALL TestPeriodicity
  CALL TestPositive370
  CALL TestNegative370
  CALL TestWindPositiveX
  CALL TestWindNegativeX
  CALL TestWindNegativeY
  CALL TestObliqueWind
  CALL TestLowWind
  CALL TestMinimumWind
  CALL TestMagnitudeMismatch
  CALL TestNonzeroZ
  CALL TestNanDelta
  CALL TestInfiniteVector
  CALL TestNanMagnitude
  CALL TestRepeatedCalls
  CALL TestInputsUnchanged
  CALL TestAngleScan

  WRITE(*, '(A)') '========================================'
  WRITE(*, '(A)') 'SAILANGLE TEST SUMMARY'
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


  SUBROUTINE Evaluate(VX, VY, VZ, MAGNITUDE, DELTA, NORMALIZED, CHORD, &
      UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    REAL(DP), INTENT(IN) :: VX, VY, VZ
    REAL(DP), INTENT(IN) :: MAGNITUDE
    REAL(DP), INTENT(IN) :: DELTA
    REAL(DP), INTENT(OUT) :: NORMALIZED
    REAL(DP), INTENT(OUT) :: CHORD(3)
    REAL(DP), INTENT(OUT) :: UPSTREAM(3)
    REAL(DP), INTENT(OUT) :: RAW_ALPHA
    REAL(DP), INTENT(OUT) :: DB_ALPHA
    INTEGER, INTENT(OUT) :: IERR
    REAL(DP) :: VECTOR(3)
    CHARACTER(LEN=512) :: MESSAGE

    VECTOR = (/ VX, VY, VZ /)
    CALL ComputeSailAngle(VECTOR, MAGNITUDE, DELTA, NORMALIZED, CHORD, &
      UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR, MESSAGE)
  END SUBROUTINE Evaluate


  LOGICAL FUNCTION NearScalar(A, B, TOLERANCE)
    REAL(DP), INTENT(IN) :: A
    REAL(DP), INTENT(IN) :: B
    REAL(DP), INTENT(IN) :: TOLERANCE

    NearScalar = ABS(A - B) <= TOLERANCE
  END FUNCTION NearScalar


  LOGICAL FUNCTION NearVector(A, B, TOLERANCE)
    REAL(DP), INTENT(IN) :: A(3)
    REAL(DP), INTENT(IN) :: B(3)
    REAL(DP), INTENT(IN) :: TOLERANCE

    NearVector = MAXVAL(ABS(A - B)) <= TOLERANCE
  END FUNCTION NearVector


  LOGICAL FUNCTION NearVectorComponents(A, X, Y, Z, TOLERANCE)
    REAL(DP), INTENT(IN) :: A(3)
    REAL(DP), INTENT(IN) :: X, Y, Z
    REAL(DP), INTENT(IN) :: TOLERANCE

    NearVectorComponents = ABS(A(IDX_X) - X) <= TOLERANCE .AND. &
      ABS(A(IDX_Y) - Y) <= TOLERANCE .AND. &
      ABS(A(IDX_Z) - Z) <= TOLERANCE
  END FUNCTION NearVectorComponents


  LOGICAL FUNCTION OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, &
      RAW_ALPHA, DB_ALPHA)
    REAL(DP), INTENT(IN) :: NORMALIZED
    REAL(DP), INTENT(IN) :: CHORD(3)
    REAL(DP), INTENT(IN) :: UPSTREAM(3)
    REAL(DP), INTENT(IN) :: RAW_ALPHA
    REAL(DP), INTENT(IN) :: DB_ALPHA

    OutputsAreZero = ABS(NORMALIZED) <= TEST_TOL .AND. &
      MAXVAL(ABS(CHORD)) <= TEST_TOL .AND. &
      MAXVAL(ABS(UPSTREAM)) <= TEST_TOL .AND. &
      ABS(RAW_ALPHA) <= TEST_TOL .AND. ABS(DB_ALPHA) <= TEST_TOL
  END FUNCTION OutputsAreZero


  SUBROUTINE TestZeroAngle
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearScalar(NORMALIZED, 0.0_DP, TEST_TOL) &
      .AND. NearVectorComponents(CHORD, 0.0_DP, -1.0_DP, 0.0_DP, VECTOR_TOL) &
      .AND. NearVectorComponents(UPSTREAM, 0.0_DP, -1.0_DP, 0.0_DP, VECTOR_TOL) &
      .AND. NearScalar(RAW_ALPHA, 0.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 0.0_DP, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  zero: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('zero-angle baseline (+y wind)', OK)
  END SUBROUTINE TestZeroAngle


  SUBROUTINE TestPositiveTen
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    REAL(DP) :: EXPECTED_CHORD(3)
    INTEGER :: IERR
    LOGICAL :: OK

    EXPECTED_CHORD = (/ -SIN(10.0_DP * DEG_TO_RAD), &
      -COS(10.0_DP * DEG_TO_RAD), 0.0_DP /)
    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, 10.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearScalar(NORMALIZED, 10.0_DP, TEST_TOL) &
      .AND. NearVector(CHORD, EXPECTED_CHORD, VECTOR_TOL) &
      .AND. NearScalar(RAW_ALPHA, 10.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 10.0_DP, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  +10: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('positive 10-degree angle', OK)
  END SUBROUTINE TestPositiveTen


  SUBROUTINE TestNegativeTen
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, -10.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearScalar(NORMALIZED, 170.0_DP, TEST_TOL) &
      .AND. NearScalar(RAW_ALPHA, 170.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 170.0_DP, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  -10: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('negative 10-degree periodic mapping', OK)
  END SUBROUTINE TestNegativeTen


  SUBROUTINE TestPeriodicity
    REAL(DP) :: N0, N180, N360
    REAL(DP) :: C0(3), C180(3), C360(3), U(3), RAW, A0, A180, A360
    INTEGER :: IERR0, IERR180, IERR360
    LOGICAL :: OK

    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      N0, C0, U, RAW, A0, IERR0)
    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, 180.0_DP, &
      N180, C180, U, RAW, A180, IERR180)
    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, 360.0_DP, &
      N360, C360, U, RAW, A360, IERR360)
    OK = IERR0 == SAIL_OK .AND. IERR180 == SAIL_OK .AND. &
      IERR360 == SAIL_OK .AND. NearScalar(N0, 0.0_DP, TEST_TOL) .AND. &
      NearScalar(N180, 0.0_DP, TEST_TOL) .AND. &
      NearScalar(N360, 0.0_DP, TEST_TOL) .AND. &
      NearVector(C0, C180, VECTOR_TOL) .AND. &
      NearVector(C0, C360, VECTOR_TOL) .AND. &
      NearScalar(A0, A180, TEST_TOL) .AND. NearScalar(A0, A360, TEST_TOL)
    CALL RecordResult('0/180/360-degree sail periodicity', OK)
  END SUBROUTINE TestPeriodicity


  SUBROUTINE TestPositive370
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, 370.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearScalar(NORMALIZED, 10.0_DP, TEST_TOL) &
      .AND. NearScalar(RAW_ALPHA, 10.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 10.0_DP, TEST_TOL)
    CALL RecordResult('positive 370-degree normalization', OK)
  END SUBROUTINE TestPositive370


  SUBROUTINE TestNegative370
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, -370.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearScalar(NORMALIZED, 170.0_DP, TEST_TOL) &
      .AND. NearScalar(RAW_ALPHA, 170.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 170.0_DP, TEST_TOL)
    CALL RecordResult('negative 370-degree normalization', OK)
  END SUBROUTINE TestNegative370


  SUBROUTINE TestWindPositiveX
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(10.0_DP, 0.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. &
      NearVectorComponents(UPSTREAM, -1.0_DP, 0.0_DP, 0.0_DP, VECTOR_TOL) &
      .AND. NearScalar(RAW_ALPHA, -90.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 90.0_DP, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  +x wind: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('relative wind along +x_b', OK)
  END SUBROUTINE TestWindPositiveX


  SUBROUTINE TestWindNegativeX
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(-10.0_DP, 0.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. &
      NearVectorComponents(UPSTREAM, 1.0_DP, 0.0_DP, 0.0_DP, VECTOR_TOL) &
      .AND. NearScalar(RAW_ALPHA, 90.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 90.0_DP, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  -x wind: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('relative wind along -x_b', OK)
  END SUBROUTINE TestWindNegativeX


  SUBROUTINE TestWindNegativeY
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(0.0_DP, -10.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. &
      NearVectorComponents(UPSTREAM, 0.0_DP, 1.0_DP, 0.0_DP, VECTOR_TOL) &
      .AND. NearScalar(ABS(RAW_ALPHA), 180.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 0.0_DP, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  -y wind: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('relative wind along -y_b', OK)
  END SUBROUTINE TestWindNegativeY


  SUBROUTINE TestObliqueWind
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    REAL(DP) :: EXPECTED_CHORD(3), EXPECTED_UPSTREAM(3)
    REAL(DP) :: EXPECTED_CROSS, EXPECTED_DOT, EXPECTED_RAW, EXPECTED_DB
    INTEGER :: IERR
    LOGICAL :: OK

    EXPECTED_CHORD = (/ -SIN(25.0_DP * DEG_TO_RAD), &
      -COS(25.0_DP * DEG_TO_RAD), 0.0_DP /)
    EXPECTED_UPSTREAM = (/ -0.6_DP, -0.8_DP, 0.0_DP /)
    EXPECTED_CROSS = EXPECTED_CHORD(IDX_X) * EXPECTED_UPSTREAM(IDX_Y) - &
      EXPECTED_CHORD(IDX_Y) * EXPECTED_UPSTREAM(IDX_X)
    EXPECTED_DOT = EXPECTED_CHORD(IDX_X) * EXPECTED_UPSTREAM(IDX_X) + &
      EXPECTED_CHORD(IDX_Y) * EXPECTED_UPSTREAM(IDX_Y)
    EXPECTED_RAW = ATAN2(EXPECTED_CROSS, EXPECTED_DOT) * RAD_TO_DEG
    EXPECTED_DB = MODULO(EXPECTED_RAW, SAIL_ANGLE_PERIOD_DEG)

    CALL Evaluate(6.0_DP, 8.0_DP, 0.0_DP, 10.0_DP, 25.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearVector(CHORD, EXPECTED_CHORD, VECTOR_TOL) &
      .AND. NearVector(UPSTREAM, EXPECTED_UPSTREAM, VECTOR_TOL) &
      .AND. NearScalar(RAW_ALPHA, EXPECTED_RAW, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, EXPECTED_DB, TEST_TOL)
    WRITE(*, '(A,2F12.6)') '  oblique: alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('independently calculated oblique wind', OK)
  END SUBROUTINE TestObliqueWind


  SUBROUTINE TestLowWind
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(0.5_DP * VREL_MIN, 0.0_DP, 0.0_DP, &
      0.5_DP * VREL_MIN, 30.0_DP, NORMALIZED, CHORD, UPSTREAM, &
      RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_ERR_LOW_WIND_SPEED .AND. &
      OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA)
    WRITE(*, '(A,I0)') '  low wind error code = ', IERR
    CALL RecordResult('low-wind zero-output handling', OK)
  END SUBROUTINE TestLowWind


  SUBROUTINE TestMinimumWind
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(VREL_MIN, 0.0_DP, 0.0_DP, VREL_MIN, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_OK .AND. NearScalar(RAW_ALPHA, -90.0_DP, TEST_TOL) &
      .AND. NearScalar(DB_ALPHA, 90.0_DP, TEST_TOL)
    WRITE(*, '(A,I0,A,2F12.6)') '  VREL_MIN code = ', IERR, &
      ', alpha_raw, alpha_db = ', RAW_ALPHA, DB_ALPHA
    CALL RecordResult('wind speed equal to VREL_MIN', OK)
  END SUBROUTINE TestMinimumWind


  SUBROUTINE TestMagnitudeMismatch
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(3.0_DP, 4.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_ERR_INVALID_INPUT .AND. &
      OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA)
    CALL RecordResult('scalar/vector wind-speed mismatch', OK)
  END SUBROUTINE TestMagnitudeMismatch


  SUBROUTINE TestNonzeroZ
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    CALL Evaluate(3.0_DP, 4.0_DP, 1.0_DP, 5.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_ERR_INVALID_INPUT .AND. &
      OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA)
    CALL RecordResult('nonzero horizontal-wind z component', OK)
  END SUBROUTINE TestNonzeroZ


  SUBROUTINE TestNanDelta
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    REAL(DP) :: NAN_VALUE
    INTEGER :: IERR
    LOGICAL :: OK

    NAN_VALUE = IEEE_VALUE(0.0_DP, IEEE_QUIET_NAN)
    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, NAN_VALUE, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_ERR_INVALID_INPUT .AND. &
      OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA)
    CALL RecordResult('NaN sail-angle rejection', OK)
  END SUBROUTINE TestNanDelta


  SUBROUTINE TestInfiniteVector
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    REAL(DP) :: INF_VALUE
    INTEGER :: IERR
    LOGICAL :: OK

    INF_VALUE = IEEE_VALUE(0.0_DP, IEEE_POSITIVE_INF)
    CALL Evaluate(INF_VALUE, 10.0_DP, 0.0_DP, 10.0_DP, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_ERR_INVALID_INPUT .AND. &
      OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA)
    CALL RecordResult('infinite wind-vector rejection', OK)
  END SUBROUTINE TestInfiniteVector


  SUBROUTINE TestNanMagnitude
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    REAL(DP) :: NAN_VALUE
    INTEGER :: IERR
    LOGICAL :: OK

    NAN_VALUE = IEEE_VALUE(0.0_DP, IEEE_QUIET_NAN)
    CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, NAN_VALUE, 0.0_DP, &
      NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
    OK = IERR == SAIL_ERR_INVALID_INPUT .AND. &
      OutputsAreZero(NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA)
    CALL RecordResult('NaN wind-magnitude rejection', OK)
  END SUBROUTINE TestNanMagnitude


  SUBROUTINE TestRepeatedCalls
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    REAL(DP) :: REF_NORMALIZED, REF_CHORD(3), REF_UPSTREAM(3)
    REAL(DP) :: REF_RAW, REF_DB
    INTEGER :: IERR, ITERATION
    LOGICAL :: OK

    CALL Evaluate(6.0_DP, 8.0_DP, 0.0_DP, 10.0_DP, 25.0_DP, &
      REF_NORMALIZED, REF_CHORD, REF_UPSTREAM, REF_RAW, REF_DB, IERR)
    OK = IERR == SAIL_OK
    DO ITERATION = 1, 1000
      CALL Evaluate(6.0_DP, 8.0_DP, 0.0_DP, 10.0_DP, 25.0_DP, &
        NORMALIZED, CHORD, UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR)
      OK = OK .AND. IERR == SAIL_OK &
        .AND. NearScalar(NORMALIZED, REF_NORMALIZED, TEST_TOL) &
        .AND. NearVector(CHORD, REF_CHORD, TEST_TOL) &
        .AND. NearVector(UPSTREAM, REF_UPSTREAM, TEST_TOL) &
        .AND. NearScalar(RAW_ALPHA, REF_RAW, TEST_TOL) &
        .AND. NearScalar(DB_ALPHA, REF_DB, TEST_TOL)
    END DO
    CALL RecordResult('1000 repeated calls without state pollution', OK)
  END SUBROUTINE TestRepeatedCalls


  SUBROUTINE TestInputsUnchanged
    REAL(DP) :: VECTOR(3), VECTOR_SAVED(3), MAGNITUDE, MAGNITUDE_SAVED
    REAL(DP) :: DELTA, DELTA_SAVED
    REAL(DP) :: NORMALIZED, CHORD(3), UPSTREAM(3), RAW_ALPHA, DB_ALPHA
    INTEGER :: IERR
    LOGICAL :: OK

    VECTOR = (/ 6.0_DP, 8.0_DP, 0.0_DP /)
    MAGNITUDE = 10.0_DP
    DELTA = 25.0_DP
    VECTOR_SAVED = VECTOR
    MAGNITUDE_SAVED = MAGNITUDE
    DELTA_SAVED = DELTA
    BLOCK
      CHARACTER(LEN=512) :: MESSAGE
      CALL ComputeSailAngle(VECTOR, MAGNITUDE, DELTA, NORMALIZED, CHORD, &
        UPSTREAM, RAW_ALPHA, DB_ALPHA, IERR, MESSAGE)
    END BLOCK
    OK = IERR == SAIL_OK .AND. NearVector(VECTOR, VECTOR_SAVED, 0.0_DP) &
      .AND. NearScalar(MAGNITUDE, MAGNITUDE_SAVED, 0.0_DP) &
      .AND. NearScalar(DELTA, DELTA_SAVED, 0.0_DP)
    CALL RecordResult('input arguments remain unchanged', OK)
  END SUBROUTINE TestInputsUnchanged


  SUBROUTINE TestAngleScan
    REAL(DP) :: DELTA, NORMALIZED1, NORMALIZED2
    REAL(DP) :: CHORD1(3), CHORD2(3), UPSTREAM1(3), UPSTREAM2(3)
    REAL(DP) :: RAW1, RAW2, DB1, DB2, CHORD_NORM, UPSTREAM_NORM
    INTEGER :: IERR1, IERR2, INDEX
    LOGICAL :: OK

    OK = .TRUE.
    DO INDEX = 0, 35
      DELTA = 10.0_DP * REAL(INDEX, DP)
      CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, DELTA, &
        NORMALIZED1, CHORD1, UPSTREAM1, RAW1, DB1, IERR1)
      CALL Evaluate(0.0_DP, 10.0_DP, 0.0_DP, 10.0_DP, &
        DELTA + 180.0_DP, NORMALIZED2, CHORD2, UPSTREAM2, RAW2, DB2, IERR2)
      CHORD_NORM = SQRT(CHORD1(IDX_X)**2 + CHORD1(IDX_Y)**2)
      UPSTREAM_NORM = SQRT(UPSTREAM1(IDX_X)**2 + UPSTREAM1(IDX_Y)**2)
      OK = OK .AND. IERR1 == SAIL_OK .AND. IERR2 == SAIL_OK &
        .AND. NORMALIZED1 >= 0.0_DP &
        .AND. NORMALIZED1 < SAIL_ANGLE_PERIOD_DEG &
        .AND. DB1 >= 0.0_DP .AND. DB1 < SAIL_ANGLE_PERIOD_DEG &
        .AND. NearScalar(CHORD_NORM, 1.0_DP, VECTOR_TOL) &
        .AND. NearScalar(UPSTREAM_NORM, 1.0_DP, VECTOR_TOL) &
        .AND. ALL(IEEE_IS_FINITE(CHORD1)) &
        .AND. ALL(IEEE_IS_FINITE(UPSTREAM1)) &
        .AND. IEEE_IS_FINITE(NORMALIZED1) .AND. IEEE_IS_FINITE(RAW1) &
        .AND. IEEE_IS_FINITE(DB1) .AND. NearScalar(DB1, DB2, TEST_TOL) &
        .AND. NearScalar(NORMALIZED1, NORMALIZED2, TEST_TOL) &
        .AND. NearVector(CHORD1, CHORD2, VECTOR_TOL) &
        .AND. NearVector(UPSTREAM1, UPSTREAM2, VECTOR_TOL) &
        .AND. IEEE_IS_FINITE(RAW2)
    END DO
    CALL RecordResult('0-to-350-degree scan and periodic consistency', OK)
  END SUBROUTINE TestAngleScan

END PROGRAM TEST_SAILANGLE
