PROGRAM TEST_SAILRELWIND
  USE SAILPARAM_MOD
  USE SAILRELWIND_MOD
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE

  REAL(DP), PARAMETER :: TEST_TOL = 1.0E-12_DP
  INTEGER :: N_PASS
  INTEGER :: N_FAIL

  N_PASS = 0
  N_FAIL = 0

  CALL TestStandardCases()
  CALL TestInvalidInputs()
  CALL TestRepeatedCalls()
  CALL TestInputsUnchanged()

  WRITE(*, '(A)') '========================================'
  WRITE(*, '(A)') 'SAILRELWIND TEST SUMMARY'
  WRITE(*, '(A,I0)') 'PASS: ', N_PASS
  WRITE(*, '(A,I0)') 'FAIL: ', N_FAIL
  WRITE(*, '(A)') '========================================'

  IF (N_FAIL > 0) STOP 1

CONTAINS

  SUBROUTINE TestStandardCases()
    REAL(DP) :: Z(3)

    Z = 0.0_DP
    CALL RunCase('15.1 Stationary vessel', &
      (/ 10.0_DP, 0.0_DP, 0.0_DP /), Z, Z, &
      (/ 2.0_DP, 1.0_DP, 3.0_DP /), Z, &
      (/ 10.0_DP, 0.0_DP, 0.0_DP /), &
      (/ 10.0_DP, 0.0_DP, 0.0_DP /), 10.0_DP, SAIL_OK)

    CALL RunCase('15.2 Pure translation', &
      (/ 10.0_DP, 4.0_DP, 1.0_DP /), &
      (/ 3.0_DP, -2.0_DP, 0.5_DP /), Z, &
      (/ 9.0_DP, -4.0_DP, 2.0_DP /), &
      (/ 3.0_DP, -2.0_DP, 0.5_DP /), &
      (/ 7.0_DP, 6.0_DP, 0.5_DP /), &
      (/ 7.0_DP, 6.0_DP, 0.0_DP /), SQRT(85.0_DP), SAIL_OK)

    CALL RunCase('15.3 Pure yaw rate', Z, Z, &
      (/ 0.0_DP, 0.0_DP, 1.0_DP /), &
      (/ 2.0_DP, 0.0_DP, 0.0_DP /), &
      (/ 0.0_DP, 2.0_DP, 0.0_DP /), &
      (/ 0.0_DP, -2.0_DP, 0.0_DP /), &
      (/ 0.0_DP, -2.0_DP, 0.0_DP /), 2.0_DP, SAIL_OK)

    CALL RunCase('15.4 Pure roll rate', Z, Z, &
      (/ 2.0_DP, 0.0_DP, 0.0_DP /), &
      (/ 0.0_DP, 3.0_DP, 4.0_DP /), &
      (/ 0.0_DP, -8.0_DP, 6.0_DP /), &
      (/ 0.0_DP, 8.0_DP, -6.0_DP /), &
      (/ 0.0_DP, 8.0_DP, 0.0_DP /), 8.0_DP, SAIL_OK)

    CALL RunCase('15.5 Pure pitch rate', Z, Z, &
      (/ 0.0_DP, 2.0_DP, 0.0_DP /), &
      (/ 3.0_DP, 0.0_DP, 4.0_DP /), &
      (/ 8.0_DP, 0.0_DP, -6.0_DP /), &
      (/ -8.0_DP, 0.0_DP, 6.0_DP /), &
      (/ -8.0_DP, 0.0_DP, 0.0_DP /), 8.0_DP, SAIL_OK)

    CALL RunCase('15.6 General cross product', &
      (/ 10.0_DP, 20.0_DP, 30.0_DP /), &
      (/ 1.0_DP, 2.0_DP, 3.0_DP /), &
      (/ 1.0_DP, 2.0_DP, 3.0_DP /), &
      (/ 4.0_DP, 5.0_DP, 6.0_DP /), &
      (/ -2.0_DP, 8.0_DP, 0.0_DP /), &
      (/ 12.0_DP, 12.0_DP, 30.0_DP /), &
      (/ 12.0_DP, 12.0_DP, 0.0_DP /), SQRT(288.0_DP), SAIL_OK)

    CALL RunCase('15.7 Sail point at CG', &
      (/ 7.0_DP, 1.0_DP, -2.0_DP /), &
      (/ 2.0_DP, -3.0_DP, 1.0_DP /), &
      (/ 4.0_DP, 5.0_DP, 6.0_DP /), Z, &
      (/ 2.0_DP, -3.0_DP, 1.0_DP /), &
      (/ 5.0_DP, 4.0_DP, -3.0_DP /), &
      (/ 5.0_DP, 4.0_DP, 0.0_DP /), SQRT(41.0_DP), SAIL_OK)

    CALL RunCase('15.8 Vertical relative wind only', &
      (/ 0.0_DP, 0.0_DP, 5.0_DP /), Z, Z, Z, Z, &
      (/ 0.0_DP, 0.0_DP, 5.0_DP /), Z, 0.0_DP, &
      SAIL_ERR_LOW_WIND_SPEED)

    CALL RunCase('15.9 Zero relative wind', &
      (/ 4.0_DP, -2.0_DP, 1.0_DP /), &
      (/ 4.0_DP, -2.0_DP, 1.0_DP /), Z, &
      (/ 8.0_DP, 2.0_DP, -1.0_DP /), &
      (/ 4.0_DP, -2.0_DP, 1.0_DP /), Z, Z, 0.0_DP, &
      SAIL_ERR_LOW_WIND_SPEED)

    CALL RunCase('15.10 Below VREL_MIN', &
      (/ 0.5_DP * VREL_MIN, 0.0_DP, 0.0_DP /), Z, Z, Z, Z, &
      (/ 0.5_DP * VREL_MIN, 0.0_DP, 0.0_DP /), &
      (/ 0.5_DP * VREL_MIN, 0.0_DP, 0.0_DP /), &
      0.5_DP * VREL_MIN, SAIL_ERR_LOW_WIND_SPEED)

    CALL RunCase('15.11 Equal to VREL_MIN', &
      (/ VREL_MIN, 0.0_DP, 0.0_DP /), Z, Z, Z, Z, &
      (/ VREL_MIN, 0.0_DP, 0.0_DP /), &
      (/ VREL_MIN, 0.0_DP, 0.0_DP /), VREL_MIN, SAIL_OK)

    CALL RunCase('15.12 Above VREL_MIN', &
      (/ 2.0_DP * VREL_MIN, 0.0_DP, 0.0_DP /), Z, Z, Z, Z, &
      (/ 2.0_DP * VREL_MIN, 0.0_DP, 0.0_DP /), &
      (/ 2.0_DP * VREL_MIN, 0.0_DP, 0.0_DP /), &
      2.0_DP * VREL_MIN, SAIL_OK)
  END SUBROUTINE TestStandardCases

  SUBROUTINE TestInvalidInputs()
    REAL(DP) :: Z(3)
    REAL(DP) :: INPUT(3)

    Z = 0.0_DP
    INPUT = Z
    INPUT(IDX_Y) = IEEE_VALUE(0.0_DP, IEEE_QUIET_NAN)
    CALL RunCase('15.13 NaN wind input', INPUT, Z, Z, Z, &
      Z, Z, Z, 0.0_DP, SAIL_ERR_INVALID_INPUT)

    INPUT = Z
    INPUT(IDX_Z) = IEEE_VALUE(0.0_DP, IEEE_POSITIVE_INF)
    CALL RunCase('15.14 Inf angular-rate input', Z, Z, INPUT, Z, &
      Z, Z, Z, 0.0_DP, SAIL_ERR_INVALID_INPUT)
  END SUBROUTINE TestInvalidInputs

  SUBROUTINE TestRepeatedCalls()
    REAL(DP) :: V_WIND(3)
    REAL(DP) :: V_CG(3)
    REAL(DP) :: OMEGA(3)
    REAL(DP) :: R_SAIL(3)
    REAL(DP) :: EXPECTED_POINT(3)
    REAL(DP) :: EXPECTED_REL(3)
    REAL(DP) :: EXPECTED_REL_H(3)
    REAL(DP) :: V_POINT(3)
    REAL(DP) :: V_REL(3)
    REAL(DP) :: V_REL_H(3)
    REAL(DP) :: V_MAG
    INTEGER :: I
    INTEGER :: IERR
    LOGICAL :: OK
    CHARACTER(LEN=512) :: MESSAGE

    V_WIND = (/ 10.0_DP, 20.0_DP, 30.0_DP /)
    V_CG = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    OMEGA = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    R_SAIL = (/ 4.0_DP, 5.0_DP, 6.0_DP /)
    EXPECTED_POINT = (/ -2.0_DP, 8.0_DP, 0.0_DP /)
    EXPECTED_REL = (/ 12.0_DP, 12.0_DP, 30.0_DP /)
    EXPECTED_REL_H = (/ 12.0_DP, 12.0_DP, 0.0_DP /)
    OK = .TRUE.
    DO I = 1, 1000
      CALL ComputeSailRelativeWind( &
        V_WIND, V_CG, OMEGA, R_SAIL, &
        V_POINT, V_REL, V_REL_H, V_MAG, IERR, MESSAGE)
      OK = OK .AND. IERR == SAIL_OK
      OK = OK .AND. VectorClose(V_POINT, EXPECTED_POINT)
      OK = OK .AND. VectorClose(V_REL, EXPECTED_REL)
      OK = OK .AND. VectorClose(V_REL_H, EXPECTED_REL_H)
      OK = OK .AND. ABS(V_MAG - SQRT(288.0_DP)) <= TEST_TOL
    END DO
    CALL ReportSimple('15.15 Repeated calls (1000)', OK, &
      'All 1000 calls returned identical expected results.')
  END SUBROUTINE TestRepeatedCalls

  SUBROUTINE TestInputsUnchanged()
    REAL(DP) :: V_WIND(3)
    REAL(DP) :: V_CG(3)
    REAL(DP) :: OMEGA(3)
    REAL(DP) :: R_SAIL(3)
    REAL(DP) :: V_WIND_COPY(3)
    REAL(DP) :: V_CG_COPY(3)
    REAL(DP) :: OMEGA_COPY(3)
    REAL(DP) :: R_SAIL_COPY(3)
    REAL(DP) :: V_POINT(3)
    REAL(DP) :: V_REL(3)
    REAL(DP) :: V_REL_H(3)
    REAL(DP) :: V_MAG
    INTEGER :: IERR
    LOGICAL :: OK
    CHARACTER(LEN=512) :: MESSAGE

    V_WIND = (/ 10.0_DP, 20.0_DP, 30.0_DP /)
    V_CG = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    OMEGA = (/ 1.0_DP, 2.0_DP, 3.0_DP /)
    R_SAIL = (/ 4.0_DP, 5.0_DP, 6.0_DP /)
    V_WIND_COPY = V_WIND
    V_CG_COPY = V_CG
    OMEGA_COPY = OMEGA
    R_SAIL_COPY = R_SAIL

    CALL ComputeSailRelativeWind(V_WIND, V_CG, OMEGA, R_SAIL, &
      V_POINT, V_REL, V_REL_H, V_MAG, IERR, MESSAGE)
    OK = VectorClose(V_WIND, V_WIND_COPY) .AND. &
      VectorClose(V_CG, V_CG_COPY) .AND. &
      VectorClose(OMEGA, OMEGA_COPY) .AND. &
      VectorClose(R_SAIL, R_SAIL_COPY)
    CALL ReportSimple('15.16 Input arrays unchanged', OK, &
      'All four input arrays match their saved copies.')
  END SUBROUTINE TestInputsUnchanged

  SUBROUTINE RunCase(NAME, V_WIND, V_CG, OMEGA, R_SAIL, &
      EXPECTED_POINT, EXPECTED_REL, EXPECTED_REL_H, &
      EXPECTED_MAG, EXPECTED_IERR)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    REAL(DP), INTENT(IN), CONTIGUOUS :: V_WIND(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: V_CG(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: OMEGA(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: R_SAIL(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: EXPECTED_POINT(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: EXPECTED_REL(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: EXPECTED_REL_H(:)
    REAL(DP), INTENT(IN) :: EXPECTED_MAG
    INTEGER, INTENT(IN) :: EXPECTED_IERR
    REAL(DP) :: ACTUAL_POINT(3)
    REAL(DP) :: ACTUAL_REL(3)
    REAL(DP) :: ACTUAL_REL_H(3)
    REAL(DP) :: ACTUAL_MAG
    INTEGER :: ACTUAL_IERR
    LOGICAL :: OK
    CHARACTER(LEN=512) :: MESSAGE

    CALL ComputeSailRelativeWind(V_WIND, V_CG, OMEGA, R_SAIL, &
      ACTUAL_POINT, ACTUAL_REL, ACTUAL_REL_H, ACTUAL_MAG, &
      ACTUAL_IERR, MESSAGE)
    OK = ACTUAL_IERR == EXPECTED_IERR .AND. &
      VectorClose(ACTUAL_POINT, EXPECTED_POINT) .AND. &
      VectorClose(ACTUAL_REL, EXPECTED_REL) .AND. &
      VectorClose(ACTUAL_REL_H, EXPECTED_REL_H) .AND. &
      ABS(ACTUAL_MAG - EXPECTED_MAG) <= TEST_TOL

    IF (OK) THEN
      N_PASS = N_PASS + 1
      WRITE(*, '(A,A)') '[PASS] ', TRIM(NAME)
      WRITE(*, '(A,I0,A,3(ES16.8,1X),A,ES16.8)') &
        '  IERR=', ACTUAL_IERR, ', V_REL=', ACTUAL_REL, &
        ', V_REL_H_MAG=', ACTUAL_MAG
    ELSE
      N_FAIL = N_FAIL + 1
      WRITE(*, '(A,A)') '[FAIL] ', TRIM(NAME)
      WRITE(*, '(A,I0,A,I0)') '  actual IERR=', ACTUAL_IERR, &
        ', expected IERR=', EXPECTED_IERR
      WRITE(*, '(A,3(ES16.8,1X))') '  actual point=', ACTUAL_POINT
      WRITE(*, '(A,3(ES16.8,1X))') '  expected point=', EXPECTED_POINT
      WRITE(*, '(A,3(ES16.8,1X))') '  actual V_REL=', ACTUAL_REL
      WRITE(*, '(A,3(ES16.8,1X))') '  expected V_REL=', EXPECTED_REL
      WRITE(*, '(A,3(ES16.8,1X))') '  actual V_REL_H=', ACTUAL_REL_H
      WRITE(*, '(A,3(ES16.8,1X))') '  expected V_REL_H=', EXPECTED_REL_H
      WRITE(*, '(A,ES16.8)') '  actual V_REL_H_MAG=', ACTUAL_MAG
      WRITE(*, '(A,ES16.8)') '  expected V_REL_H_MAG=', EXPECTED_MAG
      WRITE(*, '(A,A)') '  MESSAGE=', TRIM(MESSAGE)
    END IF
  END SUBROUTINE RunCase

  LOGICAL FUNCTION VectorClose(A, B)
    REAL(DP), INTENT(IN), CONTIGUOUS :: A(:)
    REAL(DP), INTENT(IN), CONTIGUOUS :: B(:)

    VectorClose = ALL(ABS(A - B) <= TEST_TOL)
  END FUNCTION VectorClose

  SUBROUTINE ReportSimple(NAME, OK, DETAIL)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    LOGICAL, INTENT(IN) :: OK
    CHARACTER(LEN=*), INTENT(IN) :: DETAIL

    IF (OK) THEN
      N_PASS = N_PASS + 1
      WRITE(*, '(A,A)') '[PASS] ', TRIM(NAME)
    ELSE
      N_FAIL = N_FAIL + 1
      WRITE(*, '(A,A)') '[FAIL] ', TRIM(NAME)
    END IF
    WRITE(*, '(A,A)') '  ', TRIM(DETAIL)
  END SUBROUTINE ReportSimple

END PROGRAM TEST_SAILRELWIND
