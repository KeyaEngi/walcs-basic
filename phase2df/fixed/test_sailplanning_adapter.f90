PROGRAM test_sailplanning_adapter
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE SAILPARAM_MOD, ONLY: DP, PI, SAIL_OK, SAIL_ERR_LOW_WIND_SPEED
  USE SAILMODULE_API_MOD, ONLY: ComputeSailModuleLoads
  USE SAILPLANNING_ADAPTER_MOD, ONLY: InitializeSailPlanningAdapter, &
    ComputeSailPlanningDryRun, TransformSailLoadToPlanning, &
    FinalizeSailPlanningAdapter
  USE ArrayOperations, ONLY: VectorG2L
  IMPLICIT NONE

  INTEGER, PARAMETER :: NR = 6
  REAL(DP), PARAMETER :: TOL = 1.0E-12_DP
  CHARACTER(LEN=*), PARAMETER :: DATABASE_FILE = &
    'D:\GitHub\walcs-basic\SailModule\sail_database.dat'
  REAL(DP) :: Y_STATE(2*NR), WIND(3), R_SAIL(3)
  REAL(DP) :: V_CG(3), OMEGA(3), V_WIND(3), SPEED, ALPHA, CL, CD
  REAL(DP) :: FORCE(3), MOMENT(3), LOAD(6), EXPECTED(3)
  REAL(DP) :: VSP(3), VREL(3), VRELH(3), DELTA_N, CHORD(3), UP(3)
  REAL(DP) :: ALPHA_RAW, DRAG_DIR(3), LIFT_DIR(3), QDYN
  REAL(DP) :: FORCE_DIRECT(3), MOMENT_DIRECT(3), LOAD_DIRECT(6)
  REAL(DP) :: ANGLES(3), LOAD_BODY(6), LOAD_PLANNING(6), LOAD_EXPECTED(6)
  REAL(DP) :: LOAD_RECOVERED(6)
  INTEGER :: IERR, FAILURES
  CHARACTER(LEN=2048) :: MESSAGE

  FAILURES = 0
  CALL InitializeSailPlanningAdapter(DATABASE_FILE, IERR, MESSAGE)
  CALL Check(IERR == SAIL_OK, 'adapter initialization', FAILURES)

  Y_STATE = 0.0_DP
  WIND = [5.0_DP, 2.0_DP, 1.0_DP]
  R_SAIL = [1.25_DP, -0.4_DP, 2.0_DP]
  CALL RunAdapter(10.0_DP, Y_STATE, WIND, R_SAIL, 0.0_DP)
  CALL CheckVec(V_CG, [10.0_DP, 0.0_DP, 0.0_DP], TOL, &
    'zero-attitude VectorG2L and mean speed', FAILURES)
  CALL CheckVec(V_WIND, WIND, TOL, 'zero-attitude wind VectorG2L', FAILURES)

  Y_STATE = 0.0_DP
  Y_STATE(NR+1:NR+6) = [2.0_DP, 3.0_DP, -4.0_DP, 0.1_DP, -0.2_DP, 0.3_DP]
  CALL RunAdapter(10.0_DP, Y_STATE, WIND, R_SAIL, 0.0_DP)
  CALL CheckVec(V_CG, [12.0_DP, 3.0_DP, -4.0_DP], TOL, &
    'U0 once plus surge, with sway/heave mapping', FAILURES)
  CALL CheckVec(OMEGA, [0.1_DP, -0.2_DP, 0.3_DP], TOL, &
    'zero-attitude Euler rates', FAILURES)

  CALL CheckSingleRate(1, 0.37_DP, FAILURES)
  CALL CheckSingleRate(2, -0.41_DP, FAILURES)
  CALL CheckSingleRate(3, 0.29_DP, FAILURES)

  Y_STATE = 0.0_DP
  Y_STATE(4:6) = [0.2_DP, -0.3_DP, 0.4_DP]
  Y_STATE(NR+4:NR+6) = [0.11_DP, -0.07_DP, 0.19_DP]
  CALL RunAdapter(4.0_DP, Y_STATE, WIND, R_SAIL, 12.0_DP)
  EXPECTED(1) = COS(Y_STATE(5))*COS(Y_STATE(6))*Y_STATE(NR+4) + &
                SIN(Y_STATE(6))*Y_STATE(NR+5)
  EXPECTED(2) =-COS(Y_STATE(5))*SIN(Y_STATE(6))*Y_STATE(NR+4) + &
                COS(Y_STATE(6))*Y_STATE(NR+5)
  EXPECTED(3) = SIN(Y_STATE(5))*Y_STATE(NR+4) + Y_STATE(NR+6)
  CALL CheckVec(OMEGA, EXPECTED, TOL, 'nonzero-attitude independent omega matrix', FAILURES)

  Y_STATE = 0.0_DP
  Y_STATE(6) = 0.5_DP*PI
  WIND = [1.0_DP, 0.0_DP, 0.0_DP]
  CALL RunAdapter(4.0_DP, Y_STATE, WIND, R_SAIL, 0.0_DP)
  CALL CheckVec(V_WIND, [0.0_DP, -1.0_DP, 0.0_DP], TOL, &
    'yaw +90 deg global-to-body direction', FAILURES)

  ! End-to-end comparison also proves R_SAIL_BODY passes through unchanged:
  ! the direct API receives the same relative CG-to-CE body vector.
  Y_STATE = 0.0_DP
  Y_STATE(NR+1:NR+3) = [0.5_DP, -0.25_DP, 0.1_DP]
  WIND = [1.0_DP, 8.0_DP, 0.0_DP]
  R_SAIL = [1.25_DP, -0.4_DP, 2.0_DP]
  CALL RunAdapter(4.0_DP, Y_STATE, WIND, R_SAIL, 15.0_DP)
  CALL ComputeSailModuleLoads(V_WIND, V_CG, OMEGA, R_SAIL, 15.0_DP, &
    VSP, VREL, VRELH, SPEED, DELTA_N, CHORD, UP, ALPHA_RAW, ALPHA, &
    CL, CD, DRAG_DIR, LIFT_DIR, QDYN, FORCE_DIRECT, MOMENT_DIRECT, &
    LOAD_DIRECT, IERR, MESSAGE)
  CALL Check(IERR == SAIL_OK, 'direct SailModule end-to-end status', FAILURES)
  CALL CheckVec(FORCE, FORCE_DIRECT, TOL, 'adapter/direct force', FAILURES)
  CALL CheckVec(MOMENT, MOMENT_DIRECT, TOL, 'adapter/direct moment and R_SAIL', FAILURES)
  CALL CheckVec6(LOAD, LOAD_DIRECT, TOL, 'adapter/direct LOAD_6DOF', FAILURES)
  CALL Check(ALL(IEEE_IS_FINITE(LOAD)), 'finite LOAD_6DOF', FAILURES)

  Y_STATE = 0.0_DP
  WIND = [4.0_DP, 0.0_DP, 0.0_DP]
  CALL RunAdapter(4.0_DP, Y_STATE, WIND, R_SAIL, 0.0_DP)
  CALL Check(IERR == SAIL_ERR_LOW_WIND_SPEED, 'low wind is recognizable', FAILURES)
  CALL CheckVec6(LOAD, 0.0_DP*LOAD, TOL, 'low wind clears LOAD_6DOF', FAILURES)

  ! Phase 2A body-to-planning load transformation tests.
  ANGLES = 0.0_DP
  LOAD_BODY = [1.0_DP, -2.0_DP, 3.0_DP, -4.0_DP, 5.0_DP, -6.0_DP]
  CALL TransformSailLoadToPlanning(ANGLES, LOAD_BODY, LOAD_PLANNING, IERR, MESSAGE)
  CALL Check(IERR == SAIL_OK, 'load transform zero-attitude status', FAILURES)
  CALL CheckVec6(LOAD_PLANNING, LOAD_BODY, TOL, &
    'load transform zero attitude', FAILURES)

  ANGLES = [0.0_DP, 0.0_DP, 0.5_DP*PI]
  LOAD_BODY = [100.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP]
  CALL TransformSailLoadToPlanning(ANGLES, LOAD_BODY, LOAD_PLANNING, IERR, MESSAGE)
  LOAD_EXPECTED(1:3) = IndependentL2G(LOAD_BODY(1:3), ANGLES)
  LOAD_EXPECTED(4:6) = IndependentL2G(LOAD_BODY(4:6), ANGLES)
  CALL CheckVec6(LOAD_PLANNING, LOAD_EXPECTED, TOL, &
    'pure Fx at yaw +90 degrees', FAILURES)

  LOAD_BODY = [0.0_DP, 100.0_DP, 0.0_DP, 0.0_DP, 0.0_DP, 0.0_DP]
  CALL TransformSailLoadToPlanning(ANGLES, LOAD_BODY, LOAD_PLANNING, IERR, MESSAGE)
  LOAD_EXPECTED(1:3) = IndependentL2G(LOAD_BODY(1:3), ANGLES)
  LOAD_EXPECTED(4:6) = IndependentL2G(LOAD_BODY(4:6), ANGLES)
  CALL CheckVec6(LOAD_PLANNING, LOAD_EXPECTED, TOL, &
    'pure Fy at yaw +90 degrees', FAILURES)

  LOAD_BODY = [0.0_DP, 0.0_DP, 0.0_DP, 100.0_DP, 0.0_DP, 0.0_DP]
  CALL TransformSailLoadToPlanning(ANGLES, LOAD_BODY, LOAD_PLANNING, IERR, MESSAGE)
  LOAD_EXPECTED(1:3) = IndependentL2G(LOAD_BODY(1:3), ANGLES)
  LOAD_EXPECTED(4:6) = IndependentL2G(LOAD_BODY(4:6), ANGLES)
  CALL CheckVec6(LOAD_PLANNING, LOAD_EXPECTED, TOL, &
    'pure Mx uses force rotation rule', FAILURES)

  ANGLES = [0.2_DP, -0.3_DP, 0.4_DP]
  LOAD_BODY = [10.0_DP, -20.0_DP, 30.0_DP, -40.0_DP, 50.0_DP, -60.0_DP]
  CALL TransformSailLoadToPlanning(ANGLES, LOAD_BODY, LOAD_PLANNING, IERR, MESSAGE)
  LOAD_EXPECTED(1:3) = IndependentL2G(LOAD_BODY(1:3), ANGLES)
  LOAD_EXPECTED(4:6) = IndependentL2G(LOAD_BODY(4:6), ANGLES)
  CALL CheckVec6(LOAD_PLANNING, LOAD_EXPECTED, TOL, &
    'arbitrary-attitude independent matrix', FAILURES)

  LOAD_RECOVERED(1:3) = VectorG2L(LOAD_PLANNING(1:3), ANGLES)
  LOAD_RECOVERED(4:6) = VectorG2L(LOAD_PLANNING(4:6), ANGLES)
  CALL CheckVec6(LOAD_RECOVERED, LOAD_BODY, TOL, &
    'load rotation inverse recovery', FAILURES)

  CALL FinalizeSailPlanningAdapter(IERR, MESSAGE)
  CALL Check(IERR == SAIL_OK, 'adapter finalization', FAILURES)
  IF (FAILURES /= 0) THEN
    WRITE(*,'(A,I0)') '[FAIL] test_sailplanning_adapter failures: ', FAILURES
    ERROR STOP 1
  END IF
  WRITE(*,'(A)') '[PASS] test_sailplanning_adapter: all tests passed.'

CONTAINS

  SUBROUTINE RunAdapter(U0, STATE, VG, RS, DELTA)
    REAL(DP), INTENT(IN) :: U0, STATE(2*NR), VG(3), RS(3), DELTA
    CALL ComputeSailPlanningDryRun(U0, STATE, NR, 0.0_DP, VG, RS, DELTA, &
      V_CG, OMEGA, V_WIND, SPEED, ALPHA, CL, CD, FORCE, MOMENT, LOAD, &
      IERR, MESSAGE)
    CALL Check(IERR == SAIL_OK .OR. IERR == SAIL_ERR_LOW_WIND_SPEED, &
      'adapter call status', FAILURES)
  END SUBROUTINE RunAdapter

  SUBROUTINE CheckSingleRate(INDEX, RATE, NFAIL)
    INTEGER, INTENT(IN) :: INDEX
    REAL(DP), INTENT(IN) :: RATE
    INTEGER, INTENT(INOUT) :: NFAIL
    REAL(DP) :: TARGET(3)
    Y_STATE = 0.0_DP
    Y_STATE(NR+3+INDEX) = RATE
    WIND = [5.0_DP, 2.0_DP, 0.0_DP]
    CALL RunAdapter(4.0_DP, Y_STATE, WIND, R_SAIL, 0.0_DP)
    TARGET = 0.0_DP
    TARGET(INDEX) = RATE
    CALL CheckVec(OMEGA, TARGET, TOL, 'single Euler-rate case', NFAIL)
  END SUBROUTINE CheckSingleRate

  FUNCTION IndependentL2G(V, A) RESULT(OUT)
    REAL(DP), INTENT(IN) :: V(3), A(3)
    REAL(DP) :: OUT(3), T(3,3)
    T(1,1) = COS(A(2))*COS(A(3))
    T(2,1) = SIN(A(1))*SIN(A(2))*COS(A(3)) + COS(A(1))*SIN(A(3))
    T(3,1) =-COS(A(1))*SIN(A(2))*COS(A(3)) + SIN(A(1))*SIN(A(3))
    T(1,2) =-COS(A(2))*SIN(A(3))
    T(2,2) =-SIN(A(1))*SIN(A(2))*SIN(A(3)) + COS(A(1))*COS(A(3))
    T(3,2) = COS(A(1))*SIN(A(2))*SIN(A(3)) + SIN(A(1))*COS(A(3))
    T(1,3) = SIN(A(2))
    T(2,3) =-SIN(A(1))*COS(A(2))
    T(3,3) = COS(A(1))*COS(A(2))
    OUT = MATMUL(T, V)
  END FUNCTION IndependentL2G

  SUBROUTINE Check(CONDITION, LABEL, NFAIL)
    LOGICAL, INTENT(IN) :: CONDITION
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(INOUT) :: NFAIL
    IF (CONDITION) THEN
      WRITE(*,'(A,A)') '[PASS] ', LABEL
    ELSE
      WRITE(*,'(A,A)') '[FAIL] ', LABEL
      NFAIL = NFAIL + 1
    END IF
  END SUBROUTINE Check

  SUBROUTINE CheckVec(A, B, EPS, LABEL, NFAIL)
    REAL(DP), INTENT(IN) :: A(3), B(3), EPS
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(INOUT) :: NFAIL
    CALL Check(MAXVAL(ABS(A-B)) <= EPS, LABEL, NFAIL)
  END SUBROUTINE CheckVec

  SUBROUTINE CheckVec6(A, B, EPS, LABEL, NFAIL)
    REAL(DP), INTENT(IN) :: A(6), B(6), EPS
    CHARACTER(LEN=*), INTENT(IN) :: LABEL
    INTEGER, INTENT(INOUT) :: NFAIL
    CALL Check(MAXVAL(ABS(A-B)) <= EPS, LABEL, NFAIL)
  END SUBROUTINE CheckVec6

END PROGRAM test_sailplanning_adapter
