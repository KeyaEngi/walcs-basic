PROGRAM TEST_SAILMODULE
  USE SAILPARAM_MOD
  USE SAILMODULE_API_MOD
  USE SAILRELWIND_MOD, ONLY: ComputeSailRelativeWind
  USE SAILANGLE_MOD, ONLY: ComputeSailAngle
  USE SAILINTERP_MOD, ONLY: GetSailCoeff
  USE SAILFORCE_MOD, ONLY: ComputeSailForce
  USE, INTRINSIC :: IEEE_ARITHMETIC
  IMPLICIT NONE

  REAL(DP), PARAMETER :: SCALAR_TOL = 1.0E-10_DP
  REAL(DP), PARAMETER :: VECTOR_TOL = 1.0E-10_DP
  REAL(DP) :: V_WIND_BODY(3)
  REAL(DP) :: V_CG_BODY(3)
  REAL(DP) :: OMEGA_BODY(3)
  REAL(DP) :: R_SAIL_BODY(3)
  REAL(DP) :: DELTA_S_DEG
  REAL(DP) :: V_SAIL_POINT_BODY(3)
  REAL(DP) :: V_REL_BODY(3)
  REAL(DP) :: V_REL_H_BODY(3)
  REAL(DP) :: V_REL_H_MAG
  REAL(DP) :: DELTA_NORMALIZED_DEG
  REAL(DP) :: C_CHORD_BODY(3)
  REAL(DP) :: E_UPSTREAM_BODY(3)
  REAL(DP) :: ALPHA_RAW_DEG
  REAL(DP) :: ALPHA_DB_DEG
  REAL(DP) :: CL
  REAL(DP) :: CD
  REAL(DP) :: E_DRAG_BODY(3)
  REAL(DP) :: E_LIFT_BODY(3)
  REAL(DP) :: Q_DYNAMIC
  REAL(DP) :: FORCE_BODY(3)
  REAL(DP) :: MOMENT_BODY(3)
  REAL(DP) :: LOAD_6DOF(N_DOF)
  REAL(DP) :: EXPECTED_FORCE(3)
  REAL(DP) :: EXPECTED_MOMENT(3)
  REAL(DP) :: INPUT_WIND(3)
  REAL(DP) :: INPUT_CG(3)
  REAL(DP) :: INPUT_OMEGA(3)
  REAL(DP) :: INPUT_R(3)
  REAL(DP) :: INPUT_DELTA
  REAL(DP) :: SNAPSHOT_SCALARS(7)
  REAL(DP) :: SNAPSHOT_VECTORS(24)
  REAL(DP) :: D_V_SAIL_POINT(3)
  REAL(DP) :: D_V_REL(3)
  REAL(DP) :: D_V_REL_H(3)
  REAL(DP) :: D_V_REL_H_MAG
  REAL(DP) :: D_DELTA_NORMALIZED
  REAL(DP) :: D_C_CHORD(3)
  REAL(DP) :: D_E_UPSTREAM(3)
  REAL(DP) :: D_ALPHA_RAW
  REAL(DP) :: D_ALPHA_DB
  REAL(DP) :: D_CL
  REAL(DP) :: D_CD
  REAL(DP) :: D_E_DRAG(3)
  REAL(DP) :: D_E_LIFT(3)
  REAL(DP) :: D_Q_DYNAMIC
  REAL(DP) :: D_FORCE(3)
  REAL(DP) :: D_MOMENT(3)
  REAL(DP) :: D_LOAD(N_DOF)
  REAL(DP) :: NAN_VALUE
  REAL(DP) :: INF_VALUE
  REAL(DP) :: WEIGHT
  REAL(DP) :: EXPECTED_CL
  REAL(DP) :: EXPECTED_CD
  INTEGER :: IERR
  INTEGER :: DIRECT_IERR
  INTEGER :: N_PASS
  INTEGER :: N_FAIL
  INTEGER :: ITERATION
  CHARACTER(LEN=2048) :: MESSAGE
  CHARACTER(LEN=2048) :: DIRECT_MESSAGE
  LOGICAL :: REPEAT_OK

  N_PASS = 0
  N_FAIL = 0
  NAN_VALUE = IEEE_VALUE(0.0_DP, IEEE_QUIET_NAN)
  INF_VALUE = IEEE_VALUE(0.0_DP, IEEE_POSITIVE_INF)

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)

  CALL FinalizeSailModule(IERR, MESSAGE)
  CALL CallApi()
  CALL Check('uninitialized database call', &
    IERR == SAIL_ERR_DATABASE_NOT_INITIALIZED .AND. &
    AllOutputsZero() .AND. INDEX(MESSAGE, 'not initialized') > 0, &
    'database-not-initialized error and all-zero outputs')

  CALL InitializeSailModule(IERR=IERR, MESSAGE=MESSAGE)
  CALL Check('default database initialization', &
    IERR == SAIL_OK .AND. INDEX(MESSAGE, DEFAULT_DATABASE_FILE) > 0, &
    'SAIL_OK using DEFAULT_DATABASE_FILE')

  CALL InitializeSailModule(IERR=IERR, MESSAGE=MESSAGE)
  CALL Check('repeated initialization', &
    IERR == SAIL_OK .AND. INDEX(MESSAGE, 'reload skipped') > 0, &
    'SAIL_OK with reload skipped')

  CALL FinalizeSailModule(IERR, MESSAGE)
  CALL InitializeSailModule( &
    'file_that_does_not_exist.dat', IERR, MESSAGE)
  CALL Check('missing database file', &
    IERR == SAIL_ERR_FILE_NOT_FOUND .AND. &
    INDEX(MESSAGE, 'failed during database loading') > 0, &
    'SAIL_ERR_FILE_NOT_FOUND with API prefix')
  CALL CallApi()
  CALL Check('database remains clear after failed load', &
    IERR == SAIL_ERR_DATABASE_NOT_INITIALIZED .AND. AllOutputsZero(), &
    'database-not-initialized error and zero outputs')
  CALL InitializeSailModule(IERR=IERR, MESSAGE=MESSAGE)
  CALL Check('reinitialize after failed load', IERR == SAIL_OK, &
    'SAIL_OK')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  EXPECTED_FORCE = 61.25_DP * 15.0_DP * &
    (/ -0.29395224_DP, 0.17179593_DP, 0.0_DP /)
  CALL Check('STAR baseline complete chain', &
    IERR == SAIL_OK .AND. &
    VectorClose(V_SAIL_POINT_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
    VectorClose(V_REL_BODY, (/ 0.0_DP, 10.0_DP, 0.0_DP /)) .AND. &
    VectorClose(V_REL_H_BODY, (/ 0.0_DP, 10.0_DP, 0.0_DP /)) .AND. &
    ScalarClose(V_REL_H_MAG, 10.0_DP) .AND. &
    ScalarClose(DELTA_NORMALIZED_DEG, 0.0_DP) .AND. &
    VectorClose(C_CHORD_BODY, (/ 0.0_DP, -1.0_DP, 0.0_DP /)) .AND. &
    VectorClose(E_UPSTREAM_BODY, (/ 0.0_DP, -1.0_DP, 0.0_DP /)) .AND. &
    ScalarClose(ALPHA_DB_DEG, 0.0_DP) .AND. &
    ScalarClose(CL, -0.29395224_DP) .AND. &
    ScalarClose(CD, 0.17179593_DP) .AND. &
    VectorClose(E_DRAG_BODY, (/ 0.0_DP, 1.0_DP, 0.0_DP /)) .AND. &
    VectorClose(E_LIFT_BODY, (/ 1.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
    ScalarClose(Q_DYNAMIC, 61.25_DP) .AND. &
    VectorClose(FORCE_BODY, EXPECTED_FORCE) .AND. &
    VectorClose(MOMENT_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)), &
    'published STAR baseline values')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 2.0_DP, 3.0_DP, 4.0_DP /), 0.0_DP)
  CALL CallApi()
  EXPECTED_MOMENT(1) = -4.0_DP * EXPECTED_FORCE(2)
  EXPECTED_MOMENT(2) = 4.0_DP * EXPECTED_FORCE(1)
  EXPECTED_MOMENT(3) = &
    2.0_DP * EXPECTED_FORCE(2) - 3.0_DP * EXPECTED_FORCE(1)
  CALL Check('nonzero application-point moment', &
    IERR == SAIL_OK .AND. VectorClose(FORCE_BODY, EXPECTED_FORCE) .AND. &
    VectorClose(MOMENT_BODY, EXPECTED_MOMENT), &
    'MOMENT_BODY = R_SAIL_BODY cross FORCE_BODY')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 2.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('vessel translation effect', &
    IERR == SAIL_OK .AND. &
    VectorClose(V_REL_H_BODY, (/ 0.0_DP, 8.0_DP, 0.0_DP /)) .AND. &
    ScalarClose(V_REL_H_MAG, 8.0_DP) .AND. &
    ScalarClose(Q_DYNAMIC, 0.5_DP * RHO_AIR * 64.0_DP) .AND. &
    ScalarClose(ALPHA_DB_DEG, 0.0_DP), &
    'relative speed 8 m/s and q based on 8 m/s')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 1.0_DP /), &
    (/ 2.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('angular velocity and application-point effect', &
    IERR == SAIL_OK .AND. &
    VectorClose(V_SAIL_POINT_BODY, (/ 0.0_DP, 2.0_DP, 0.0_DP /)) .AND. &
    VectorClose(V_REL_H_BODY, (/ 0.0_DP, 8.0_DP, 0.0_DP /)), &
    'OMEGA cross R = (0,2,0) and relative wind = (0,8,0)')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 2.5_DP)
  CALL CallApi()
  WEIGHT = (2.5_DP - 0.0_DP) / (5.0_DP - 0.0_DP)
  EXPECTED_CL = -0.29395224_DP + &
    WEIGHT * (-0.59753960_DP + 0.29395224_DP)
  EXPECTED_CD = 0.17179593_DP + &
    WEIGHT * (0.21124380_DP - 0.17179593_DP)
  EXPECTED_FORCE = Q_DYNAMIC * SAIL_AREA * &
    (EXPECTED_CD * E_DRAG_BODY + EXPECTED_CL * E_LIFT_BODY)
  CALL Check('non-node interpolation complete chain', &
    IERR == SAIL_OK .AND. ScalarClose(ALPHA_DB_DEG, 2.5_DP) .AND. &
    ScalarClose(CL, EXPECTED_CL) .AND. ScalarClose(CD, EXPECTED_CD) .AND. &
    VectorClose(FORCE_BODY, EXPECTED_FORCE), &
    'independent linear interpolation between 0 and 5 degrees')

  CALL SetInputs( &
    (/ 8.0_DP, 12.0_DP, 1.0_DP /), &
    (/ 2.0_DP, -1.0_DP, 0.5_DP /), &
    (/ 0.1_DP, -0.2_DP, 0.3_DP /), &
    (/ 3.0_DP, 2.0_DP, 5.0_DP /), 25.0_DP)
  CALL CallApi()
  CALL ComputeSailRelativeWind( &
    V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, &
    D_V_SAIL_POINT, D_V_REL, D_V_REL_H, D_V_REL_H_MAG, &
    DIRECT_IERR, DIRECT_MESSAGE)
  IF (DIRECT_IERR == SAIL_OK) THEN
    CALL ComputeSailAngle( &
      D_V_REL_H, D_V_REL_H_MAG, DELTA_S_DEG, D_DELTA_NORMALIZED, &
      D_C_CHORD, D_E_UPSTREAM, D_ALPHA_RAW, D_ALPHA_DB, &
      DIRECT_IERR, DIRECT_MESSAGE)
  END IF
  IF (DIRECT_IERR == SAIL_OK) THEN
    CALL GetSailCoeff( &
      D_ALPHA_DB, D_CL, D_CD, DIRECT_IERR, DIRECT_MESSAGE)
  END IF
  IF (DIRECT_IERR == SAIL_OK) THEN
    CALL ComputeSailForce( &
      D_V_REL_H, D_V_REL_H_MAG, D_CL, D_CD, R_SAIL_BODY, &
      D_E_DRAG, D_E_LIFT, D_Q_DYNAMIC, D_FORCE, D_MOMENT, D_LOAD, &
      DIRECT_IERR, DIRECT_MESSAGE)
  END IF
  CALL Check('API and direct-chain consistency', &
    IERR == SAIL_OK .AND. DIRECT_IERR == SAIL_OK .AND. &
    VectorClose(V_SAIL_POINT_BODY, D_V_SAIL_POINT) .AND. &
    VectorClose(V_REL_BODY, D_V_REL) .AND. &
    VectorClose(V_REL_H_BODY, D_V_REL_H) .AND. &
    ScalarClose(V_REL_H_MAG, D_V_REL_H_MAG) .AND. &
    ScalarClose(DELTA_NORMALIZED_DEG, D_DELTA_NORMALIZED) .AND. &
    VectorClose(C_CHORD_BODY, D_C_CHORD) .AND. &
    VectorClose(E_UPSTREAM_BODY, D_E_UPSTREAM) .AND. &
    ScalarClose(ALPHA_RAW_DEG, D_ALPHA_RAW) .AND. &
    ScalarClose(ALPHA_DB_DEG, D_ALPHA_DB) .AND. &
    ScalarClose(CL, D_CL) .AND. ScalarClose(CD, D_CD) .AND. &
    VectorClose(E_DRAG_BODY, D_E_DRAG) .AND. &
    VectorClose(E_LIFT_BODY, D_E_LIFT) .AND. &
    ScalarClose(Q_DYNAMIC, D_Q_DYNAMIC) .AND. &
    VectorClose(FORCE_BODY, D_FORCE) .AND. &
    VectorClose(MOMENT_BODY, D_MOMENT) .AND. &
    VectorCloseN(LOAD_6DOF, D_LOAD), &
    'all API outputs equal sequential public-module calls')

  CALL SetInputs( &
    (/ 4.0_DP, -2.0_DP, 1.0_DP /), &
    (/ 4.0_DP, -2.0_DP, 1.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('low horizontal relative wind', &
    IERR == SAIL_ERR_LOW_WIND_SPEED .AND. AllLoadsZero() .AND. &
    INDEX(MESSAGE, 'Nonfatal') > 0, &
    'recognizable low-wind error and zero downstream loads')

  CALL SetInputs( &
    (/ 0.0_DP, 0.0_DP, 5.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('vertical-only relative wind', &
    IERR == SAIL_ERR_LOW_WIND_SPEED .AND. &
    VectorClose(V_REL_BODY, (/ 0.0_DP, 0.0_DP, 5.0_DP /)) .AND. &
    VectorClose(V_REL_H_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
    ScalarClose(V_REL_H_MAG, 0.0_DP) .AND. AllLoadsZero(), &
    'vertical relative wind retained and aerodynamic loads zero')

  CALL SetInputs( &
    (/ NAN_VALUE, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('NaN wind error propagation', &
    IERR == SAIL_ERR_INVALID_INPUT .AND. &
    INDEX(MESSAGE, 'Relative-wind stage failed') > 0 .AND. &
    AllOutputsZero(), 'relative-wind invalid-input error and zero outputs')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, INF_VALUE /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('infinite angular-velocity error propagation', &
    IERR == SAIL_ERR_INVALID_INPUT .AND. &
    INDEX(MESSAGE, 'Relative-wind stage failed') > 0 .AND. &
    AllOutputsZero(), 'relative-wind invalid-input error and zero outputs')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), NAN_VALUE)
  CALL CallApi()
  CALL Check('NaN sail-angle error propagation', &
    IERR == SAIL_ERR_INVALID_INPUT .AND. &
    INDEX(MESSAGE, 'Sail-angle stage failed') > 0 .AND. AllLoadsZero(), &
    'sail-angle invalid-input error and zero downstream loads')

  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ INF_VALUE, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('infinite sail-position error propagation', &
    IERR == SAIL_ERR_INVALID_INPUT .AND. &
    INDEX(MESSAGE, 'Relative-wind stage failed') > 0 .AND. &
    AllOutputsZero(), 'relative-wind invalid-input error and zero outputs')

  CALL SetInputs( &
    (/ 8.0_DP, 12.0_DP, 1.0_DP /), &
    (/ 2.0_DP, -1.0_DP, 0.5_DP /), &
    (/ 0.1_DP, -0.2_DP, 0.3_DP /), &
    (/ 3.0_DP, 2.0_DP, 5.0_DP /), 25.0_DP)
  INPUT_WIND = V_WIND_BODY
  INPUT_CG = V_CG_BODY
  INPUT_OMEGA = OMEGA_BODY
  INPUT_R = R_SAIL_BODY
  INPUT_DELTA = DELTA_S_DEG
  CALL CallApi()
  SNAPSHOT_SCALARS = (/ V_REL_H_MAG, DELTA_NORMALIZED_DEG, &
    ALPHA_RAW_DEG, ALPHA_DB_DEG, CL, CD, Q_DYNAMIC /)
  SNAPSHOT_VECTORS(1:3) = V_SAIL_POINT_BODY
  SNAPSHOT_VECTORS(4:6) = V_REL_BODY
  SNAPSHOT_VECTORS(7:9) = V_REL_H_BODY
  SNAPSHOT_VECTORS(10:12) = FORCE_BODY
  SNAPSHOT_VECTORS(13:15) = MOMENT_BODY
  SNAPSHOT_VECTORS(16:21) = LOAD_6DOF
  SNAPSHOT_VECTORS(22:24) = C_CHORD_BODY
  REPEAT_OK = IERR == SAIL_OK
  DO ITERATION = 1, 1000
    CALL CallApi()
    REPEAT_OK = REPEAT_OK .AND. IERR == SAIL_OK
    REPEAT_OK = REPEAT_OK .AND. VectorCloseN(SNAPSHOT_SCALARS, &
      (/ V_REL_H_MAG, DELTA_NORMALIZED_DEG, ALPHA_RAW_DEG, &
         ALPHA_DB_DEG, CL, CD, Q_DYNAMIC /))
    REPEAT_OK = REPEAT_OK .AND. &
      VectorClose(SNAPSHOT_VECTORS(1:3), V_SAIL_POINT_BODY)
    REPEAT_OK = REPEAT_OK .AND. &
      VectorClose(SNAPSHOT_VECTORS(4:6), V_REL_BODY)
    REPEAT_OK = REPEAT_OK .AND. &
      VectorClose(SNAPSHOT_VECTORS(7:9), V_REL_H_BODY)
    REPEAT_OK = REPEAT_OK .AND. &
      VectorClose(SNAPSHOT_VECTORS(10:12), FORCE_BODY)
    REPEAT_OK = REPEAT_OK .AND. &
      VectorClose(SNAPSHOT_VECTORS(13:15), MOMENT_BODY)
    REPEAT_OK = REPEAT_OK .AND. &
      VectorCloseN(SNAPSHOT_VECTORS(16:21), LOAD_6DOF)
    REPEAT_OK = REPEAT_OK .AND. &
      VectorClose(SNAPSHOT_VECTORS(22:24), C_CHORD_BODY)
  END DO
  CALL Check('1000 repeated complete-chain calls', REPEAT_OK, &
    'all 1000 results equal the first result without reinitialization')

  CALL Check('input arguments remain unchanged', &
    VectorClose(V_WIND_BODY, INPUT_WIND) .AND. &
    VectorClose(V_CG_BODY, INPUT_CG) .AND. &
    VectorClose(OMEGA_BODY, INPUT_OMEGA) .AND. &
    VectorClose(R_SAIL_BODY, INPUT_R) .AND. &
    ScalarClose(DELTA_S_DEG, INPUT_DELTA), &
    'all five input arguments unchanged')

  CALL FinalizeSailModule(IERR, MESSAGE)
  CALL CallApi()
  CALL Check('calculation after finalization', &
    IERR == SAIL_ERR_DATABASE_NOT_INITIALIZED .AND. AllOutputsZero(), &
    'database-not-initialized error and zero outputs')

  CALL InitializeSailModule(IERR=IERR, MESSAGE=MESSAGE)
  CALL SetInputs( &
    (/ 0.0_DP, 10.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), &
    (/ 0.0_DP, 0.0_DP, 0.0_DP /), 0.0_DP)
  CALL CallApi()
  CALL Check('reinitialize lifecycle and repeat baseline', &
    IERR == SAIL_OK .AND. ScalarClose(CL, -0.29395224_DP) .AND. &
    ScalarClose(CD, 0.17179593_DP), &
    'successful reinitialization and original baseline coefficients')

  CALL FinalizeSailModule(IERR, MESSAGE)
  CALL Check('final cleanup', IERR == SAIL_OK, 'SAIL_OK')

  WRITE(*, '(A)') '========================================'
  WRITE(*, '(A)') 'SAILMODULE END-TO-END TEST SUMMARY'
  WRITE(*, '(A,I0)') 'PASS: ', N_PASS
  WRITE(*, '(A,I0)') 'FAIL: ', N_FAIL
  WRITE(*, '(A)') '========================================'
  IF (N_FAIL > 0) STOP 1

CONTAINS

  SUBROUTINE SetInputs(WIND, CG, OMEGA, R_SAIL, DELTA)
    REAL(DP), INTENT(IN) :: WIND(:)
    REAL(DP), INTENT(IN) :: CG(:)
    REAL(DP), INTENT(IN) :: OMEGA(:)
    REAL(DP), INTENT(IN) :: R_SAIL(:)
    REAL(DP), INTENT(IN) :: DELTA

    V_WIND_BODY = WIND
    V_CG_BODY = CG
    OMEGA_BODY = OMEGA
    R_SAIL_BODY = R_SAIL
    DELTA_S_DEG = DELTA
  END SUBROUTINE SetInputs


  SUBROUTINE CallApi()
    CALL ComputeSailModuleLoads( &
      V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, DELTA_S_DEG, &
      V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, V_REL_H_MAG, &
      DELTA_NORMALIZED_DEG, C_CHORD_BODY, E_UPSTREAM_BODY, &
      ALPHA_RAW_DEG, ALPHA_DB_DEG, CL, CD, E_DRAG_BODY, &
      E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, LOAD_6DOF, &
      IERR, MESSAGE)
  END SUBROUTINE CallApi


  SUBROUTINE Check(NAME, CONDITION, EXPECTED)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    LOGICAL, INTENT(IN) :: CONDITION
    CHARACTER(LEN=*), INTENT(IN) :: EXPECTED

    IF (CONDITION) THEN
      N_PASS = N_PASS + 1
      WRITE(*, '(A,A)') '[PASS] ', TRIM(NAME)
    ELSE
      N_FAIL = N_FAIL + 1
      WRITE(*, '(A,A)') '[FAIL] ', TRIM(NAME)
      WRITE(*, '(A,I0)') '  IERR: ', IERR
      WRITE(*, '(A,A)') '  MESSAGE: ', TRIM(MESSAGE)
      WRITE(*, '(A,3(ES16.8,1X))') '  V_REL_BODY: ', V_REL_BODY
      WRITE(*, '(A,3(ES16.8,1X))') '  V_REL_H_BODY: ', V_REL_H_BODY
      WRITE(*, '(A,ES16.8)') '  V_REL_H_MAG: ', V_REL_H_MAG
      WRITE(*, '(A,ES16.8)') '  ALPHA_RAW_DEG: ', ALPHA_RAW_DEG
      WRITE(*, '(A,ES16.8)') '  ALPHA_DB_DEG: ', ALPHA_DB_DEG
      WRITE(*, '(A,ES16.8)') '  CL: ', CL
      WRITE(*, '(A,ES16.8)') '  CD: ', CD
      WRITE(*, '(A,3(ES16.8,1X))') '  FORCE_BODY: ', FORCE_BODY
      WRITE(*, '(A,3(ES16.8,1X))') '  MOMENT_BODY: ', MOMENT_BODY
      WRITE(*, '(A,6(ES16.8,1X))') '  LOAD_6DOF: ', LOAD_6DOF
      WRITE(*, '(A,A)') '  EXPECTED: ', TRIM(EXPECTED)
    END IF
  END SUBROUTINE Check


  LOGICAL FUNCTION ScalarClose(A, B)
    REAL(DP), INTENT(IN) :: A
    REAL(DP), INTENT(IN) :: B

    ScalarClose = ABS(A - B) <= &
      SCALAR_TOL * MAX(1.0_DP, ABS(A), ABS(B))
  END FUNCTION ScalarClose


  LOGICAL FUNCTION VectorClose(A, B)
    REAL(DP), INTENT(IN) :: A(:)
    REAL(DP), INTENT(IN) :: B(:)

    VectorClose = ALL(ABS(A - B) <= &
      VECTOR_TOL * MAX(1.0_DP, MAXVAL(ABS(A)), MAXVAL(ABS(B))))
  END FUNCTION VectorClose


  LOGICAL FUNCTION VectorCloseN(A, B)
    REAL(DP), INTENT(IN) :: A(:)
    REAL(DP), INTENT(IN) :: B(:)

    VectorCloseN = SIZE(A) == SIZE(B)
    IF (VectorCloseN) THEN
      VectorCloseN = ALL(ABS(A - B) <= &
        VECTOR_TOL * MAX(1.0_DP, MAXVAL(ABS(A)), MAXVAL(ABS(B))))
    END IF
  END FUNCTION VectorCloseN


  LOGICAL FUNCTION AllLoadsZero()
    AllLoadsZero = &
      ScalarClose(DELTA_NORMALIZED_DEG, 0.0_DP) .AND. &
      VectorClose(C_CHORD_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      VectorClose(E_UPSTREAM_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      ScalarClose(ALPHA_RAW_DEG, 0.0_DP) .AND. &
      ScalarClose(ALPHA_DB_DEG, 0.0_DP) .AND. &
      ScalarClose(CL, 0.0_DP) .AND. ScalarClose(CD, 0.0_DP) .AND. &
      VectorClose(E_DRAG_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      VectorClose(E_LIFT_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      ScalarClose(Q_DYNAMIC, 0.0_DP) .AND. &
      VectorClose(FORCE_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      VectorClose(MOMENT_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      VectorCloseN(LOAD_6DOF, (/ 0.0_DP, 0.0_DP, 0.0_DP, &
        0.0_DP, 0.0_DP, 0.0_DP /))
  END FUNCTION AllLoadsZero


  LOGICAL FUNCTION AllOutputsZero()
    AllOutputsZero = &
      VectorClose(V_SAIL_POINT_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      VectorClose(V_REL_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      VectorClose(V_REL_H_BODY, (/ 0.0_DP, 0.0_DP, 0.0_DP /)) .AND. &
      ScalarClose(V_REL_H_MAG, 0.0_DP) .AND. AllLoadsZero()
  END FUNCTION AllOutputsZero

END PROGRAM TEST_SAILMODULE
