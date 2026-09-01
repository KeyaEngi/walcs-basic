!===============================================================================
! Public orchestration API for the complete single-sail calculation chain.
!===============================================================================
MODULE SAILMODULE_API_MOD
  USE SAILPARAM_MOD, ONLY: &
    DP, N_DOF, DEFAULT_DATABASE_FILE, SAIL_OK, &
    SAIL_ERR_DATABASE_NOT_INITIALIZED, SAIL_ERR_LOW_WIND_SPEED
  USE SAILDATABASE_MOD, ONLY: &
    ReadSailDatabase, ClearSailDatabase, IsSailDatabaseInitialized
  USE SAILRELWIND_MOD, ONLY: ComputeSailRelativeWind
  USE SAILANGLE_MOD, ONLY: ComputeSailAngle
  USE SAILINTERP_MOD, ONLY: GetSailCoeff
  USE SAILFORCE_MOD, ONLY: ComputeSailForce
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSailModule
  PUBLIC :: ComputeSailModuleLoads
  PUBLIC :: FinalizeSailModule

CONTAINS

  SUBROUTINE InitializeSailModule(DATABASE_FILE, IERR, MESSAGE)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: DATABASE_FILE
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    CHARACTER(LEN=2048) :: LOCAL_MESSAGE
    CHARACTER(LEN=2048) :: API_MESSAGE
    CHARACTER(LEN=1024) :: ACTUAL_FILE

    IERR = SAIL_OK
    MESSAGE = ''
    LOCAL_MESSAGE = ''
    API_MESSAGE = ''

    IF (PRESENT(DATABASE_FILE)) THEN
      ACTUAL_FILE = TRIM(DATABASE_FILE)
    ELSE
      ACTUAL_FILE = DEFAULT_DATABASE_FILE
    END IF

    CALL ReadSailDatabase(TRIM(ACTUAL_FILE), IERR, LOCAL_MESSAGE)
    IF (IERR /= SAIL_OK) THEN
      WRITE(API_MESSAGE, '(A,A)') &
        'InitializeSailModule failed during database loading: ', &
        TRIM(LOCAL_MESSAGE)
      MESSAGE = TRIM(API_MESSAGE)
      RETURN
    END IF

    WRITE(API_MESSAGE, '(A,A,A,A)') &
      'Sail module initialized successfully: ', TRIM(ACTUAL_FILE), &
      '. ', TRIM(LOCAL_MESSAGE)
    MESSAGE = TRIM(API_MESSAGE)
  END SUBROUTINE InitializeSailModule


  SUBROUTINE FinalizeSailModule(IERR, MESSAGE)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    CALL ClearSailDatabase()
    IERR = SAIL_OK
    MESSAGE = 'Sail module finalized successfully; database storage cleared.'
  END SUBROUTINE FinalizeSailModule


  SUBROUTINE ComputeSailModuleLoads( &
      V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, DELTA_S_DEG, &
      V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, V_REL_H_MAG, &
      DELTA_NORMALIZED_DEG, C_CHORD_BODY, E_UPSTREAM_BODY, &
      ALPHA_RAW_DEG, ALPHA_DB_DEG, CL, CD, E_DRAG_BODY, &
      E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, LOAD_6DOF, &
      IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: V_WIND_BODY(3)
    REAL(DP), INTENT(IN) :: V_CG_BODY(3)
    REAL(DP), INTENT(IN) :: OMEGA_BODY(3)
    REAL(DP), INTENT(IN) :: R_SAIL_BODY(3)
    REAL(DP), INTENT(IN) :: DELTA_S_DEG
    REAL(DP), INTENT(OUT) :: V_SAIL_POINT_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_H_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_H_MAG
    REAL(DP), INTENT(OUT) :: DELTA_NORMALIZED_DEG
    REAL(DP), INTENT(OUT) :: C_CHORD_BODY(3)
    REAL(DP), INTENT(OUT) :: E_UPSTREAM_BODY(3)
    REAL(DP), INTENT(OUT) :: ALPHA_RAW_DEG
    REAL(DP), INTENT(OUT) :: ALPHA_DB_DEG
    REAL(DP), INTENT(OUT) :: CL
    REAL(DP), INTENT(OUT) :: CD
    REAL(DP), INTENT(OUT) :: E_DRAG_BODY(3)
    REAL(DP), INTENT(OUT) :: E_LIFT_BODY(3)
    REAL(DP), INTENT(OUT) :: Q_DYNAMIC
    REAL(DP), INTENT(OUT) :: FORCE_BODY(3)
    REAL(DP), INTENT(OUT) :: MOMENT_BODY(3)
    REAL(DP), INTENT(OUT) :: LOAD_6DOF(N_DOF)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    INTEGER :: LOCAL_IERR
    CHARACTER(LEN=2048) :: LOCAL_MESSAGE
    CHARACTER(LEN=2048) :: API_MESSAGE

    V_SAIL_POINT_BODY = 0.0_DP
    V_REL_BODY = 0.0_DP
    V_REL_H_BODY = 0.0_DP
    V_REL_H_MAG = 0.0_DP
    DELTA_NORMALIZED_DEG = 0.0_DP
    C_CHORD_BODY = 0.0_DP
    E_UPSTREAM_BODY = 0.0_DP
    ALPHA_RAW_DEG = 0.0_DP
    ALPHA_DB_DEG = 0.0_DP
    CL = 0.0_DP
    CD = 0.0_DP
    E_DRAG_BODY = 0.0_DP
    E_LIFT_BODY = 0.0_DP
    Q_DYNAMIC = 0.0_DP
    FORCE_BODY = 0.0_DP
    MOMENT_BODY = 0.0_DP
    LOAD_6DOF = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''
    LOCAL_MESSAGE = ''
    API_MESSAGE = ''

    IF (.NOT. IsSailDatabaseInitialized()) THEN
      IERR = SAIL_ERR_DATABASE_NOT_INITIALIZED
      MESSAGE = &
        'ComputeSailModuleLoads failed: database is not initialized.'
      RETURN
    END IF

    CALL ComputeSailRelativeWind( &
      V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, &
      V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, V_REL_H_MAG, &
      LOCAL_IERR, LOCAL_MESSAGE)
    IF (LOCAL_IERR == SAIL_ERR_LOW_WIND_SPEED) THEN
      IERR = LOCAL_IERR
      WRITE(API_MESSAGE, '(A,A)') &
        'Low horizontal relative wind: angle of attack is undefined; ', &
        'sail loads were set to zero. Nonfatal host-recognizable status.'
      MESSAGE = TRIM(API_MESSAGE)
      RETURN
    ELSE IF (LOCAL_IERR /= SAIL_OK) THEN
      IERR = LOCAL_IERR
      WRITE(API_MESSAGE, '(A,A)') &
        'Relative-wind stage failed: ', TRIM(LOCAL_MESSAGE)
      MESSAGE = TRIM(API_MESSAGE)
      RETURN
    END IF

    CALL ComputeSailAngle( &
      V_REL_H_BODY, V_REL_H_MAG, DELTA_S_DEG, &
      DELTA_NORMALIZED_DEG, C_CHORD_BODY, E_UPSTREAM_BODY, &
      ALPHA_RAW_DEG, ALPHA_DB_DEG, LOCAL_IERR, LOCAL_MESSAGE)
    IF (LOCAL_IERR /= SAIL_OK) THEN
      IERR = LOCAL_IERR
      WRITE(API_MESSAGE, '(A,A)') &
        'Sail-angle stage failed: ', TRIM(LOCAL_MESSAGE)
      MESSAGE = TRIM(API_MESSAGE)
      RETURN
    END IF

    CALL GetSailCoeff( &
      ALPHA_DB_DEG, CL, CD, LOCAL_IERR, LOCAL_MESSAGE)
    IF (LOCAL_IERR /= SAIL_OK) THEN
      CL = 0.0_DP
      CD = 0.0_DP
      IERR = LOCAL_IERR
      WRITE(API_MESSAGE, '(A,A)') &
        'Interpolation stage failed: ', TRIM(LOCAL_MESSAGE)
      MESSAGE = TRIM(API_MESSAGE)
      RETURN
    END IF

    CALL ComputeSailForce( &
      V_REL_H_BODY, V_REL_H_MAG, CL, CD, R_SAIL_BODY, &
      E_DRAG_BODY, E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, &
      LOAD_6DOF, LOCAL_IERR, LOCAL_MESSAGE)
    IF (LOCAL_IERR /= SAIL_OK) THEN
      IERR = LOCAL_IERR
      WRITE(API_MESSAGE, '(A,A)') &
        'Force stage failed: ', TRIM(LOCAL_MESSAGE)
      MESSAGE = TRIM(API_MESSAGE)
      RETURN
    END IF

    IERR = SAIL_OK
    WRITE(API_MESSAGE, &
      '(A,ES13.5,A,ES13.5,A,ES13.5,A,ES13.5,A,ES13.5,A,ES13.5,' // &
      'A,3(ES13.5,1X),A,3(ES13.5,1X),A)') &
      'Sail module calculation succeeded; speed=', V_REL_H_MAG, &
      ', delta=', DELTA_NORMALIZED_DEG, ', alpha_raw=', ALPHA_RAW_DEG, &
      ', alpha_db=', ALPHA_DB_DEG, ', CL=', CL, ', CD=', CD, &
      ', force=(', FORCE_BODY, '), moment=(', MOMENT_BODY, ').'
    MESSAGE = TRIM(API_MESSAGE)
  END SUBROUTINE ComputeSailModuleLoads

END MODULE SAILMODULE_API_MOD
