!===============================================================================
! SAILINTERP_MOD
!
! Piecewise-linear interpolation of the initialized sail aerodynamic database.
! The caller is responsible for mapping attack angles into the database range.
!===============================================================================
MODULE SAILINTERP_MOD
  USE SAILPARAM_MOD, ONLY: &
    DP, ANGLE_TOL_DEG, DATABASE_TOL, ALPHA_DB_MIN_DEG, &
    ALPHA_DB_MAX_DEG, SAIL_OK, SAIL_ERR_DATABASE_NOT_INITIALIZED, &
    SAIL_ERR_INVALID_INPUT, SAIL_ERR_INTERPOLATION
  USE SAILDATABASE_MOD, ONLY: &
    IsSailDatabaseInitialized, GetSailDatabaseSize, &
    GetSailDatabaseNode, GetSailDatabaseBounds
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: GetSailCoeff

CONTAINS

  SUBROUTINE GetSailCoeff(ALPHA_QUERY_DEG, CL_OUT, CD_OUT, IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: ALPHA_QUERY_DEG
    REAL(DP), INTENT(OUT) :: CL_OUT
    REAL(DP), INTENT(OUT) :: CD_OUT
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    INTEGER :: N_DATABASE
    INTEGER :: LOW
    INTEGER :: HIGH
    INTEGER :: MID
    INTEGER :: DATABASE_IERR
    REAL(DP) :: ALPHA_MIN
    REAL(DP) :: ALPHA_MAX
    REAL(DP) :: ALPHA_EFFECTIVE
    REAL(DP) :: ALPHA_1
    REAL(DP) :: ALPHA_2
    REAL(DP) :: CL_1
    REAL(DP) :: CL_2
    REAL(DP) :: CD_1
    REAL(DP) :: CD_2
    REAL(DP) :: DELTA_ALPHA
    REAL(DP) :: WEIGHT
    CHARACTER(LEN=2048) :: DATABASE_MESSAGE
    CHARACTER(LEN=2048) :: LOCAL_MESSAGE

    CL_OUT = 0.0_DP
    CD_OUT = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''

    IF (.NOT. IsSailDatabaseInitialized()) THEN
      IERR = SAIL_ERR_DATABASE_NOT_INITIALIZED
      MESSAGE = 'Sail database is not initialized.'
      RETURN
    END IF

    IF (.NOT. IEEE_IS_FINITE(ALPHA_QUERY_DEG)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Attack-angle query must be finite.'
      RETURN
    END IF

    CALL GetSailDatabaseBounds(ALPHA_MIN, ALPHA_MAX, DATABASE_IERR, &
      DATABASE_MESSAGE)
    IF (DATABASE_IERR /= SAIL_OK) THEN
      IERR = DATABASE_IERR
      WRITE(LOCAL_MESSAGE, '(A,A)') &
        'Database bounds query failed during interpolation: ', &
        TRIM(DATABASE_MESSAGE)
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    IF (.NOT. IEEE_IS_FINITE(ALPHA_MIN) .OR. &
        .NOT. IEEE_IS_FINITE(ALPHA_MAX)) THEN
      IERR = SAIL_ERR_INTERPOLATION
      MESSAGE = 'Database bounds are non-finite during interpolation.'
      RETURN
    END IF

    IF (ABS(ALPHA_MIN - ALPHA_DB_MIN_DEG) > ANGLE_TOL_DEG .OR. &
        ABS(ALPHA_MAX - ALPHA_DB_MAX_DEG) > ANGLE_TOL_DEG) THEN
      IERR = SAIL_ERR_INTERPOLATION
      MESSAGE = 'Database bounds do not match the required alpha range.'
      RETURN
    END IF

    IF (ALPHA_QUERY_DEG < ALPHA_MIN - ANGLE_TOL_DEG .OR. &
        ALPHA_QUERY_DEG > ALPHA_MAX + ANGLE_TOL_DEG) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      WRITE(LOCAL_MESSAGE, '(A,ES24.16,A,ES24.16,A,ES24.16,A)') &
        'Attack-angle query ', ALPHA_QUERY_DEG, &
        ' deg is outside database range [', ALPHA_MIN, ', ', &
        ALPHA_MAX, '] deg; extrapolation is not permitted.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    ! Clamp only the tolerance-sized roundoff admitted by the range check.
    ALPHA_EFFECTIVE = MAX(ALPHA_MIN, MIN(ALPHA_MAX, ALPHA_QUERY_DEG))

    N_DATABASE = GetSailDatabaseSize()
    IF (N_DATABASE < 2) THEN
      IERR = SAIL_ERR_INTERPOLATION
      WRITE(LOCAL_MESSAGE, '(A,I0,A)') &
        'At least two database nodes are required for interpolation; found ', &
        N_DATABASE, '.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    LOW = 1
    HIGH = N_DATABASE

    CALL ReadNode(LOW, ALPHA_1, CL_1, CD_1, IERR, MESSAGE)
    IF (IERR /= SAIL_OK) RETURN
    IF (ABS(ALPHA_EFFECTIVE - ALPHA_1) <= ANGLE_TOL_DEG) THEN
      CALL ReturnExactNode(ALPHA_QUERY_DEG, ALPHA_1, CL_1, CD_1, &
        CL_OUT, CD_OUT, IERR, MESSAGE)
      RETURN
    END IF

    CALL ReadNode(HIGH, ALPHA_2, CL_2, CD_2, IERR, MESSAGE)
    IF (IERR /= SAIL_OK) RETURN
    IF (ABS(ALPHA_EFFECTIVE - ALPHA_2) <= ANGLE_TOL_DEG) THEN
      CALL ReturnExactNode(ALPHA_QUERY_DEG, ALPHA_2, CL_2, CD_2, &
        CL_OUT, CD_OUT, IERR, MESSAGE)
      RETURN
    END IF

    ! Binary search uses actual node angles and therefore makes no assumption
    ! that the database spacing is uniform.
    DO WHILE (HIGH - LOW > 1)
      MID = LOW + (HIGH - LOW) / 2
      CALL ReadNode(MID, ALPHA_2, CL_2, CD_2, IERR, MESSAGE)
      IF (IERR /= SAIL_OK) RETURN

      IF (ABS(ALPHA_EFFECTIVE - ALPHA_2) <= ANGLE_TOL_DEG) THEN
        CALL ReturnExactNode(ALPHA_QUERY_DEG, ALPHA_2, CL_2, CD_2, &
          CL_OUT, CD_OUT, IERR, MESSAGE)
        RETURN
      ELSE IF (ALPHA_EFFECTIVE < ALPHA_2) THEN
        HIGH = MID
      ELSE
        LOW = MID
      END IF
    END DO

    ! Re-read both final neighbors: each accessor call is checked and no
    ! database arrays are copied or accessed directly.
    CALL ReadNode(LOW, ALPHA_1, CL_1, CD_1, IERR, MESSAGE)
    IF (IERR /= SAIL_OK) RETURN
    IF (ABS(ALPHA_EFFECTIVE - ALPHA_1) <= ANGLE_TOL_DEG) THEN
      CALL ReturnExactNode(ALPHA_QUERY_DEG, ALPHA_1, CL_1, CD_1, &
        CL_OUT, CD_OUT, IERR, MESSAGE)
      RETURN
    END IF

    CALL ReadNode(HIGH, ALPHA_2, CL_2, CD_2, IERR, MESSAGE)
    IF (IERR /= SAIL_OK) RETURN
    IF (ABS(ALPHA_EFFECTIVE - ALPHA_2) <= ANGLE_TOL_DEG) THEN
      CALL ReturnExactNode(ALPHA_QUERY_DEG, ALPHA_2, CL_2, CD_2, &
        CL_OUT, CD_OUT, IERR, MESSAGE)
      RETURN
    END IF

    IF (HIGH /= LOW + 1) THEN
      IERR = SAIL_ERR_INTERPOLATION
      MESSAGE = 'Binary search did not produce adjacent database nodes.'
      RETURN
    END IF

    DELTA_ALPHA = ALPHA_2 - ALPHA_1
    IF (.NOT. IEEE_IS_FINITE(DELTA_ALPHA) .OR. &
        DELTA_ALPHA <= DATABASE_TOL) THEN
      IERR = SAIL_ERR_INTERPOLATION
      WRITE(LOCAL_MESSAGE, '(A,ES24.16,A,I0,A,I0,A)') &
        'Invalid attack-angle interval width ', DELTA_ALPHA, &
        ' between database nodes ', LOW, ' and ', HIGH, '.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    WEIGHT = (ALPHA_EFFECTIVE - ALPHA_1) / DELTA_ALPHA
    IF (.NOT. IEEE_IS_FINITE(WEIGHT) .OR. &
        WEIGHT < -DATABASE_TOL .OR. WEIGHT > 1.0_DP + DATABASE_TOL) THEN
      IERR = SAIL_ERR_INTERPOLATION
      WRITE(LOCAL_MESSAGE, '(A,ES24.16,A)') &
        'Interpolation weight is outside the admissible range: ', &
        WEIGHT, '.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF
    WEIGHT = MAX(0.0_DP, MIN(1.0_DP, WEIGHT))

    CL_OUT = CL_1 + WEIGHT * (CL_2 - CL_1)
    CD_OUT = CD_1 + WEIGHT * (CD_2 - CD_1)
    IF (.NOT. IEEE_IS_FINITE(CL_OUT) .OR. &
        .NOT. IEEE_IS_FINITE(CD_OUT)) THEN
      CL_OUT = 0.0_DP
      CD_OUT = 0.0_DP
      IERR = SAIL_ERR_INTERPOLATION
      MESSAGE = 'Linear interpolation produced a non-finite coefficient.'
      RETURN
    END IF

    IERR = SAIL_OK
    WRITE(LOCAL_MESSAGE, '(A,ES16.8,A,ES16.8,A,ES16.8,A,ES16.8,A)') &
      'Linear interpolation completed for alpha = ', ALPHA_QUERY_DEG, &
      ' deg between ', ALPHA_1, ' and ', ALPHA_2, &
      ' deg; weight = ', WEIGHT, '.'
    MESSAGE = TRIM(LOCAL_MESSAGE)
  END SUBROUTINE GetSailCoeff


  SUBROUTINE ReadNode(INDEX_NODE, ALPHA_NODE, CL_NODE, CD_NODE, &
      IERR, MESSAGE)
    INTEGER, INTENT(IN) :: INDEX_NODE
    REAL(DP), INTENT(OUT) :: ALPHA_NODE
    REAL(DP), INTENT(OUT) :: CL_NODE
    REAL(DP), INTENT(OUT) :: CD_NODE
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    INTEGER :: DATABASE_IERR
    CHARACTER(LEN=2048) :: DATABASE_MESSAGE
    CHARACTER(LEN=2048) :: LOCAL_MESSAGE

    ALPHA_NODE = 0.0_DP
    CL_NODE = 0.0_DP
    CD_NODE = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''

    CALL GetSailDatabaseNode(INDEX_NODE, ALPHA_NODE, CL_NODE, CD_NODE, &
      DATABASE_IERR, DATABASE_MESSAGE)
    IF (DATABASE_IERR /= SAIL_OK) THEN
      IERR = DATABASE_IERR
      WRITE(LOCAL_MESSAGE, '(A,I0,A,A)') &
        'Database node query failed during interpolation at index ', &
        INDEX_NODE, ': ', TRIM(DATABASE_MESSAGE)
      MESSAGE = TRIM(LOCAL_MESSAGE)
    END IF
  END SUBROUTINE ReadNode


  SUBROUTINE ReturnExactNode(ALPHA_QUERY, ALPHA_NODE, CL_NODE, CD_NODE, &
      CL_OUT, CD_OUT, IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: ALPHA_QUERY
    REAL(DP), INTENT(IN) :: ALPHA_NODE
    REAL(DP), INTENT(IN) :: CL_NODE
    REAL(DP), INTENT(IN) :: CD_NODE
    REAL(DP), INTENT(OUT) :: CL_OUT
    REAL(DP), INTENT(OUT) :: CD_OUT
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    CHARACTER(LEN=2048) :: LOCAL_MESSAGE

    CL_OUT = CL_NODE
    CD_OUT = CD_NODE
    IERR = SAIL_OK
    WRITE(LOCAL_MESSAGE, '(A,ES16.8,A,ES16.8,A)') &
      'Exact database node returned for query alpha = ', ALPHA_QUERY, &
      ' deg at node alpha = ', ALPHA_NODE, ' deg.'
    MESSAGE = TRIM(LOCAL_MESSAGE)
  END SUBROUTINE ReturnExactNode

END MODULE SAILINTERP_MOD
