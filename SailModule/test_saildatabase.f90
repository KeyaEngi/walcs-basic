PROGRAM TEST_SAILDATABASE
  USE SAILPARAM_MOD
  USE SAILDATABASE_MOD
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE

  INTEGER :: I
  INTEGER :: IERR
  INTEGER :: N_PASS
  INTEGER :: N_FAIL
  REAL(DP) :: ALPHA
  REAL(DP) :: CL
  REAL(DP) :: CD
  REAL(DP) :: ALPHA_MIN
  REAL(DP) :: ALPHA_MAX
  REAL(DP) :: PREVIOUS_ALPHA
  LOGICAL :: OK
  CHARACTER(LEN=2048) :: MESSAGE
  CHARACTER(LEN=2048) :: DETAIL

  N_PASS = 0
  N_FAIL = 0

  CALL ClearSailDatabase()
  OK = .NOT. IsSailDatabaseInitialized() .AND. &
    GetSailDatabaseSize() == 0
  WRITE(DETAIL, '(A,L1,A,I0)') 'initialized=', &
    IsSailDatabaseInitialized(), ', size=', GetSailDatabaseSize()
  CALL ReportResult(1, 'Initial state', OK, DETAIL)

  CALL GetSailDatabaseNode(1, ALPHA, CL, CD, IERR, MESSAGE)
  OK = IERR == SAIL_ERR_DATABASE_NOT_INITIALIZED .AND. &
    ALPHA == 0.0_DP .AND. CL == 0.0_DP .AND. CD == 0.0_DP
  WRITE(DETAIL, '(A,I0,A,3(ES14.6,1X))') 'IERR=', IERR, &
    ', outputs=', ALPHA, CL, CD
  CALL ReportResult(2, 'Uninitialized node query', OK, DETAIL)

  CALL GetSailDatabaseBounds(ALPHA_MIN, ALPHA_MAX, IERR, MESSAGE)
  OK = IERR == SAIL_ERR_DATABASE_NOT_INITIALIZED .AND. &
    ALPHA_MIN == 0.0_DP .AND. ALPHA_MAX == 0.0_DP
  WRITE(DETAIL, '(A,I0,A,2(ES14.6,1X))') 'IERR=', IERR, &
    ', bounds=', ALPHA_MIN, ALPHA_MAX
  CALL ReportResult(3, 'Uninitialized bounds query', OK, DETAIL)

  CALL ReadSailDatabase('sail_database.dat', IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. IsSailDatabaseInitialized() .AND. &
    GetSailDatabaseSize() == 33
  WRITE(DETAIL, '(A,I0,A,L1,A,I0,A,A)') 'IERR=', IERR, &
    ', initialized=', IsSailDatabaseInitialized(), ', size=', &
    GetSailDatabaseSize(), ', message=', TRIM(MESSAGE)
  CALL ReportResult(4, 'Normal database read', OK, DETAIL)

  CALL GetSailDatabaseBounds(ALPHA_MIN, ALPHA_MAX, IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. &
    ABS(ALPHA_MIN - 0.0_DP) <= ANGLE_TOL_DEG .AND. &
    ABS(ALPHA_MAX - 180.0_DP) <= ANGLE_TOL_DEG
  WRITE(DETAIL, '(A,I0,A,2(ES14.6,1X))') 'IERR=', IERR, &
    ', bounds=', ALPHA_MIN, ALPHA_MAX
  CALL ReportResult(5, 'Database bounds', OK, DETAIL)

  CALL GetSailDatabaseNode(1, ALPHA, CL, CD, IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. ABS(ALPHA) <= 1.0E-10_DP .AND. &
    ABS(CL - (-0.29395224_DP)) <= 1.0E-10_DP .AND. &
    ABS(CD - 0.17179593_DP) <= 1.0E-10_DP
  WRITE(DETAIL, '(A,I0,A,3(ES14.6,1X))') 'IERR=', IERR, &
    ', node=', ALPHA, CL, CD
  CALL ReportResult(6, 'First database node', OK, DETAIL)

  CALL GetSailDatabaseNode(33, ALPHA, CL, CD, IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. &
    ABS(ALPHA - 180.0_DP) <= 1.0E-10_DP .AND. &
    ABS(CL - 0.26555131_DP) <= 1.0E-10_DP .AND. &
    ABS(CD - 0.17291675_DP) <= 1.0E-10_DP
  WRITE(DETAIL, '(A,I0,A,3(ES14.6,1X))') 'IERR=', IERR, &
    ', node=', ALPHA, CL, CD
  CALL ReportResult(7, 'Last database node', OK, DETAIL)

  OK = .TRUE.
  PREVIOUS_ALPHA = -HUGE(0.0_DP)
  DO I = 1, 33
    CALL GetSailDatabaseNode(I, ALPHA, CL, CD, IERR, MESSAGE)
    IF (IERR /= SAIL_OK .OR. ALPHA <= PREVIOUS_ALPHA) OK = .FALSE.
    PREVIOUS_ALPHA = ALPHA
  END DO
  WRITE(DETAIL, '(A,L1,A,ES14.6)') 'strictly increasing=', OK, &
    ', final alpha=', PREVIOUS_ALPHA
  CALL ReportResult(8, 'Strictly increasing attack angles', OK, DETAIL)

  OK = .TRUE.
  DO I = 1, 33
    CALL GetSailDatabaseNode(I, ALPHA, CL, CD, IERR, MESSAGE)
    IF (IERR /= SAIL_OK .OR. .NOT. IEEE_IS_FINITE(ALPHA) .OR. &
        .NOT. IEEE_IS_FINITE(CL) .OR. .NOT. IEEE_IS_FINITE(CD) .OR. &
        CD < -DATABASE_TOL) OK = .FALSE.
  END DO
  WRITE(DETAIL, '(A,L1,A,I0)') 'all finite and CD valid=', OK, &
    ', nodes checked=', 33
  CALL ReportResult(9, 'Finite values and drag coefficients', OK, DETAIL)

  OK = .TRUE.
  CALL GetSailDatabaseNode(0, ALPHA, CL, CD, IERR, MESSAGE)
  IF (IERR /= SAIL_ERR_INVALID_INPUT .OR. ALPHA /= 0.0_DP .OR. &
      CL /= 0.0_DP .OR. CD /= 0.0_DP) OK = .FALSE.
  CALL GetSailDatabaseNode(34, ALPHA, CL, CD, IERR, MESSAGE)
  IF (IERR /= SAIL_ERR_INVALID_INPUT .OR. ALPHA /= 0.0_DP .OR. &
      CL /= 0.0_DP .OR. CD /= 0.0_DP) OK = .FALSE.
  WRITE(DETAIL, '(A,L1,A,I0,A,3(ES14.6,1X))') &
    'both invalid indices rejected=', OK, ', last IERR=', IERR, &
    ', last outputs=', ALPHA, CL, CD
  CALL ReportResult(10, 'Invalid node indices', OK, DETAIL)

  CALL ReadSailDatabase('sail_database.dat', IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. IsSailDatabaseInitialized() .AND. &
    GetSailDatabaseSize() == 33 .AND. &
    (INDEX(MESSAGE, 'already initialized') > 0 .OR. &
     INDEX(MESSAGE, 'reload skipped') > 0)
  WRITE(DETAIL, '(A,I0,A,I0,A,A)') 'IERR=', IERR, ', size=', &
    GetSailDatabaseSize(), ', message=', TRIM(MESSAGE)
  CALL ReportResult(11, 'Same-file reload', OK, DETAIL)

  CALL ClearSailDatabase()
  CALL ReadSailDatabase(IERR=IERR, MESSAGE=MESSAGE)
  OK = IERR == SAIL_OK .AND. IsSailDatabaseInitialized() .AND. &
    GetSailDatabaseSize() == 33
  WRITE(DETAIL, '(A,I0,A,I0,A,A)') 'IERR=', IERR, ', size=', &
    GetSailDatabaseSize(), ', message=', TRIM(MESSAGE)
  CALL ReportResult(12, 'Default database file', OK, DETAIL)

  CALL ClearSailDatabase()
  CALL ReadSailDatabase('file_that_does_not_exist.dat', IERR, MESSAGE)
  OK = IERR == SAIL_ERR_FILE_NOT_FOUND .AND. &
    .NOT. IsSailDatabaseInitialized() .AND. GetSailDatabaseSize() == 0
  WRITE(DETAIL, '(A,I0,A,L1,A,I0,A,A)') 'IERR=', IERR, &
    ', initialized=', IsSailDatabaseInitialized(), ', size=', &
    GetSailDatabaseSize(), ', message=', TRIM(MESSAGE)
  CALL ReportResult(13, 'Missing database file', OK, DETAIL)

  CALL ReadSailDatabase('sail_database.dat', IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. IsSailDatabaseInitialized() .AND. &
    GetSailDatabaseSize() == 33
  CALL ClearSailDatabase()
  CALL ClearSailDatabase()
  OK = OK .AND. .NOT. IsSailDatabaseInitialized() .AND. &
    GetSailDatabaseSize() == 0
  WRITE(DETAIL, '(A,L1,A,I0)') 'initialized after two clears=', &
    IsSailDatabaseInitialized(), ', size=', GetSailDatabaseSize()
  CALL ReportResult(14, 'Repeated database cleanup', OK, DETAIL)

  CALL ReadSailDatabase('sail_database.dat', IERR, MESSAGE)
  OK = IERR == SAIL_OK .AND. GetSailDatabaseSize() == 33
  CALL GetSailDatabaseNode(1, ALPHA, CL, CD, IERR, MESSAGE)
  OK = OK .AND. IERR == SAIL_OK .AND. ABS(ALPHA) <= 1.0E-10_DP .AND. &
    ABS(CL - (-0.29395224_DP)) <= 1.0E-10_DP .AND. &
    ABS(CD - 0.17179593_DP) <= 1.0E-10_DP
  CALL GetSailDatabaseNode(33, ALPHA, CL, CD, IERR, MESSAGE)
  OK = OK .AND. IERR == SAIL_OK .AND. &
    ABS(ALPHA - 180.0_DP) <= 1.0E-10_DP .AND. &
    ABS(CL - 0.26555131_DP) <= 1.0E-10_DP .AND. &
    ABS(CD - 0.17291675_DP) <= 1.0E-10_DP
  WRITE(DETAIL, '(A,L1,A,I0,A,3(ES14.6,1X))') 'reload valid=', OK, &
    ', size=', GetSailDatabaseSize(), ', last node=', ALPHA, CL, CD
  CALL ReportResult(15, 'Reload after cleanup', OK, DETAIL)

  WRITE(*, '(A)') '========================================'
  WRITE(*, '(A)') 'SAILDATABASE TEST SUMMARY'
  WRITE(*, '(A,I0)') 'PASS: ', N_PASS
  WRITE(*, '(A,I0)') 'FAIL: ', N_FAIL
  WRITE(*, '(A)') '========================================'

  IF (N_FAIL > 0) THEN
    STOP 1
  ELSE
    STOP 0
  END IF

CONTAINS

  SUBROUTINE ReportResult(TEST_NUMBER, TEST_NAME, PASSED, TEST_DETAIL)
    INTEGER, INTENT(IN) :: TEST_NUMBER
    CHARACTER(LEN=*), INTENT(IN) :: TEST_NAME
    LOGICAL, INTENT(IN) :: PASSED
    CHARACTER(LEN=*), INTENT(IN) :: TEST_DETAIL

    WRITE(*, '(A)') '--------------------------------------------------'
    WRITE(*, '(A,I2.2,A,A)') 'TEST ', TEST_NUMBER, ': ', TRIM(TEST_NAME)
    IF (PASSED) THEN
      N_PASS = N_PASS + 1
      WRITE(*, '(A)') '[PASS]'
    ELSE
      N_FAIL = N_FAIL + 1
      WRITE(*, '(A)') '[FAIL]'
    END IF
    WRITE(*, '(A)') TRIM(TEST_DETAIL)
  END SUBROUTINE ReportResult

END PROGRAM TEST_SAILDATABASE
