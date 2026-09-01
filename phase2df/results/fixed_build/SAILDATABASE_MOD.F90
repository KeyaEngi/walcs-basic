!===============================================================================
! SAILDATABASE_MOD
!
! Read-only management of the one-dimensional, three-column sail aerodynamic
! database.  Angles are stored in degrees; CL and CD are dimensionless.
!===============================================================================
MODULE SAILDATABASE_MOD
  USE SAILPARAM_MOD, ONLY: &
    DP, DEFAULT_DATABASE_FILE, ALPHA_DB_MIN_DEG, ALPHA_DB_MAX_DEG, &
    ANGLE_TOL_DEG, DATABASE_TOL, SAIL_OK, SAIL_ERR_FILE_NOT_FOUND, &
    SAIL_ERR_FILE_OPEN, SAIL_ERR_FILE_READ, SAIL_ERR_DATABASE_EMPTY, &
    SAIL_ERR_DATABASE_ORDER, SAIL_ERR_DATABASE_DUPLICATE, &
    SAIL_ERR_DATABASE_RANGE, SAIL_ERR_DATABASE_NOT_INITIALIZED, &
    SAIL_ERR_INVALID_INPUT
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  INTEGER, SAVE :: N_DB = 0
  REAL(DP), ALLOCATABLE, SAVE :: ALPHA_DB_DEG(:)
  REAL(DP), ALLOCATABLE, SAVE :: CL_DB(:)
  REAL(DP), ALLOCATABLE, SAVE :: CD_DB(:)
  LOGICAL, SAVE :: DATABASE_INITIALIZED = .FALSE.
  CHARACTER(LEN=1024), SAVE :: LOADED_FILE = ''

  PUBLIC :: ReadSailDatabase
  PUBLIC :: ValidateSailDatabase
  PUBLIC :: ClearSailDatabase
  PUBLIC :: IsSailDatabaseInitialized
  PUBLIC :: GetSailDatabaseSize
  PUBLIC :: GetSailDatabaseNode
  PUBLIC :: GetSailDatabaseBounds

CONTAINS

  ! Read and validate a three-column alpha_deg, CL, CD text database.
  SUBROUTINE ReadSailDatabase(FILE_NAME, IERR, MESSAGE)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: FILE_NAME
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    INTEGER :: IU
    INTEGER :: IOS
    INTEGER :: LINE_NUMBER
    INTEGER :: DATA_COUNT
    INTEGER :: DATA_INDEX
    INTEGER :: ALLOC_STAT
    REAL(DP) :: ALPHA_VALUE
    REAL(DP) :: CL_VALUE
    REAL(DP) :: CD_VALUE
    LOGICAL :: FILE_EXISTS
    LOGICAL :: FILE_OPEN
    LOGICAL :: FIRST_CANDIDATE
    LOGICAL :: VALID_DATA
    CHARACTER(LEN=1024) :: REQUESTED_FILE
    CHARACTER(LEN=1024) :: LINE
    CHARACTER(LEN=2048) :: LOCAL_MESSAGE
    CHARACTER(LEN=2048) :: VALIDATION_MESSAGE

    IERR = SAIL_OK
    MESSAGE = ''
    REQUESTED_FILE = ''
    FILE_OPEN = .FALSE.

    IF (PRESENT(FILE_NAME)) THEN
      IF (LEN_TRIM(FILE_NAME) == 0) THEN
        CALL FailLoad(SAIL_ERR_INVALID_INPUT, &
          'Database file name is empty.')
        RETURN
      END IF
      IF (LEN_TRIM(FILE_NAME) > LEN(REQUESTED_FILE)) THEN
        CALL FailLoad(SAIL_ERR_INVALID_INPUT, &
          'Database file name exceeds the supported path length.')
        RETURN
      END IF
      REQUESTED_FILE = TRIM(FILE_NAME)
    ELSE
      REQUESTED_FILE = DEFAULT_DATABASE_FILE
    END IF

    IF (DATABASE_INITIALIZED) THEN
      IF (TRIM(REQUESTED_FILE) == TRIM(LOADED_FILE)) THEN
        WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
          'Sail database is already initialized from "', &
          TRIM(REQUESTED_FILE), '"; reload skipped, ', N_DB, &
          ' nodes available.'
        MESSAGE = TRIM(LOCAL_MESSAGE)
        RETURN
      END IF
    END IF

    ! First-version reload policy: requesting a different file clears the old
    ! database before the new file is read.  If the new load fails, the old
    ! database is intentionally not retained.
    CALL ClearSailDatabase()

    INQUIRE(FILE=TRIM(REQUESTED_FILE), EXIST=FILE_EXISTS)
    IF (.NOT. FILE_EXISTS) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A)') &
        'Database file not found during existence check: "', &
        TRIM(REQUESTED_FILE), '".'
      CALL FailLoad(SAIL_ERR_FILE_NOT_FOUND, LOCAL_MESSAGE)
      RETURN
    END IF

    OPEN(NEWUNIT=IU, FILE=TRIM(REQUESTED_FILE), STATUS='OLD', &
      ACTION='READ', FORM='FORMATTED', IOSTAT=IOS)
    IF (IOS /= 0) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
        'Could not open database file "', TRIM(REQUESTED_FILE), &
        '" for reading; IOSTAT=', IOS, '.'
      CALL FailLoad(SAIL_ERR_FILE_OPEN, LOCAL_MESSAGE)
      RETURN
    END IF
    FILE_OPEN = .TRUE.

    ! First pass: count valid rows and reject malformed non-comment content.
    DATA_COUNT = 0
    LINE_NUMBER = 0
    FIRST_CANDIDATE = .TRUE.
    DO
      READ(IU, '(A)', IOSTAT=IOS) LINE
      IF (IOS < 0) EXIT
      IF (IOS > 0) THEN
        WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A,I0,A)') &
          'I/O error while counting rows in "', TRIM(REQUESTED_FILE), &
          '" at line ', LINE_NUMBER + 1, '; IOSTAT=', IOS, '.'
        CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
        RETURN
      END IF

      LINE_NUMBER = LINE_NUMBER + 1
      IF (IsBlankOrComment(LINE)) CYCLE

      CALL ParseDataLine(LINE, ALPHA_VALUE, CL_VALUE, CD_VALUE, &
        VALID_DATA)
      IF (VALID_DATA) THEN
        DATA_COUNT = DATA_COUNT + 1
        FIRST_CANDIDATE = .FALSE.
      ELSE IF (FIRST_CANDIDATE .AND. IsHeaderLine(LINE)) THEN
        FIRST_CANDIDATE = .FALSE.
      ELSE
        WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A,A,A)') &
          'Invalid three-column row in "', TRIM(REQUESTED_FILE), &
          '" at line ', LINE_NUMBER, ': "', &
          TRIM(LINE(1:MIN(LEN_TRIM(LINE), 160))), '".'
        CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
        RETURN
      END IF
    END DO

    IF (DATA_COUNT == 0) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A)') &
        'No valid data rows were found in database file "', &
        TRIM(REQUESTED_FILE), '".'
      CALL FailLoad(SAIL_ERR_DATABASE_EMPTY, LOCAL_MESSAGE)
      RETURN
    END IF

    ALLOCATE(ALPHA_DB_DEG(DATA_COUNT), CL_DB(DATA_COUNT), &
      CD_DB(DATA_COUNT), STAT=ALLOC_STAT)
    IF (ALLOC_STAT /= 0) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
        'Memory allocation failed while loading "', &
        TRIM(REQUESTED_FILE), '"; STAT=', ALLOC_STAT, '.'
      CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
      RETURN
    END IF
    N_DB = DATA_COUNT

    REWIND(IU, IOSTAT=IOS)
    IF (IOS /= 0) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
        'Could not rewind database file "', TRIM(REQUESTED_FILE), &
        '"; IOSTAT=', IOS, '.'
      CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
      RETURN
    END IF

    ! Second pass: parse and store exactly the rows accepted above.
    DATA_INDEX = 0
    LINE_NUMBER = 0
    FIRST_CANDIDATE = .TRUE.
    DO
      READ(IU, '(A)', IOSTAT=IOS) LINE
      IF (IOS < 0) EXIT
      IF (IOS > 0) THEN
        WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A,I0,A)') &
          'I/O error while loading "', TRIM(REQUESTED_FILE), &
          '" at line ', LINE_NUMBER + 1, '; IOSTAT=', IOS, '.'
        CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
        RETURN
      END IF

      LINE_NUMBER = LINE_NUMBER + 1
      IF (IsBlankOrComment(LINE)) CYCLE

      CALL ParseDataLine(LINE, ALPHA_VALUE, CL_VALUE, CD_VALUE, &
        VALID_DATA)
      IF (VALID_DATA) THEN
        DATA_INDEX = DATA_INDEX + 1
        IF (DATA_INDEX > N_DB) THEN
          WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A,I0,A)') &
            'Database row count increased during the second read of "', &
            TRIM(REQUESTED_FILE), '" at line ', LINE_NUMBER, &
            '; first-pass node count=', N_DB, '.'
          CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
          RETURN
        END IF
        ALPHA_DB_DEG(DATA_INDEX) = ALPHA_VALUE
        CL_DB(DATA_INDEX) = CL_VALUE
        CD_DB(DATA_INDEX) = CD_VALUE
        FIRST_CANDIDATE = .FALSE.
      ELSE IF (FIRST_CANDIDATE .AND. IsHeaderLine(LINE)) THEN
        FIRST_CANDIDATE = .FALSE.
      ELSE
        WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
          'Database content changed or became invalid while loading "', &
          TRIM(REQUESTED_FILE), '" at line ', LINE_NUMBER, '.'
        CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
        RETURN
      END IF
    END DO

    CLOSE(IU, IOSTAT=IOS)
    FILE_OPEN = .FALSE.
    IF (IOS /= 0) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
        'Could not close database file "', TRIM(REQUESTED_FILE), &
        '" after reading; IOSTAT=', IOS, '.'
      CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
      RETURN
    END IF

    IF (DATA_INDEX /= N_DB) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A,I0,A)') &
        'Database row count decreased during the second read of "', &
        TRIM(REQUESTED_FILE), '"; expected ', N_DB, &
        ' nodes but read ', DATA_INDEX, '.'
      CALL FailLoad(SAIL_ERR_FILE_READ, LOCAL_MESSAGE)
      RETURN
    END IF

    CALL ValidateSailDatabase(IERR, VALIDATION_MESSAGE)
    IF (IERR /= SAIL_OK) THEN
      WRITE(LOCAL_MESSAGE, '(A,A,A,A)') &
        'Database validation failed for "', TRIM(REQUESTED_FILE), &
        '": ', TRIM(VALIDATION_MESSAGE)
      CALL FailLoad(IERR, LOCAL_MESSAGE)
      RETURN
    END IF

    DATABASE_INITIALIZED = .TRUE.
    LOADED_FILE = TRIM(REQUESTED_FILE)
    IERR = SAIL_OK
    WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A,ES16.8,A,ES16.8,A)') &
      'Sail database loaded successfully from "', &
      TRIM(REQUESTED_FILE), '": ', N_DB, ' nodes, alpha range ', &
      ALPHA_DB_DEG(1), ' to ', ALPHA_DB_DEG(N_DB), ' deg.'
    MESSAGE = TRIM(LOCAL_MESSAGE)

  CONTAINS

    SUBROUTINE FailLoad(ERROR_CODE, DETAIL)
      INTEGER, INTENT(IN) :: ERROR_CODE
      CHARACTER(LEN=*), INTENT(IN) :: DETAIL
      INTEGER :: CLOSE_IOS
      CHARACTER(LEN=2048) :: FAILURE_MESSAGE

      CLOSE_IOS = 0
      FAILURE_MESSAGE = TRIM(DETAIL)
      IF (FILE_OPEN) THEN
        CLOSE(IU, IOSTAT=CLOSE_IOS)
        FILE_OPEN = .FALSE.
        IF (CLOSE_IOS /= 0) THEN
          WRITE(FAILURE_MESSAGE, '(A,A,A,I0,A,A)') &
            'Could not close database file "', TRIM(REQUESTED_FILE), &
            '" while cleaning up a failed load; IOSTAT=', CLOSE_IOS, &
            '. Original failure: ', TRIM(DETAIL)
        END IF
      END IF
      CALL ClearSailDatabase()
      IF (CLOSE_IOS /= 0) THEN
        IERR = SAIL_ERR_FILE_READ
      ELSE
        IERR = ERROR_CODE
      END IF
      MESSAGE = TRIM(FAILURE_MESSAGE)
    END SUBROUTINE FailLoad

  END SUBROUTINE ReadSailDatabase


  ! Validate the database currently held in the module's private storage.
  SUBROUTINE ValidateSailDatabase(IERR, MESSAGE)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    INTEGER :: I
    REAL(DP) :: DELTA_ALPHA
    CHARACTER(LEN=1024) :: LOCAL_MESSAGE

    IERR = SAIL_OK
    MESSAGE = ''

    IF (.NOT. ALLOCATED(ALPHA_DB_DEG) .OR. &
        .NOT. ALLOCATED(CL_DB) .OR. .NOT. ALLOCATED(CD_DB)) THEN
      IERR = SAIL_ERR_DATABASE_NOT_INITIALIZED
      MESSAGE = 'Sail database storage is not fully allocated.'
      RETURN
    END IF

    IF (N_DB < 2) THEN
      IERR = SAIL_ERR_DATABASE_EMPTY
      MESSAGE = 'Sail database must contain at least two nodes.'
      RETURN
    END IF

    IF (SIZE(ALPHA_DB_DEG) /= N_DB .OR. SIZE(CL_DB) /= N_DB .OR. &
        SIZE(CD_DB) /= N_DB) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Sail database array lengths do not match N_DB.'
      RETURN
    END IF

    DO I = 1, N_DB
      IF (.NOT. IEEE_IS_FINITE(ALPHA_DB_DEG(I)) .OR. &
          .NOT. IEEE_IS_FINITE(CL_DB(I)) .OR. &
          .NOT. IEEE_IS_FINITE(CD_DB(I))) THEN
        IERR = SAIL_ERR_INVALID_INPUT
        WRITE(LOCAL_MESSAGE, '(A,I0,A)') &
          'Non-finite alpha, CL, or CD value at database node ', I, '.'
        MESSAGE = TRIM(LOCAL_MESSAGE)
        RETURN
      END IF

      IF (CD_DB(I) < -DATABASE_TOL) THEN
        IERR = SAIL_ERR_INVALID_INPUT
        WRITE(LOCAL_MESSAGE, '(A,I0,A,ES16.8,A)') &
          'Negative CD at database node ', I, ': ', CD_DB(I), '.'
        MESSAGE = TRIM(LOCAL_MESSAGE)
        RETURN
      END IF
    END DO

    DO I = 1, N_DB - 1
      DELTA_ALPHA = ALPHA_DB_DEG(I + 1) - ALPHA_DB_DEG(I)
      IF (DELTA_ALPHA < -DATABASE_TOL) THEN
        IERR = SAIL_ERR_DATABASE_ORDER
        WRITE(LOCAL_MESSAGE, '(A,I0,A,I0,A)') &
          'Attack angles are out of order between nodes ', I, &
          ' and ', I + 1, '.'
        MESSAGE = TRIM(LOCAL_MESSAGE)
        RETURN
      END IF
      IF (ABS(DELTA_ALPHA) <= DATABASE_TOL) THEN
        IERR = SAIL_ERR_DATABASE_DUPLICATE
        WRITE(LOCAL_MESSAGE, '(A,I0,A,I0,A)') &
          'Duplicate attack angle at nodes ', I, ' and ', I + 1, '.'
        MESSAGE = TRIM(LOCAL_MESSAGE)
        RETURN
      END IF
    END DO

    IF (ABS(ALPHA_DB_DEG(1) - ALPHA_DB_MIN_DEG) > ANGLE_TOL_DEG .OR. &
        ABS(ALPHA_DB_DEG(N_DB) - ALPHA_DB_MAX_DEG) > ANGLE_TOL_DEG) THEN
      IERR = SAIL_ERR_DATABASE_RANGE
      WRITE(LOCAL_MESSAGE, '(A,ES16.8,A,ES16.8,A)') &
        'Database endpoints do not cover the required range: ', &
        ALPHA_DB_DEG(1), ' to ', ALPHA_DB_DEG(N_DB), ' deg.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    MESSAGE = 'Sail database validation completed successfully.'
  END SUBROUTINE ValidateSailDatabase


  ! Release all private database storage. Repeated calls are safe.
  SUBROUTINE ClearSailDatabase()
    IF (ALLOCATED(ALPHA_DB_DEG)) DEALLOCATE(ALPHA_DB_DEG)
    IF (ALLOCATED(CL_DB)) DEALLOCATE(CL_DB)
    IF (ALLOCATED(CD_DB)) DEALLOCATE(CD_DB)
    N_DB = 0
    DATABASE_INITIALIZED = .FALSE.
    LOADED_FILE = ''
  END SUBROUTINE ClearSailDatabase


  ! Report whether a database has been loaded and validated successfully.
  LOGICAL FUNCTION IsSailDatabaseInitialized()
    IsSailDatabaseInitialized = DATABASE_INITIALIZED
  END FUNCTION IsSailDatabaseInitialized


  ! Return the number of nodes, or zero when no database is initialized.
  INTEGER FUNCTION GetSailDatabaseSize()
    IF (DATABASE_INITIALIZED) THEN
      GetSailDatabaseSize = N_DB
    ELSE
      GetSailDatabaseSize = 0
    END IF
  END FUNCTION GetSailDatabaseSize


  ! Copy one database node to caller-owned scalar outputs.
  SUBROUTINE GetSailDatabaseNode(INDEX_NODE, ALPHA_DEG, CL, CD, &
      IERR, MESSAGE)
    INTEGER, INTENT(IN) :: INDEX_NODE
    REAL(DP), INTENT(OUT) :: ALPHA_DEG
    REAL(DP), INTENT(OUT) :: CL
    REAL(DP), INTENT(OUT) :: CD
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    CHARACTER(LEN=1024) :: LOCAL_MESSAGE

    ALPHA_DEG = 0.0_DP
    CL = 0.0_DP
    CD = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''

    IF (.NOT. DATABASE_INITIALIZED) THEN
      IERR = SAIL_ERR_DATABASE_NOT_INITIALIZED
      MESSAGE = 'Sail database is not initialized.'
      RETURN
    END IF

    IF (INDEX_NODE < 1 .OR. INDEX_NODE > N_DB) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      WRITE(LOCAL_MESSAGE, '(A,I0,A,I0,A)') &
        'Database node index ', INDEX_NODE, &
        ' is outside the valid range 1 to ', N_DB, '.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    ALPHA_DEG = ALPHA_DB_DEG(INDEX_NODE)
    CL = CL_DB(INDEX_NODE)
    CD = CD_DB(INDEX_NODE)
    MESSAGE = 'Sail database node returned successfully.'
  END SUBROUTINE GetSailDatabaseNode


  ! Return the first and last attack-angle nodes in degrees.
  SUBROUTINE GetSailDatabaseBounds(ALPHA_MIN, ALPHA_MAX, IERR, MESSAGE)
    REAL(DP), INTENT(OUT) :: ALPHA_MIN
    REAL(DP), INTENT(OUT) :: ALPHA_MAX
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    ALPHA_MIN = 0.0_DP
    ALPHA_MAX = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''

    IF (.NOT. DATABASE_INITIALIZED) THEN
      IERR = SAIL_ERR_DATABASE_NOT_INITIALIZED
      MESSAGE = 'Sail database is not initialized.'
      RETURN
    END IF

    ALPHA_MIN = ALPHA_DB_DEG(1)
    ALPHA_MAX = ALPHA_DB_DEG(N_DB)
    MESSAGE = 'Sail database bounds returned successfully.'
  END SUBROUTINE GetSailDatabaseBounds


  LOGICAL FUNCTION IsBlankOrComment(LINE)
    CHARACTER(LEN=*), INTENT(IN) :: LINE
    CHARACTER(LEN=LEN(LINE)) :: ADJUSTED_LINE

    ADJUSTED_LINE = ADJUSTL(LINE)
    IsBlankOrComment = LEN_TRIM(ADJUSTED_LINE) == 0
    IF (.NOT. IsBlankOrComment) THEN
      IsBlankOrComment = ADJUSTED_LINE(1:1) == '!' .OR. &
        ADJUSTED_LINE(1:1) == '#'
    END IF
  END FUNCTION IsBlankOrComment


  LOGICAL FUNCTION IsHeaderLine(LINE)
    CHARACTER(LEN=*), INTENT(IN) :: LINE
    CHARACTER(LEN=LEN(LINE)) :: LOWER_LINE

    LOWER_LINE = ToLowerAscii(LINE)
    IsHeaderLine = INDEX(LOWER_LINE, 'alpha') > 0 .AND. &
      INDEX(LOWER_LINE, 'cl') > 0 .AND. INDEX(LOWER_LINE, 'cd') > 0
  END FUNCTION IsHeaderLine


  SUBROUTINE ParseDataLine(LINE, ALPHA_VALUE, CL_VALUE, CD_VALUE, VALID)
    CHARACTER(LEN=*), INTENT(IN) :: LINE
    REAL(DP), INTENT(OUT) :: ALPHA_VALUE
    REAL(DP), INTENT(OUT) :: CL_VALUE
    REAL(DP), INTENT(OUT) :: CD_VALUE
    LOGICAL, INTENT(OUT) :: VALID

    INTEGER :: IOS
    INTEGER :: EXTRA_IOS
    CHARACTER(LEN=64) :: EXTRA_TOKEN

    ALPHA_VALUE = 0.0_DP
    CL_VALUE = 0.0_DP
    CD_VALUE = 0.0_DP
    EXTRA_TOKEN = ''

    ! Comments are accepted only as whole lines by IsBlankOrComment.  Reject
    ! inline comment text here so that every data row contains exactly three
    ! fields and no trailing content.
    IF (INDEX(LINE, '!') > 0 .OR. INDEX(LINE, '#') > 0) THEN
      VALID = .FALSE.
      RETURN
    END IF

    READ(LINE, *, IOSTAT=IOS) ALPHA_VALUE, CL_VALUE, CD_VALUE
    IF (IOS /= 0) THEN
      VALID = .FALSE.
      RETURN
    END IF

    READ(LINE, *, IOSTAT=EXTRA_IOS) ALPHA_VALUE, CL_VALUE, CD_VALUE, &
      EXTRA_TOKEN
    VALID = EXTRA_IOS < 0
  END SUBROUTINE ParseDataLine


  FUNCTION ToLowerAscii(TEXT) RESULT(LOWER_TEXT)
    CHARACTER(LEN=*), INTENT(IN) :: TEXT
    CHARACTER(LEN=LEN(TEXT)) :: LOWER_TEXT
    INTEGER :: I
    INTEGER :: CODE

    LOWER_TEXT = TEXT
    DO I = 1, LEN(TEXT)
      CODE = IACHAR(TEXT(I:I))
      IF (CODE >= IACHAR('A') .AND. CODE <= IACHAR('Z')) THEN
        LOWER_TEXT(I:I) = ACHAR(CODE + IACHAR('a') - IACHAR('A'))
      END IF
    END DO
  END FUNCTION ToLowerAscii

END MODULE SAILDATABASE_MOD
