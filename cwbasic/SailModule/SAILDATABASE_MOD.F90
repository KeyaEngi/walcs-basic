MODULE SAILDATABASE_MOD
!------------------------------------------------------------------------------
! Sail aerodynamic coefficient database.
!
! Expected file format:
!   alpha_deg   CL   CD
! Lines beginning with # or ! are ignored. Blank lines are ignored.
!------------------------------------------------------------------------------
USE SAILPARAM_MOD, ONLY: SAIL_DP, SAIL_MAX_DB_POINTS

IMPLICIT NONE
SAVE
PRIVATE

INTEGER :: SailDbN = 0
REAL(SAIL_DP), DIMENSION(SAIL_MAX_DB_POINTS) :: SailDbAlphaDeg = 0.0_SAIL_DP
REAL(SAIL_DP), DIMENSION(SAIL_MAX_DB_POINTS) :: SailDbCL = 0.0_SAIL_DP
REAL(SAIL_DP), DIMENSION(SAIL_MAX_DB_POINTS) :: SailDbCD = 0.0_SAIL_DP

PUBLIC :: SailDbN
PUBLIC :: SailDbAlphaDeg
PUBLIC :: SailDbCL
PUBLIC :: SailDbCD
PUBLIC :: ReadSailDatabase

CONTAINS

SUBROUTINE ReadSailDatabase(fileName, ierr)
IMPLICIT NONE
CHARACTER(LEN=*), INTENT(IN) :: fileName
INTEGER, INTENT(OUT) :: ierr

INTEGER :: ios
INTEGER :: unitNo
REAL(SAIL_DP) :: alphaDeg
REAL(SAIL_DP) :: clValue
REAL(SAIL_DP) :: cdValue
CHARACTER(LEN=256) :: line

ierr = 0
unitNo = 771
SailDbN = 0
SailDbAlphaDeg = 0.0_SAIL_DP
SailDbCL = 0.0_SAIL_DP
SailDbCD = 0.0_SAIL_DP

OPEN(unitNo, FILE=fileName, STATUS='OLD', ACTION='READ', IOSTAT=ios)
IF (ios /= 0) THEN
    ierr = 1
    RETURN
END IF

DO
    READ(unitNo, '(A)', IOSTAT=ios) line
    IF (ios /= 0) EXIT
    IF (LEN_TRIM(ADJUSTL(line)) == 0) CYCLE
    IF (line(1:1) == '#' .OR. line(1:1) == '!') CYCLE

    READ(line, *, IOSTAT=ios) alphaDeg, clValue, cdValue
    IF (ios /= 0) THEN
        ierr = 2
        EXIT
    END IF

    IF (SailDbN >= SAIL_MAX_DB_POINTS) THEN
        ierr = 3
        EXIT
    END IF

    SailDbN = SailDbN + 1
    SailDbAlphaDeg(SailDbN) = alphaDeg
    SailDbCL(SailDbN) = clValue
    SailDbCD(SailDbN) = cdValue
END DO

CLOSE(unitNo)

IF (ierr == 0 .AND. SailDbN <= 0) ierr = 4

END SUBROUTINE ReadSailDatabase

END MODULE SAILDATABASE_MOD
