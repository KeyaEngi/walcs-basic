!===============================================================================
! Relative wind at the sail aerodynamic center in body-fixed coordinates.
!===============================================================================
MODULE SAILRELWIND_MOD
  USE SAILPARAM_MOD, ONLY: &
    DP, IDX_X, IDX_Y, IDX_Z, VREL_MIN, SAIL_OK, &
    SAIL_ERR_INVALID_INPUT, SAIL_ERR_LOW_WIND_SPEED
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeSailRelativeWind

CONTAINS

  SUBROUTINE ComputeSailRelativeWind( &
      V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, &
      V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, &
      V_REL_H_MAG, IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: V_WIND_BODY(3)
    REAL(DP), INTENT(IN) :: V_CG_BODY(3)
    REAL(DP), INTENT(IN) :: OMEGA_BODY(3)
    REAL(DP), INTENT(IN) :: R_SAIL_BODY(3)
    REAL(DP), INTENT(OUT) :: V_SAIL_POINT_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_H_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_H_MAG
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    REAL(DP) :: OMEGA_CROSS_R(3)
    CHARACTER(LEN=512) :: LOCAL_MESSAGE
    INTEGER :: COMPONENT

    V_SAIL_POINT_BODY = 0.0_DP
    V_REL_BODY = 0.0_DP
    V_REL_H_BODY = 0.0_DP
    V_REL_H_MAG = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''
    LOCAL_MESSAGE = ''

    COMPONENT = FirstNonFinite(V_WIND_BODY)
    IF (COMPONENT > 0) THEN
      CALL SetInvalidInput('V_WIND_BODY', COMPONENT, IERR, MESSAGE)
      RETURN
    END IF

    COMPONENT = FirstNonFinite(V_CG_BODY)
    IF (COMPONENT > 0) THEN
      CALL SetInvalidInput('V_CG_BODY', COMPONENT, IERR, MESSAGE)
      RETURN
    END IF

    COMPONENT = FirstNonFinite(OMEGA_BODY)
    IF (COMPONENT > 0) THEN
      CALL SetInvalidInput('OMEGA_BODY', COMPONENT, IERR, MESSAGE)
      RETURN
    END IF

    COMPONENT = FirstNonFinite(R_SAIL_BODY)
    IF (COMPONENT > 0) THEN
      CALL SetInvalidInput('R_SAIL_BODY', COMPONENT, IERR, MESSAGE)
      RETURN
    END IF

    CALL CrossProduct3D(OMEGA_BODY, R_SAIL_BODY, OMEGA_CROSS_R)
    V_SAIL_POINT_BODY = V_CG_BODY + OMEGA_CROSS_R
    V_REL_BODY = V_WIND_BODY - V_SAIL_POINT_BODY

    V_REL_H_BODY(IDX_X) = V_REL_BODY(IDX_X)
    V_REL_H_BODY(IDX_Y) = V_REL_BODY(IDX_Y)
    V_REL_H_BODY(IDX_Z) = 0.0_DP
    V_REL_H_MAG = SQRT( &
      V_REL_H_BODY(IDX_X)**2 + V_REL_H_BODY(IDX_Y)**2)

    IF (.NOT. OutputsAreFinite(V_SAIL_POINT_BODY, V_REL_BODY, &
        V_REL_H_BODY, V_REL_H_MAG)) THEN
      V_SAIL_POINT_BODY = 0.0_DP
      V_REL_BODY = 0.0_DP
      V_REL_H_BODY = 0.0_DP
      V_REL_H_MAG = 0.0_DP
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Computed relative-wind output is not finite.'
      RETURN
    END IF

    IF (V_REL_H_MAG < VREL_MIN) THEN
      IERR = SAIL_ERR_LOW_WIND_SPEED
      WRITE(LOCAL_MESSAGE, '(A,ES16.8,A)') &
        'Horizontal relative wind is below VREL_MIN; speed = ', &
        V_REL_H_MAG, ' m/s.'
      MESSAGE = TRIM(LOCAL_MESSAGE)
      RETURN
    END IF

    WRITE(LOCAL_MESSAGE, '(A,ES16.8,A,3(ES16.8,1X),A,3(ES16.8,1X))') &
      'Sail relative wind computed successfully; horizontal speed = ', &
      V_REL_H_MAG, ' m/s; V_REL_BODY = ', V_REL_BODY, &
      '; V_SAIL_POINT_BODY = ', V_SAIL_POINT_BODY
    MESSAGE = TRIM(LOCAL_MESSAGE)
  END SUBROUTINE ComputeSailRelativeWind

  SUBROUTINE CrossProduct3D(A, B, C)
    REAL(DP), INTENT(IN) :: A(3)
    REAL(DP), INTENT(IN) :: B(3)
    REAL(DP), INTENT(OUT) :: C(3)

    C(IDX_X) = A(IDX_Y) * B(IDX_Z) - A(IDX_Z) * B(IDX_Y)
    C(IDX_Y) = A(IDX_Z) * B(IDX_X) - A(IDX_X) * B(IDX_Z)
    C(IDX_Z) = A(IDX_X) * B(IDX_Y) - A(IDX_Y) * B(IDX_X)
  END SUBROUTINE CrossProduct3D

  INTEGER FUNCTION FirstNonFinite(VECTOR)
    REAL(DP), INTENT(IN) :: VECTOR(3)
    INTEGER :: I

    FirstNonFinite = 0
    DO I = 1, 3
      IF (.NOT. IEEE_IS_FINITE(VECTOR(I))) THEN
        FirstNonFinite = I
        RETURN
      END IF
    END DO
  END FUNCTION FirstNonFinite

  SUBROUTINE SetInvalidInput(NAME, COMPONENT, IERR, MESSAGE)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    INTEGER, INTENT(IN) :: COMPONENT
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE
    CHARACTER(LEN=256) :: LOCAL_MESSAGE

    IERR = SAIL_ERR_INVALID_INPUT
    WRITE(LOCAL_MESSAGE, '(A,A,A,I0,A)') &
      'Invalid input: ', TRIM(NAME), ' component ', COMPONENT, &
      ' is not finite.'
    MESSAGE = TRIM(LOCAL_MESSAGE)
  END SUBROUTINE SetInvalidInput

  LOGICAL FUNCTION OutputsAreFinite( &
      V_POINT, V_REL, V_REL_H, V_REL_H_MAG)
    REAL(DP), INTENT(IN) :: V_POINT(3)
    REAL(DP), INTENT(IN) :: V_REL(3)
    REAL(DP), INTENT(IN) :: V_REL_H(3)
    REAL(DP), INTENT(IN) :: V_REL_H_MAG

    OutputsAreFinite = FirstNonFinite(V_POINT) == 0 .AND. &
      FirstNonFinite(V_REL) == 0 .AND. &
      FirstNonFinite(V_REL_H) == 0 .AND. &
      IEEE_IS_FINITE(V_REL_H_MAG)
  END FUNCTION OutputsAreFinite

END MODULE SAILRELWIND_MOD
