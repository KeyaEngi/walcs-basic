!===============================================================================
! SAILFORCE_MOD
!
! Convert horizontal relative wind and aerodynamic coefficients into body-axis
! sail force, moment about the vessel CG, and six-degree-of-freedom loads.
!===============================================================================
MODULE SAILFORCE_MOD
  USE SAILPARAM_MOD, ONLY: &
    DP, IDX_X, IDX_Y, IDX_Z, IDX_FX, IDX_FY, IDX_FZ, &
    IDX_MX, IDX_MY, IDX_MZ, N_DOF, RHO_AIR, SAIL_AREA, VREL_MIN, &
    DATABASE_TOL, SAIL_OK, SAIL_ERR_INVALID_INPUT, &
    SAIL_ERR_LOW_WIND_SPEED
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: ComputeSailForce

CONTAINS

  SUBROUTINE ComputeSailForce( &
      V_REL_H_BODY, V_REL_H_MAG, CL, CD, R_SAIL_BODY, &
      E_DRAG_BODY, E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, &
      LOAD_6DOF, IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: V_REL_H_BODY(3)
    REAL(DP), INTENT(IN) :: V_REL_H_MAG
    REAL(DP), INTENT(IN) :: CL
    REAL(DP), INTENT(IN) :: CD
    REAL(DP), INTENT(IN) :: R_SAIL_BODY(3)
    REAL(DP), INTENT(OUT) :: E_DRAG_BODY(3)
    REAL(DP), INTENT(OUT) :: E_LIFT_BODY(3)
    REAL(DP), INTENT(OUT) :: Q_DYNAMIC
    REAL(DP), INTENT(OUT) :: FORCE_BODY(3)
    REAL(DP), INTENT(OUT) :: MOMENT_BODY(3)
    REAL(DP), INTENT(OUT) :: LOAD_6DOF(N_DOF)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    REAL(DP), PARAMETER :: DIRECTION_TOL = 1.0E-10_DP
    REAL(DP) :: V_MAG_CHECK
    REAL(DP) :: MAG_TOL
    REAL(DP) :: DRAG_NORM
    REAL(DP) :: LIFT_NORM
    REAL(DP) :: DIRECTION_DOT
    CHARACTER(LEN=512) :: DIAGNOSTIC

    CALL ZeroOutputs(E_DRAG_BODY, E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, &
      MOMENT_BODY, LOAD_6DOF)
    IERR = SAIL_OK
    MESSAGE = ''

    IF (.NOT. ALL(IEEE_IS_FINITE(V_REL_H_BODY))) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind vector contains NaN or infinity.'
      RETURN
    END IF
    IF (.NOT. IEEE_IS_FINITE(V_REL_H_MAG)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind speed is NaN or infinity.'
      RETURN
    END IF
    IF (.NOT. IEEE_IS_FINITE(CL)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: lift coefficient contains NaN or infinity.'
      RETURN
    END IF
    IF (.NOT. IEEE_IS_FINITE(CD)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: drag coefficient contains NaN or infinity.'
      RETURN
    END IF
    IF (.NOT. ALL(IEEE_IS_FINITE(R_SAIL_BODY))) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: sail position vector contains NaN or infinity.'
      RETURN
    END IF
    IF (V_REL_H_MAG < 0.0_DP) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind speed must be nonnegative.'
      RETURN
    END IF

    V_MAG_CHECK = SQRT(V_REL_H_BODY(IDX_X)**2 + &
                       V_REL_H_BODY(IDX_Y)**2)
    IF (.NOT. IEEE_IS_FINITE(V_MAG_CHECK)) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind vector magnitude overflowed.'
      RETURN
    END IF
    MAG_TOL = MAX(1.0E-12_DP, &
                  1.0E-10_DP * MAX(1.0_DP, V_MAG_CHECK))

    IF (ABS(V_REL_H_BODY(IDX_Z)) > MAG_TOL) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: horizontal relative wind has a nonzero z component.'
      RETURN
    END IF
    IF (ABS(V_REL_H_MAG - V_MAG_CHECK) > MAG_TOL) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: relative-wind scalar and vector magnitudes disagree.'
      RETURN
    END IF
    IF (CD < -DATABASE_TOL) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid input: drag coefficient is negative beyond tolerance.'
      RETURN
    END IF
    IF (V_REL_H_MAG < VREL_MIN) THEN
      IERR = SAIL_ERR_LOW_WIND_SPEED
      MESSAGE = 'Low horizontal relative wind; all sail loads were set to zero.'
      RETURN
    END IF

    E_DRAG_BODY(IDX_X) = V_REL_H_BODY(IDX_X) / V_REL_H_MAG
    E_DRAG_BODY(IDX_Y) = V_REL_H_BODY(IDX_Y) / V_REL_H_MAG
    E_DRAG_BODY(IDX_Z) = 0.0_DP

    ! E_LIFT_BODY = E_DRAG_BODY cross E_Z_BODY.
    E_LIFT_BODY(IDX_X) = E_DRAG_BODY(IDX_Y)
    E_LIFT_BODY(IDX_Y) = -E_DRAG_BODY(IDX_X)
    E_LIFT_BODY(IDX_Z) = 0.0_DP

    DRAG_NORM = VectorNorm2D(E_DRAG_BODY)
    LIFT_NORM = VectorNorm2D(E_LIFT_BODY)
    DIRECTION_DOT = DOT_PRODUCT(E_DRAG_BODY, E_LIFT_BODY)
    IF (ABS(DRAG_NORM - 1.0_DP) > DIRECTION_TOL .OR. &
        ABS(LIFT_NORM - 1.0_DP) > DIRECTION_TOL .OR. &
        ABS(DIRECTION_DOT) > DIRECTION_TOL) THEN
      CALL ZeroOutputs(E_DRAG_BODY, E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, &
        MOMENT_BODY, LOAD_6DOF)
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid result: drag/lift direction construction failed.'
      RETURN
    END IF

    Q_DYNAMIC = 0.5_DP * RHO_AIR * V_REL_H_MAG**2
    FORCE_BODY = Q_DYNAMIC * SAIL_AREA * &
      (CD * E_DRAG_BODY + CL * E_LIFT_BODY)
    FORCE_BODY(IDX_Z) = 0.0_DP

    CALL CrossProduct3D(R_SAIL_BODY, FORCE_BODY, MOMENT_BODY)

    LOAD_6DOF(IDX_FX) = FORCE_BODY(IDX_X)
    LOAD_6DOF(IDX_FY) = FORCE_BODY(IDX_Y)
    LOAD_6DOF(IDX_FZ) = FORCE_BODY(IDX_Z)
    LOAD_6DOF(IDX_MX) = MOMENT_BODY(IDX_X)
    LOAD_6DOF(IDX_MY) = MOMENT_BODY(IDX_Y)
    LOAD_6DOF(IDX_MZ) = MOMENT_BODY(IDX_Z)

    IF (.NOT. OutputsAreFinite(E_DRAG_BODY, E_LIFT_BODY, Q_DYNAMIC, &
        FORCE_BODY, MOMENT_BODY, LOAD_6DOF)) THEN
      CALL ZeroOutputs(E_DRAG_BODY, E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, &
        MOMENT_BODY, LOAD_6DOF)
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Invalid result: sail-force computation overflowed or is non-finite.'
      RETURN
    END IF

    WRITE(DIAGNOSTIC, '(A,ES14.6,A,ES14.6,A,ES14.6,A,ES14.6,' // &
      'A,3(ES14.6,1X),A,3(ES14.6,1X),A)') &
      'Sail force computed successfully; speed=', V_REL_H_MAG, &
      ', CL=', CL, ', CD=', CD, ', q=', Q_DYNAMIC, &
      ', force=(', FORCE_BODY, '), moment=(', MOMENT_BODY, ').'
    MESSAGE = TRIM(DIAGNOSTIC)
  END SUBROUTINE ComputeSailForce


  PURE SUBROUTINE CrossProduct3D(A, B, C)
    REAL(DP), INTENT(IN) :: A(3)
    REAL(DP), INTENT(IN) :: B(3)
    REAL(DP), INTENT(OUT) :: C(3)

    C(IDX_X) = A(IDX_Y) * B(IDX_Z) - A(IDX_Z) * B(IDX_Y)
    C(IDX_Y) = A(IDX_Z) * B(IDX_X) - A(IDX_X) * B(IDX_Z)
    C(IDX_Z) = A(IDX_X) * B(IDX_Y) - A(IDX_Y) * B(IDX_X)
  END SUBROUTINE CrossProduct3D


  PURE FUNCTION VectorNorm2D(VECTOR) RESULT(NORM_VALUE)
    REAL(DP), INTENT(IN) :: VECTOR(3)
    REAL(DP) :: NORM_VALUE

    NORM_VALUE = SQRT(VECTOR(IDX_X)**2 + VECTOR(IDX_Y)**2)
  END FUNCTION VectorNorm2D


  LOGICAL FUNCTION OutputsAreFinite(E_DRAG, E_LIFT, Q, FORCE, MOMENT, LOAD)
    REAL(DP), INTENT(IN) :: E_DRAG(3)
    REAL(DP), INTENT(IN) :: E_LIFT(3)
    REAL(DP), INTENT(IN) :: Q
    REAL(DP), INTENT(IN) :: FORCE(3)
    REAL(DP), INTENT(IN) :: MOMENT(3)
    REAL(DP), INTENT(IN) :: LOAD(N_DOF)

    OutputsAreFinite = ALL(IEEE_IS_FINITE(E_DRAG)) .AND. &
      ALL(IEEE_IS_FINITE(E_LIFT)) .AND. IEEE_IS_FINITE(Q) .AND. &
      ALL(IEEE_IS_FINITE(FORCE)) .AND. ALL(IEEE_IS_FINITE(MOMENT)) .AND. &
      ALL(IEEE_IS_FINITE(LOAD))
  END FUNCTION OutputsAreFinite


  PURE SUBROUTINE ZeroOutputs(E_DRAG, E_LIFT, Q, FORCE, MOMENT, LOAD)
    REAL(DP), INTENT(OUT) :: E_DRAG(3)
    REAL(DP), INTENT(OUT) :: E_LIFT(3)
    REAL(DP), INTENT(OUT) :: Q
    REAL(DP), INTENT(OUT) :: FORCE(3)
    REAL(DP), INTENT(OUT) :: MOMENT(3)
    REAL(DP), INTENT(OUT) :: LOAD(N_DOF)

    E_DRAG = 0.0_DP
    E_LIFT = 0.0_DP
    Q = 0.0_DP
    FORCE = 0.0_DP
    MOMENT = 0.0_DP
    LOAD = 0.0_DP
  END SUBROUTINE ZeroOutputs

END MODULE SAILFORCE_MOD
