! Minimal, non-intrusive adapter between planning and the accepted SailModule.
MODULE SAILPLANNING_ADAPTER_MOD
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE SAILPARAM_MOD, ONLY: DP, SAIL_OK, SAIL_ERR_LOW_WIND_SPEED, &
    SAIL_ERR_INVALID_INPUT
  USE SAILMODULE_API_MOD, ONLY: InitializeSailModule, ComputeSailModuleLoads, &
    FinalizeSailModule
  USE ArrayOperations, ONLY: VectorG2L, VectorL2G
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: InitializeSailPlanningAdapter
  PUBLIC :: ComputeSailPlanningDryRun
  PUBLIC :: TransformSailLoadToPlanning
  PUBLIC :: FinalizeSailPlanningAdapter

  LOGICAL :: ADAPTER_INITIALIZED = .FALSE.

CONTAINS

  SUBROUTINE TransformSailLoadToPlanning( &
      ANGLES, LOAD_6DOF_BODY, LOAD_6DOF_PLANNING, IERR, MESSAGE)
    REAL(DP), INTENT(IN) :: ANGLES(3)
    REAL(DP), INTENT(IN) :: LOAD_6DOF_BODY(6)
    REAL(DP), INTENT(OUT) :: LOAD_6DOF_PLANNING(6)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    LOAD_6DOF_PLANNING = 0.0_DP
    IF (.NOT. ALL(IEEE_IS_FINITE(ANGLES)) .OR. &
        .NOT. ALL(IEEE_IS_FINITE(LOAD_6DOF_BODY))) THEN
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Non-finite Sail load transformation input.'
      RETURN
    END IF

    ! planning rigid modes 1:3 are translations in its translating frame;
    ! modes 4:6 are moments about CG in that same frame. SailModule returns
    ! both vectors in body axes and its moment is already about CG, so rotate
    ! force and moment independently and do not add another r cross F term.
    LOAD_6DOF_PLANNING(1:3) = VectorL2G(LOAD_6DOF_BODY(1:3), ANGLES)
    LOAD_6DOF_PLANNING(4:6) = VectorL2G(LOAD_6DOF_BODY(4:6), ANGLES)

    IF (.NOT. ALL(IEEE_IS_FINITE(LOAD_6DOF_PLANNING))) THEN
      LOAD_6DOF_PLANNING = 0.0_DP
      IERR = SAIL_ERR_INVALID_INPUT
      MESSAGE = 'Non-finite Sail load transformation output.'
      RETURN
    END IF
    IERR = SAIL_OK
    MESSAGE = ''
  END SUBROUTINE TransformSailLoadToPlanning
  SUBROUTINE InitializeSailPlanningAdapter(DATABASE_FILE, IERR, MESSAGE)
    CHARACTER(LEN=*), INTENT(IN) :: DATABASE_FILE
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    IF (ADAPTER_INITIALIZED) THEN
      IERR = SAIL_OK
      MESSAGE = 'Sail planning adapter was already initialized.'
      RETURN
    END IF

    CALL InitializeSailModule(TRIM(DATABASE_FILE), IERR, MESSAGE)
    IF (IERR == SAIL_OK) ADAPTER_INITIALIZED = .TRUE.
  END SUBROUTINE InitializeSailPlanningAdapter


  SUBROUTINE FinalizeSailPlanningAdapter(IERR, MESSAGE)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE

    IF (.NOT. ADAPTER_INITIALIZED) THEN
      IERR = SAIL_OK
      MESSAGE = 'Sail planning adapter was not initialized; nothing to finalize.'
      RETURN
    END IF

    CALL FinalizeSailModule(IERR, MESSAGE)
    IF (IERR == SAIL_OK) ADAPTER_INITIALIZED = .FALSE.
  END SUBROUTINE FinalizeSailPlanningAdapter


  SUBROUTINE ComputeSailPlanningDryRun( &
      U0, Y_STATE, NR, STAGE_TIME, &
      V_WIND_GLOBAL, R_SAIL_BODY, DELTA_S_DEG, &
      V_CG_BODY, OMEGA_BODY, V_WIND_BODY, &
      V_REL_H_MAG, ALPHA_DB_DEG, CL, CD, &
      FORCE_BODY, MOMENT_BODY, LOAD_6DOF, &
      IERR, MESSAGE, RK_STAGE, PRINT_DIAGNOSTICS, &
      V_REL_BODY_OUT, E_UPSTREAM_BODY_OUT, C_CHORD_BODY_OUT, &
      ALPHA_RAW_DEG_OUT, Q_DYNAMIC_OUT)
    INTEGER, INTENT(IN) :: NR
    REAL(DP), INTENT(IN) :: U0
    REAL(DP), INTENT(IN) :: Y_STATE(2*NR)
    REAL(DP), INTENT(IN) :: STAGE_TIME
    REAL(DP), INTENT(IN) :: V_WIND_GLOBAL(3)
    REAL(DP), INTENT(IN) :: R_SAIL_BODY(3)
    REAL(DP), INTENT(IN) :: DELTA_S_DEG
    REAL(DP), INTENT(OUT) :: V_CG_BODY(3)
    REAL(DP), INTENT(OUT) :: OMEGA_BODY(3)
    REAL(DP), INTENT(OUT) :: V_WIND_BODY(3)
    REAL(DP), INTENT(OUT) :: V_REL_H_MAG
    REAL(DP), INTENT(OUT) :: ALPHA_DB_DEG
    REAL(DP), INTENT(OUT) :: CL
    REAL(DP), INTENT(OUT) :: CD
    REAL(DP), INTENT(OUT) :: FORCE_BODY(3)
    REAL(DP), INTENT(OUT) :: MOMENT_BODY(3)
    REAL(DP), INTENT(OUT) :: LOAD_6DOF(6)
    INTEGER, INTENT(OUT) :: IERR
    CHARACTER(LEN=*), INTENT(OUT) :: MESSAGE
    INTEGER, INTENT(IN), OPTIONAL :: RK_STAGE
    LOGICAL, INTENT(IN), OPTIONAL :: PRINT_DIAGNOSTICS
    REAL(DP), INTENT(OUT), OPTIONAL :: V_REL_BODY_OUT(3)
    REAL(DP), INTENT(OUT), OPTIONAL :: E_UPSTREAM_BODY_OUT(3)
    REAL(DP), INTENT(OUT), OPTIONAL :: C_CHORD_BODY_OUT(3)
    REAL(DP), INTENT(OUT), OPTIONAL :: ALPHA_RAW_DEG_OUT, Q_DYNAMIC_OUT

    REAL(DP) :: ANGLES(3), EULER_RATES(3), V_CG_TRANSLATING(3)
    REAL(DP) :: V_SAIL_POINT_BODY(3), V_REL_BODY(3), V_REL_H_BODY(3)
    REAL(DP) :: DELTA_NORMALIZED_DEG, C_CHORD_BODY(3)
    REAL(DP) :: E_UPSTREAM_BODY(3), ALPHA_RAW_DEG
    REAL(DP) :: E_DRAG_BODY(3), E_LIFT_BODY(3), Q_DYNAMIC
    LOGICAL :: DO_PRINT
    INTEGER :: STAGE

    V_CG_BODY = 0.0_DP
    OMEGA_BODY = 0.0_DP
    V_WIND_BODY = 0.0_DP
    V_REL_H_MAG = 0.0_DP
    ALPHA_DB_DEG = 0.0_DP
    CL = 0.0_DP
    CD = 0.0_DP
    FORCE_BODY = 0.0_DP
    MOMENT_BODY = 0.0_DP
    LOAD_6DOF = 0.0_DP
    IERR = SAIL_OK
    MESSAGE = ''
    IF (PRESENT(V_REL_BODY_OUT)) V_REL_BODY_OUT = 0.0_DP
    IF (PRESENT(E_UPSTREAM_BODY_OUT)) E_UPSTREAM_BODY_OUT = 0.0_DP
    IF (PRESENT(C_CHORD_BODY_OUT)) C_CHORD_BODY_OUT = 0.0_DP
    IF (PRESENT(ALPHA_RAW_DEG_OUT)) ALPHA_RAW_DEG_OUT = 0.0_DP
    IF (PRESENT(Q_DYNAMIC_OUT)) Q_DYNAMIC_OUT = 0.0_DP

    IF (NR < 6) THEN
      IERR = -1
      MESSAGE = 'ComputeSailPlanningDryRun requires NR >= 6.'
      RETURN
    END IF

    ANGLES = Y_STATE(4:6)
    EULER_RATES = Y_STATE(NR+4:NR+6)

    ! planning stores the mean forward speed separately from surge disturbance.
    V_CG_TRANSLATING(1) = U0 + Y_STATE(NR+1)
    V_CG_TRANSLATING(2) =      Y_STATE(NR+2)
    V_CG_TRANSLATING(3) =      Y_STATE(NR+3)

    ! VectorG2L uses planning's exact roll-pitch-yaw matrix and does not modify
    ! either input. Both translating/global vectors are expressed in the same
    ! axes before this single transformation to SailModule body axes.
    V_CG_BODY = VectorG2L(V_CG_TRANSLATING, ANGLES)
    V_WIND_BODY = VectorG2L(V_WIND_GLOBAL, ANGLES)

    CALL EulerRatesToBodyOmega(ANGLES, EULER_RATES, OMEGA_BODY)

    CALL ComputeSailModuleLoads( &
      V_WIND_BODY, V_CG_BODY, OMEGA_BODY, R_SAIL_BODY, DELTA_S_DEG, &
      V_SAIL_POINT_BODY, V_REL_BODY, V_REL_H_BODY, V_REL_H_MAG, &
      DELTA_NORMALIZED_DEG, C_CHORD_BODY, E_UPSTREAM_BODY, &
      ALPHA_RAW_DEG, ALPHA_DB_DEG, CL, CD, E_DRAG_BODY, &
      E_LIFT_BODY, Q_DYNAMIC, FORCE_BODY, MOMENT_BODY, LOAD_6DOF, &
      IERR, MESSAGE)

    ! Low wind is a nonfatal host condition. The SailModule API already clears
    ! every output; clear the load again here so no previous RK-stage load can
    ! ever survive at the host boundary.
    IF (IERR == SAIL_ERR_LOW_WIND_SPEED) LOAD_6DOF = 0.0_DP
    IF (PRESENT(V_REL_BODY_OUT)) V_REL_BODY_OUT = V_REL_BODY
    IF (PRESENT(E_UPSTREAM_BODY_OUT)) E_UPSTREAM_BODY_OUT = E_UPSTREAM_BODY
    IF (PRESENT(C_CHORD_BODY_OUT)) C_CHORD_BODY_OUT = C_CHORD_BODY
    IF (PRESENT(ALPHA_RAW_DEG_OUT)) ALPHA_RAW_DEG_OUT = ALPHA_RAW_DEG
    IF (PRESENT(Q_DYNAMIC_OUT)) Q_DYNAMIC_OUT = Q_DYNAMIC

    DO_PRINT = .FALSE.
    IF (PRESENT(PRINT_DIAGNOSTICS)) DO_PRINT = PRINT_DIAGNOSTICS
    STAGE = 0
    IF (PRESENT(RK_STAGE)) STAGE = RK_STAGE
    IF (DO_PRINT) THEN
      WRITE(*,'(A)') '--- Sail planning dry-run diagnostic ---'
      WRITE(*,'(A,ES24.16)') 'time = ', STAGE_TIME
      WRITE(*,'(A,I0)') 'RK stage = K', STAGE
      WRITE(*,'(A,3(ES24.16,1X))') 'roll pitch yaw = ', ANGLES
      WRITE(*,'(A,3(ES24.16,1X))') 'roll_dot pitch_dot yaw_dot = ', EULER_RATES
      WRITE(*,'(A,ES24.16)') 'U0 = ', U0
      WRITE(*,'(A,3(ES24.16,1X))') 'V_CG_TRANSLATING = ', V_CG_TRANSLATING
      WRITE(*,'(A,3(ES24.16,1X))') 'V_CG_BODY = ', V_CG_BODY
      WRITE(*,'(A,3(ES24.16,1X))') 'V_WIND_GLOBAL = ', V_WIND_GLOBAL
      WRITE(*,'(A,3(ES24.16,1X))') 'V_WIND_BODY = ', V_WIND_BODY
      WRITE(*,'(A,3(ES24.16,1X))') 'OMEGA_BODY = ', OMEGA_BODY
      WRITE(*,'(A,3(ES24.16,1X))') 'R_SAIL_BODY = ', R_SAIL_BODY
      WRITE(*,'(A,ES24.16)') 'DELTA_S_DEG = ', DELTA_S_DEG
      WRITE(*,'(A,ES24.16)') 'V_REL_H_MAG = ', V_REL_H_MAG
      WRITE(*,'(A,ES24.16)') 'ALPHA_DB_DEG = ', ALPHA_DB_DEG
      WRITE(*,'(A,ES24.16)') 'CL = ', CL
      WRITE(*,'(A,ES24.16)') 'CD = ', CD
      WRITE(*,'(A,3(ES24.16,1X))') 'FORCE_BODY = ', FORCE_BODY
      WRITE(*,'(A,3(ES24.16,1X))') 'MOMENT_BODY = ', MOMENT_BODY
      WRITE(*,'(A,6(ES24.16,1X))') 'LOAD_6DOF = ', LOAD_6DOF
      WRITE(*,'(A,I0)') 'IERR = ', IERR
      WRITE(*,'(A,A)') 'MESSAGE = ', TRIM(MESSAGE)
    END IF

    ! This adapter computes diagnostics only; host feedback is explicit in Timmotion.
  END SUBROUTINE ComputeSailPlanningDryRun


  SUBROUTINE EulerRatesToBodyOmega(ANGLES, EULER_RATES, OMEGA_BODY)
    REAL(DP), INTENT(IN) :: ANGLES(3)
    REAL(DP), INTENT(IN) :: EULER_RATES(3)
    REAL(DP), INTENT(OUT) :: OMEGA_BODY(3)
    REAL(DP) :: PITCH, YAW

    PITCH = ANGLES(2)
    YAW = ANGLES(3)

    ! Exact planning kinematics (the same 3-2-1 convention encoded by
    ! RotateG2L): [p q r]^T = E(pitch,yaw)[roll_dot pitch_dot yaw_dot]^T.
    ! p = cos(pitch)cos(yaw) roll_dot + sin(yaw) pitch_dot
    ! q =-cos(pitch)sin(yaw) roll_dot + cos(yaw) pitch_dot
    ! r = sin(pitch)          roll_dot + yaw_dot
    ! This forward map is finite at pitch=+/-pi/2; its inverse Euler-rate map
    ! is singular there. Roll does not appear in planning's E matrix.
    OMEGA_BODY(1) = COS(PITCH)*COS(YAW)*EULER_RATES(1) + &
                    SIN(YAW)*EULER_RATES(2)
    OMEGA_BODY(2) =-COS(PITCH)*SIN(YAW)*EULER_RATES(1) + &
                    COS(YAW)*EULER_RATES(2)
    OMEGA_BODY(3) = SIN(PITCH)*EULER_RATES(1) + EULER_RATES(3)
  END SUBROUTINE EulerRatesToBodyOmega

END MODULE SAILPLANNING_ADAPTER_MOD
