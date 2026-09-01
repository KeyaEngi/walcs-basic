!===============================================================================
! SAILPARAM_MOD
!
! Central parameter definitions for the independent sail aerodynamics module.
!
! Body-fixed coordinate system:
!   +x_b points toward the bow.
!   +y_b points toward port.
!   +z_b points upward.
!
! Six-degree-of-freedom load order: (Fx, Fy, Fz, Mx, My, Mz).
! The current aerodynamic model considers horizontal relative wind only.
! The sail chord line is an undirected axis with a 180-degree period.
! This module has no dependency on WALCS, WALCS-LE, planning, or any other
! host solver program.
!===============================================================================
MODULE SAILPARAM_MOD
  IMPLICIT NONE
  PRIVATE

  ! Numerical precision
  INTEGER, PARAMETER, PUBLIC :: DP = SELECTED_REAL_KIND(15, 307)

  ! Mathematical constants
  REAL(DP), PARAMETER, PUBLIC :: PI         = ACOS(-1.0_DP)
  REAL(DP), PARAMETER, PUBLIC :: DEG_TO_RAD = PI / 180.0_DP
  REAL(DP), PARAMETER, PUBLIC :: RAD_TO_DEG = 180.0_DP / PI

  ! Air and sail properties
  REAL(DP), PARAMETER, PUBLIC :: RHO_AIR    = 1.225_DP                  ! kg/m^3
  REAL(DP), PARAMETER, PUBLIC :: SAIL_CHORD = 2.5_DP                    ! m
  REAL(DP), PARAMETER, PUBLIC :: SAIL_SPAN  = 6.0_DP                    ! m
  REAL(DP), PARAMETER, PUBLIC :: SAIL_AREA  = SAIL_CHORD * SAIL_SPAN   ! m^2

  ! Numerical tolerances
  ! Horizontal relative wind speed below VREL_MIN is treated as low wind by
  ! downstream modules to avoid division by zero; no check is performed here.
  REAL(DP), PARAMETER, PUBLIC :: VREL_MIN     = 1.0E-8_DP    ! m/s
  REAL(DP), PARAMETER, PUBLIC :: ANGLE_TOL_DEG = 1.0E-10_DP  ! deg
  ! Numerical tolerance used in general database comparisons.
  REAL(DP), PARAMETER, PUBLIC :: DATABASE_TOL  = 1.0E-12_DP

  ! Coordinate indices for three-dimensional vectors
  INTEGER, PARAMETER, PUBLIC :: IDX_X = 1
  INTEGER, PARAMETER, PUBLIC :: IDX_Y = 2
  INTEGER, PARAMETER, PUBLIC :: IDX_Z = 3

  ! Six-degree-of-freedom indices
  ! Body axes: +x_b bow, +y_b port, +z_b upward.
  ! Load order: (Fx, Fy, Fz, Mx, My, Mz).
  INTEGER, PARAMETER, PUBLIC :: N_DOF     = 6
  INTEGER, PARAMETER, PUBLIC :: IDX_SURGE = 1
  INTEGER, PARAMETER, PUBLIC :: IDX_SWAY  = 2
  INTEGER, PARAMETER, PUBLIC :: IDX_HEAVE = 3
  INTEGER, PARAMETER, PUBLIC :: IDX_ROLL  = 4
  INTEGER, PARAMETER, PUBLIC :: IDX_PITCH = 5
  INTEGER, PARAMETER, PUBLIC :: IDX_YAW   = 6

  ! Force and moment index aliases
  INTEGER, PARAMETER, PUBLIC :: IDX_FX = 1
  INTEGER, PARAMETER, PUBLIC :: IDX_FY = 2
  INTEGER, PARAMETER, PUBLIC :: IDX_FZ = 3
  INTEGER, PARAMETER, PUBLIC :: IDX_MX = 4
  INTEGER, PARAMETER, PUBLIC :: IDX_MY = 5
  INTEGER, PARAMETER, PUBLIC :: IDX_MZ = 6

  ! Error codes
  INTEGER, PARAMETER, PUBLIC :: SAIL_OK                           = 0
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_FILE_NOT_FOUND           = 1
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_FILE_OPEN                = 2
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_FILE_READ                = 3
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_DATABASE_EMPTY           = 4
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_DATABASE_ORDER           = 5
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_DATABASE_DUPLICATE       = 6
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_DATABASE_RANGE           = 7
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_DATABASE_NOT_INITIALIZED = 8
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_INVALID_INPUT            = 9
  ! Recognizable low-wind state; the total API may return zero loads for it.
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_LOW_WIND_SPEED            = 10
  INTEGER, PARAMETER, PUBLIC :: SAIL_ERR_INTERPOLATION             = 11

  ! Database angle limits
  ! The database covers 0 to 180 degrees. Sail chord orientations separated
  ! by 180 degrees are equivalent because the chord is an undirected axis.
  REAL(DP), PARAMETER, PUBLIC :: ALPHA_DB_MIN_DEG       = 0.0_DP    ! deg
  REAL(DP), PARAMETER, PUBLIC :: ALPHA_DB_MAX_DEG       = 180.0_DP  ! deg
  REAL(DP), PARAMETER, PUBLIC :: SAIL_ANGLE_PERIOD_DEG  = 180.0_DP  ! deg

  ! Default file names
  CHARACTER(LEN=17), PARAMETER, PUBLIC :: DEFAULT_DATABASE_FILE = &
    'sail_database.dat'

END MODULE SAILPARAM_MOD
