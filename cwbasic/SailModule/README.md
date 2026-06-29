# SailModule

This directory contains an independent Fortran90 sail aerodynamic load framework.
It is not connected to the WALCS main program yet.

## Files

- `SAILPARAM_MOD.F90`: shared precision, constants, and air density.
- `SAILDATABASE_MOD.F90`: database storage and `ReadSailDatabase()`.
- `SAILINTERP_MOD.F90`: linear interpolation through `GetSailCoeff()`.
- `SAILFORCE_MOD.F90`: six-degree-of-freedom external force through `GetSailForce()`.
- `sail_database.dat`: sample `alpha_deg, CL, CD` table.

## Call Order

```text
ReadSailDatabase("sail_database.dat", ierr)
GetSailCoeff(alphaDeg, CL, CD, ierr)
GetSailForce(CL, CD, apparentWindSpeed, apparentWindAngleDeg,
             sailArea, xCE, yCE, zCE, Fsail, ierr)
```

`Fsail(1,1:6)` is the real part of the external force and moment vector.
`Fsail(2,1:6)` is reserved for the imaginary part and is currently zero.

The intended future WALCS integration point is immediately before `shipmotion`
is called in `MAINCAL_MOD.F90`:

```text
F(:,1:6,9) = F(:,1:6,9) + Fsail(:,1:6)
```

## Coordinate Convention

- X forward
- Y to port
- Z upward
- `apparentWindAngleDeg` is measured in the XY plane from +X toward +Y
- Drag acts along the apparent wind velocity direction
- Lift is 90 degrees counter-clockwise from drag in the XY plane
