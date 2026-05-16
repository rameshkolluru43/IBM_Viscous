#ifndef SOLVEROPTIONS_H
#define SOLVEROPTIONS_H

/*****************************************************************************************/
/** Physics / discretization toggles (included before solvers and from Config.h).          */
/*****************************************************************************************/

/* Navier–Stokes (laminar viscosity + IBM no-slip wall) is ON by default. Set to 0 for     */
/* legacy inviscid Euler with slip IBM only.                                               */
#ifndef USE_NAVIER_STOKES
#define USE_NAVIER_STOKES 1
#endif

#ifndef USE_WENO5
#define USE_WENO5 0
#endif

#ifndef LAMINAR_VISCOSITY
#define LAMINAR_VISCOSITY 1.0e-4
#endif

#ifndef PRANDTL_NUMBER
#define PRANDTL_NUMBER 0.72
#endif

#ifndef VISCOUS_CFL_FRACTION
#define VISCOUS_CFL_FRACTION 0.2
#endif

#ifndef WENO_EPSILON
#define WENO_EPSILON 1.0e-6
#endif

/*****************************************************************************************/
/** IBM wall thermal BC (only used when USE_NAVIER_STOKES == 1).                          */
/* Set IBM_WALL_FIXED_TW to 1 and define WALL_TEMPERATURE_TW for isothermal wall.         */
/* Otherwise IBM_WALL_ADIABATIC (default): ∂T/∂n = 0 at the wall.                         */
/*****************************************************************************************/
#ifndef IBM_WALL_FIXED_TW
#define IBM_WALL_FIXED_TW 0
#endif

#ifndef IBM_WALL_ADIABATIC
#define IBM_WALL_ADIABATIC 1
#endif

#ifndef WALL_TEMPERATURE_TW
#define WALL_TEMPERATURE_TW 1.0
#endif

#if IBM_WALL_FIXED_TW && IBM_WALL_ADIABATIC
#undef IBM_WALL_ADIABATIC
#define IBM_WALL_ADIABATIC 0
#endif

#endif // SOLVEROPTIONS_H
