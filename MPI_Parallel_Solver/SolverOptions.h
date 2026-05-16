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

/*****************************************************************************************/
/** Wall function on IBM (with Navier–Stokes): match wall shear via log-law instead of     */
/** resolving the viscous sublayer on the grid. Set to 0 for resolved no-slip (u=-u_im). */
/*****************************************************************************************/
#ifndef USE_WALL_FUNCTION
#define USE_WALL_FUNCTION 0
#endif

#ifndef WALL_KAPPA
#define WALL_KAPPA 0.41
#endif

#ifndef WALL_E
#define WALL_E 9.793
#endif

#ifndef WALL_FUNCTION_NEWTON_ITERS
#define WALL_FUNCTION_NEWTON_ITERS 25
#endif

/** Optional override for channel / low-speed wall-function tests (0 = use Config.h). */
#ifndef TEST_CHANNEL_MACH
#define TEST_CHANNEL_MACH 0.0
#endif

#endif // SOLVEROPTIONS_H
