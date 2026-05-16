#ifndef CONFIG_H
#define CONFIG_H

#include "SolverOptions.h"

/*****************************************************************************************/
/** Default: USE_NAVIER_STOKES=1 → laminar viscosity + IBM no-slip + adiabatic wall.      */
/**   USE_WALL_FUNCTION=1 — log-law IBM wall (low-speed; see cases/viscous_channel).       */
/**   IBM_WALL_FIXED_TW / WALL_TEMPERATURE_TW — isothermal wall (see SolverOptions.h).   */
/**   USE_WENO5=1 — WENO-5 inviscid interface (optional).                                */
/** Set USE_NAVIER_STOKES=0 for legacy Euler + slip IBM.                                 */
/*****************************************************************************************/

/*****************************************************************************************/
/** CONFIGURATION PARAMETERS FOR MPI CFD SOLVER */
/*****************************************************************************************/

// Grid and Memory Configuration
#define MAX_GRID_X 1500
#define MAX_GRID_Y 1500
#define MAX_GHOST_LAYERS 2

// Physical Constants
#define GAMMA 1.4
#define GAS_CONSTANT 1.0
#define PI 3.141592653589793

// Derived Constants
#define GAMMA_MINUS_1 (GAMMA - 1.0)
#define GAMMA_OVER_GAMMA_MINUS_1 (GAMMA / (GAMMA - 1.0))
#define INV_GAMMA_MINUS_1 (1.0 / (GAMMA - 1.0))

// Solver Parameters
#define CFL_NUMBER 0.05
#define MAX_ITERATIONS 10000
#define OUTPUT_FREQUENCY 500
#define VALIDATION_FREQUENCY 100

// Flow Conditions (Default Values)
#define MACH_INFINITY 2.0
#define DENSITY_INFINITY 1.4
#define PRESSURE_INFINITY 1.0

// Convergence Parameters
#define CONVERGENCE_TOLERANCE 1e-6
#define MIN_RESIDUAL 1e-12
#define MAX_RESIDUAL 1e12

// MPI Configuration
#define DEFAULT_GHOST_LAYERS 2
#define MAX_MPI_PROCESSES 1024

// File Names
#define GRID_FILE "GridPoints.dat"
#define IBM_FILE "ImmersedBoundary.dat"
#define RESIDUAL_FILE "Residues_MPI.dat"
#define OUTPUT_PREFIX "Aerospike_IBM_MPI_Solution"

// Debug and Validation Flags
#ifdef DEBUG
#define DEBUG_OUTPUT 1
#define VALIDATE_EVERY_STEP 1
#else
#define DEBUG_OUTPUT 0
#define VALIDATE_EVERY_STEP 0
#endif

// Performance Optimization Flags
#define USE_COMPILER_OPTIMIZATIONS 1
#define CACHE_OPTIMIZE 1
#define VECTORIZE_LOOPS 1
// USE_OPENMP is detected automatically via the _OPENMP macro defined when -fopenmp is in use.
// If you want to *force* serial loops even when compiling with -fopenmp, you can put
// `#undef _OPENMP` here, but normally just `make NO_OPENMP=1` is what you want.

// Domain Decomposition Options
#define DECOMP_1D_I 1 // 1D decomposition in i-direction
#define DECOMP_1D_J 0 // 1D decomposition in j-direction (not implemented)
#define DECOMP_2D 0   // 2D decomposition (future enhancement)

// Output Control
#define OUTPUT_SOLUTION_FILES 1
#define OUTPUT_BOUNDARY_SOLUTION 1
#define OUTPUT_RESIDUAL_DETAILS 1

// Memory Management
#define USE_DYNAMIC_ALLOCATION 1
#define PREALLOCATE_BUFFERS 1

#endif // CONFIG_H
