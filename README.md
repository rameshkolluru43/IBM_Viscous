# Shrinaths Immersed Boundary Project

## Overview
This project implements numerical solvers for **2D compressible flow** using the **Immersed Boundary Method (IBM)** on structured grids. The default build solves the **laminar Navier–Stokes equations** (MOVERS-N inviscid flux + Newtonian viscous terms) with **no-slip, adiabatic IBM walls**. An **inviscid Euler + slip-wall** mode is available for legacy runs. **MPI** domain decomposition and optional **OpenMP** / **WENO-5** are supported.

## Directory Structure

```
Shrinaths_Immersed_Boundary/
├── README.md                 # This file
├── docs/                     # Project documentation
│   ├── CodeReview.md
│   ├── CODE_IMPROVEMENTS_SUMMARY.md
│   ├── UPGRADE_AND_PERFORMANCE.md
│   └── WhatsModified.md      # Changelog of applied changes
├── geometry/                 # Geometry definition files (.txt)
│   ├── README.md
│   ├── BluntBody.txt, BluntBody_withAeroSpike.txt, DBump.txt, ...
│   └── TriangularPrism.txt, Domain_TriangularPrism.txt
├── output/                   # Archived simulation outputs
│   ├── README.md
│   ├── geometries/           # Serial solver outputs
│   └── mpi_solver/           # MPI solution .dat, .vtk, logs
├── Shaders/                  # Metal shaders for flow visualization
│   ├── README.md
│   ├── FlowVisualization.metal
│   └── ContourOverlay.metal
├── Geometries/               # Serial IBM solver source
│   ├── MOVERSNEulerSolverIBMV4.2.cpp, LocatorV4.1.cpp
│   ├── LocatorFunctions.h, MathFunctions.h, SolverFunctions.h
│   └── README.md (see MPI_Parallel_Solver for build/run)
└── MPI_Parallel_Solver/      # Parallel solver, utilities, build
    ├── Makefile, test_mpi.sh, README.md, SolverOptions.h
    ├── MOVERSNEulerSolverMPI.cpp, MOVERSNEulerSolverIBMV4.2.cpp
    ├── convert_to_vtk.cpp, LocatorV4.1.cpp, Config.h
    ├── Physics.h, IbmBoundary.h, ViscousFlux.h, WENO5.h, HighOrderFlux.h
    └── LocatorFunctions.h, MathFunctions.h, SolverFunctions.h
```

### Key folders

- **geometry/** – Single place for geometry definitions (`.txt`). Use `../geometry/<file>.txt` when running from a solver directory.
- **output/** – Archived results (`.dat`, `.vtk`, `.log`). New runs still write to the solver directory; move files here to keep it tidy.
- **docs/** – Documentation (e.g. code improvements summary).
- **Shaders/** – Metal shaders for GPU visualization of flow fields.
- **Geometries/** – Serial solver source only (no geometry files; those live in `geometry/`).
- **MPI_Parallel_Solver/** – MPI solver, locator, VTK converter, and build (Makefile, scripts).

## Key Features
- **Immersed Boundary Method (IBM)**: Complex geometries (aerospike, blunt body, bumps, etc.) on a Cartesian grid without body-fitted meshes.
- **Navier–Stokes (default)**: Laminar viscosity, viscous flux on fluid and fluid–wall faces, explicit viscous CFL limit.
- **IBM wall model**: No-slip velocity (`u_ghost = −u_image` after bilinear interpolation); **adiabatic** or **fixed T_w** thermal BC (`SolverOptions.h`).
- **MOVERS-N scheme**: Central inviscid interface flux plus MOVERS-N dissipation; optional **WENO-5** for higher-order inviscid reconstruction.
- **MPI + OpenMP**: 1D domain decomposition in *i*; shared-memory parallelism in hot loops when built with OpenMP.
- **VTK output**: `.dat` and `.vtk` snapshots for ParaView / VisIt.
- **Legacy Euler mode**: `USE_NAVIER_STOKES 0` → inviscid flow with **slip** IBM (no viscous flux).

## Physics configuration

Edit `MPI_Parallel_Solver/SolverOptions.h` (included by `Config.h`):

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `USE_NAVIER_STOKES` | `1` | Navier–Stokes + no-slip IBM; set `0` for Euler + slip |
| `LAMINAR_VISCOSITY` | `1e-4` | Dynamic viscosity μ |
| `PRANDTL_NUMBER` | `0.72` | Prandtl number |
| `IBM_WALL_ADIABATIC` | `1` | Adiabatic wall (∂T/∂n = 0) |
| `IBM_WALL_FIXED_TW` | `0` | Set to `1` for isothermal wall at `WALL_TEMPERATURE_TW` |
| `USE_WENO5` | `0` | Jiang–Shu WENO-5 on (ρ, u, v, P) at interior faces |

```sh
# Legacy inviscid Euler + slip IBM only
make NO_OPENMP=1 CXXFLAGS='-O3 -Wall -std=c++11 -march=native -mtune=native -DUSE_NAVIER_STOKES=0' \
  MPIFLAGS='-O3 -Wall -std=c++11 -march=native -mtune=native -DUSE_NAVIER_STOKES=0'
```

## Getting Started

### Prerequisites
- C++ compiler (e.g., g++, clang++)
- MPI library (e.g., OpenMPI, MPICH)
- VTK tools for visualization (optional)

### Build Instructions
1. Go to `MPI_Parallel_Solver`.
2. Run `make` to build the MPI solver, serial solver, locator, and VTK converter.

```sh
cd MPI_Parallel_Solver
make
```

Targets: `MOVERSNEulerSolverMPI`, `MOVERSNEulerSolverIBM`, `LocatorV4.1`, `convert_to_vtk`. Use `make help` for the full list.

### Running Simulations

Solver executables and input files live in `MPI_Parallel_Solver/`. The solvers expect `GridPoints.dat` and `ImmersedBoundary.dat` in the current directory (copy from `output/mpi_solver/` if you archived them there).

- **Serial solver** (single-core):

```sh
cd MPI_Parallel_Solver
./MOVERSNEulerSolverIBM
```

- **Parallel solver** (MPI):

```sh
cd MPI_Parallel_Solver
mpirun -np <num_procs> ./MOVERSNEulerSolverMPI
```

- **Convert solution to VTK**:

```sh
cd MPI_Parallel_Solver
./convert_to_vtk <solution.dat> <output.vtk>
```

### Example Test Cases
- `BluntBody.txt`, `BluntBody_withAeroSpike.txt`, `DBump.txt`, `TriangularPrism.txt` are provided for benchmarking and validation.

## File Descriptions
- **MOVERSNEulerSolverIBMV4.2.cpp** / **MOVERSNEulerSolverMPI.cpp**: Serial and MPI solvers (shared physics).
- **SolverOptions.h**: Physics toggles (NS / Euler, WENO, μ, Pr, wall T_w).
- **IbmBoundary.h**: No-slip / slip IBM velocity; adiabatic / fixed T_w; viscous face mask.
- **ViscousFlux.h**, **WENO5.h**, **HighOrderFlux.h**, **Physics.h**: Viscous flux, WENO-5, Euler flux helper, auxiliary variables.
- **LocatorV4.1.cpp**: Grid and IBM preprocessing.
- **convert_to_vtk.cpp**: `.dat` → VTK conversion.
- **Makefile**, **test_mpi.sh**: Build and MPI smoke tests.

## Documentation
- **`docs/CodeReview.md`** – Detailed code review with prioritized bugs and modification suggestions (correctness, performance, structure).
- **`docs/WhatsModified.md`** – Summary of what was changed (folder layout, shaders, upgrades, warnings, build, run).
- **`docs/CODE_IMPROVEMENTS_SUMMARY.md`** – Earlier code improvements and refactors.
- **`docs/UPGRADE_AND_PERFORMANCE.md`** – Upgrades (MathFunctions, constants, residual) and further performance ideas (OpenMP, MPI overlap, I/O, GPU).
- **`MPI_Parallel_Solver/README.md`** – Parallel solver usage and configuration.
- **`geometry/README.md`** and **`output/README.md`** – Data layout and where to put geometry and results.

## Contributing
1. Fork the repository.
2. Create a feature branch.
3. Commit your changes.
4. Submit a pull request.

## License
This project is released under the MIT License.

## Contact
For questions or collaboration, please contact the project maintainer.
