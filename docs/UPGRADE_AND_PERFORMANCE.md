# Immersed Boundary Code: Upgrades and Performance

This document summarizes upgrades applied to the MOVERS-N / Immersed Boundary solvers and gives guidance for further performance improvements. The **default build** is **laminar Navier–Stokes** with **IBM no-slip** walls; see `SolverOptions.h` and `docs/WhatsModified.md` §9.

---

## 1. Upgrades Applied

### 1.1 MathFunctions.h

- **Fixed missing Max5/Min5**: `Max6`/`Min6` referred to undefined `Max5`/`Min5`; all are now implemented via `std::max`/`std::min`.
- **Standard library use**: Replaced hand-written min/max with `std::max`, `std::min`, and `std::max({...})`/`std::min({...})` for clarity and compiler optimization.
- **Inline and const correctness**: All small template helpers are `inline`; eigenvalue helpers use local variables and `std::abs` for better inlining.
- **MoversN**: Simplified branch structure for clearer codegen.

### 1.2 Constants and Maintainability

- **No magic numbers**: Replaced literal `0.4` (γ−1) and `1.4` (γ) with `GAMMA_MINUS_1`, `GAMMA`, and `GAMMA_OVER_GAMMA_MINUS_1` in:
  - Boundary conditions (inlet, farfield, outflow, ghost)
  - Primitive update (P, a, H from U)
  - Boundary output
- **Single source of truth**: Physical constants live in the top-level constants or `Config.h` (MPI build).

### 1.3 Residual Calculation

- **Faster residual loop**: Replaced `pow(diff, 2)` with `diff * diff` and `pow(Resn*, 0.5)` with `std::sqrt(Resn*)`; compute sqrt and inverse count once per step.
- **Cleaner accumulation**: Use `+=` and local `d1..d4` for readability and cache-friendly access.

### 1.4 Navier–Stokes and IBM wall (`SolverOptions.h`, `IbmBoundary.h`, `ViscousFlux.h`)

- **Default `USE_NAVIER_STOKES 1`**: Laminar Newtonian viscous flux on interior and fluid–IBM faces; explicit viscous CFL cap combined with convective limit.
- **IBM no-slip** (automatic with NS): \(\mathbf{u}_\mathrm{ghost} = -\mathbf{u}_\mathrm{image}\) after bilinear interpolation at `Gst==1` cells.
- **Wall thermal**: adiabatic (`T_ghost = T_image`) or fixed `T_w` (`T_ghost = 2 T_w - T_image`); configure via `IBM_WALL_ADIABATIC` / `IBM_WALL_FIXED_TW` and `WALL_TEMPERATURE_TW`.
- **Euler legacy**: `USE_NAVIER_STOKES 0` restores inviscid slip IBM and disables viscous flux.

### 1.5 Optional WENO-5 (`USE_WENO5`, `WENO5.h`)

- Component-wise fifth-order reconstruction on \((\rho, u, v, P)\) for inviscid interface flux; MOVERS-N dissipation unchanged. Off by default.

### 1.6 MPI and OpenMP (round 2)

- True 1D work partitioning, 9- or 10-field ghost exchange (adds `T` when NS), `MPI_Gatherv` output gather, OpenMP in hot loops, locator O(1) image lookup.

---

## 2. Performance Tips (Already Reflected in Code)

- **Loop order**: Inner loop over `i` (row index) matches row-major layout and improves cache use.
- **Local caching**: Interface flux blocks cache `Rho[i][j]`, `u[i][j]`, etc., in locals to cut repeated loads.
- **Compiler flags**: `-O3 -march=native -mtune=native` (and `-ffast-math` if acceptable) improve throughput; see `Makefile` and `#pragma GCC optimize`.

---

## 3. Further Performance Improvements (Recommended)

### 3.1 Build and Run

- **Add `convert_to_vtk` to Makefile**: So VTK conversion is built with the same flags as the solver.
- **Optional OpenMP**: For the serial solver, add `#pragma omp parallel for` on the main interior loops (flux, update, residual) and link with `-fopenmp`. Use `collapse(2)` only if the loop body is cheap and load-balanced.
- **MPI + OpenMP hybrid**: Use OpenMP within each MPI rank for the local domain; reduce MPI size and increase threads per node to lower communication and improve scaling on many-core nodes.

### 3.2 Memory and Cache

- **Dynamic allocation**: Replace fixed `double A[MAX_GRID_X][MAX_GRID_Y]` with 1D `std::vector<double>` or aligned buffers of size `(Nx+pad)*(Ny+pad)` and index as `A[i*(Ny+pad)+j]`. Reduces footprint for small runs and keeps layout predictable.
- **SoA for conserved variables**: Store `U1,U2,U3,U4` as four separate contiguous arrays (Structure of Arrays) to improve vectorization and cache use in flux and update kernels.
- **Loop tiling**: For very large grids, tile the i–j loops (e.g. 64×64 or 128×128) to improve L2/L3 reuse in flux and primitive update.

### 3.3 MPI

- **Overlap communication and computation**: After boundary conditions, post `MPI_Isend`/`MPI_Irecv` for ghost layers, then compute fluxes in the interior; `MPI_Waitall` before updating ghost cells. This hides part of the communication cost.
- **Derived datatypes**: Use `MPI_Type_create_subarray` for 2D ghost regions to avoid manual pack/unpack and reduce copies.
- **Residual**: You already use `MPI_Allreduce` for global residuals; keep residual and other global reductions to a minimum per step.

### 3.4 I/O

- **Less frequent residual file write**: Write `Residues.dat` every N steps (e.g. 10–100) instead of every step to reduce I/O.
- **Binary solution output**: Optionally write solution `.dat` in raw binary (e.g. `fwrite`) and document the format; convert to VTK in a post step. Speeds up large runs and reduces disk use.
- **Async I/O**: Offload solution write to a separate thread or use non-blocking I/O so the main loop is not blocked.

### 3.5 Algorithmic (Future)

- **Multi-stage time integration**: Replace single-stage Euler with a low-storage Runge–Kutta (e.g. RK2–SSP or RK3–SSP) for better stability and possibly larger CFL.
- **Local time stepping**: You already have local `dt[i][j]`; using it for update (instead of global `dtg`) can speed up steady-state convergence in regions that allow a larger step.
- **GPU/Metal compute**: The costliest kernels (interface flux, dissipative flux, primitive update) are good candidates for Metal compute shaders or CUDA/HIP; keep grid and IBM data on the host and transfer only conserved/primitive arrays and metrics.

---

## 4. Validation After Changes

- **Regression**: Run a fixed grid (e.g. 120×42) for a few hundred steps and compare residuals and a few solution values (e.g. at one cell) before/after changes.
- **Conservation**: Integrate total mass, momentum, and energy over the fluid domain and check drift over time (especially after changing flux or BCs).
- **MPI vs serial**: For the same grid and time steps, compare serial and MPI (e.g. 2–4 ranks) solutions; norms of the difference should be at the level of round-off.

---

## 5. File and Config Summary

| Item | Location |
|------|----------|
| Physics toggles | `SolverOptions.h` (`USE_NAVIER_STOKES`, `USE_WENO5`, μ, Pr, wall T_w) |
| IBM wall BC | `IbmBoundary.h` |
| Viscous flux | `ViscousFlux.h` |
| WENO / Euler flux | `WENO5.h`, `HighOrderFlux.h` |
| Auxiliary variables | `Physics.h` (`updateAuxiliary`) |
| Constants (MPI) | `Config.h` (includes `SolverOptions.h`) |
| Constants (serial) | `MOVERSNEulerSolverIBMV4.2.cpp` top |
| Math helpers | `MathFunctions.h` |
| Solvers | `MOVERSNEulerSolverMPI.cpp`, `MOVERSNEulerSolverIBMV4.2.cpp` |
| MPI ghost exchange | `exchangeGhostCells()` (9 or 10 fields) |
| VTK conversion | `convert_to_vtk.cpp` (Makefile target) |

### Navier–Stokes tuning notes

- **Viscous runs** need smaller \(\Delta t\) when \(\mu\) is large; `VISCOUS_CFL_FRACTION` limits the step. Reduce `LAMINAR_VISCOSITY` or CFL if the case is stiff.
- **Re** is not a single input — set \(\mu\) relative to grid spacing and freestream scales for your geometry.
- For **isothermal walls**, set `IBM_WALL_FIXED_TW 1` and `WALL_TEMPERATURE_TW`; verify `P = \rho R T` at ghost cells is consistent with your freestream nondimensionalisation.

Applying the items in §3 incrementally and re-running your usual tests will help keep the code stable while improving performance and maintainability.
