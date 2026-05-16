# What’s Modified

Summary of changes made to the Shrinaths Immersed Boundary project (MOVERS-N Euler / IBM solvers).

---

## 1. Folder organization

- **docs/** – Central documentation. `CODE_IMPROVEMENTS_SUMMARY.md` moved here from `Geometries/`.
- **geometry/** – Single location for all geometry definition files (`.txt`). Duplicates removed from `Geometries/` and `MPI_Parallel_Solver/`. Added `geometry/README.md`.
- **output/** – Archived results: `output/geometries/` (serial) and `output/mpi_solver/` (MPI `.dat`, `.vtk`, logs). Added `output/README.md`. New runs still write to the solver directory; move files here to archive.
- **Shaders/** – Metal shaders for flow visualization (see below).
- **Geometries/** – Now contains only serial solver source (no geometry `.txt`, no generated `.dat`).
- **MPI_Parallel_Solver/** – Contains only source, Makefile, scripts (no geometry `.txt`, no archived outputs).
- **.gitignore** – Added for build artifacts, optional solver outputs, editors, and Metal build products.

---

## 2. Metal shaders (Shaders/)

- **FlowVisualization.metal** – Vertex shader `flowVertex` (full-screen quad) and fragment shaders `flowFragment`, `flowFragmentVelocityMag` for 2D scalar flow fields. Uses scalar texture (0), mask texture (1), and `FlowUniforms` (buffer 0). Colormaps: jet, viridis, coolwarm, grayscale.
- **ContourOverlay.metal** – Fragment shader `contourFragment` for iso-contour overlay on the same data layout; use with `flowVertex` and alpha blending.
- **Shaders/README.md** – Data layout, `FlowUniforms` / `ContourUniforms`, entry points, and integration notes.

---

## 3. Code upgrades and performance

- **MathFunctions.h** (both `MPI_Parallel_Solver` and `Geometries`):
  - Fixed undefined `Max5`/`Min5` used by `Max6`/`Min6`.
  - Replaced custom min/max with `std::max` / `std::min` and `std::max({...})` / `std::min({...})`.
  - All small helpers marked `inline`; eigenvalue and `MoversN` use `std::abs` and simpler control flow.
- **Constants** – Replaced magic numbers `0.4` (γ−1) and `1.4` (γ) with `GAMMA_MINUS_1`, `GAMMA`, and `GAMMA_OVER_GAMMA_MINUS_1` in boundary conditions, primitive update, and boundary output (serial and MPI solver files).
- **Residual calculation** – Replaced `pow(diff, 2)` with `diff * diff` and `pow(Resn*, 0.5)` with `std::sqrt(Resn*)`; sqrt and inverse count computed once per step.
- **Config.h** – Added `USE_OPENMP` (default 0) for future hybrid MPI+OpenMP.
- **docs/UPGRADE_AND_PERFORMANCE.md** – Describes these upgrades and suggests further improvements (OpenMP, MPI overlap, I/O, GPU/Metal).

---

## 4. Build and Makefile

- **convert_to_vtk** – Added as a Makefile target; `make` now builds MPI solver, serial solver, locator, and VTK converter. `make help` documents it.
- **clean** – Removes `convert_to_vtk` as well as other executables.

---

## 5. Compiler warning fixes

- **Pragmas** – GCC-only pragmas (`#pragma GCC optimize`, `#pragma GCC target`) wrapped in `#if defined(__GNUC__) && !defined(__clang__)` to avoid Clang “unknown pragma” warnings.
- **LocatorFunctions.h** – Removed unused variables in `FootX` and `FootY` (`PtonLnSgmt`, `PtonLn`, `dist`, `qy` where unused).
- **LocatorV4.1.cpp** – Fixed “`/*` within block comment” (section header); removed unused variables `junk`, `Px`, `Py`, `Area`, `OutPt`, `X1`, `X2`, `Y1`, `Y2`.
- **MOVERSNEulerSolverMPI.cpp** – Marked unused globals with `__attribute__((unused))` (e.g. `F1`–`F4`, `G1L`–`G4R`, `a_max*`/`a_min*`, `pg`, `ug`, `rhog`, etc., `c1`).
- **MOVERSNEulerSolverIBMV4.2.cpp** (MPI_Parallel_Solver and Geometries) – Same for unused globals; fixed block-comment style (section dividers changed from `/** ... /***` to `/** ... */` and `/* ... */` borders) to remove “`/*` within block comment” warnings.
- **Geometries/LocatorFunctions.h** and **Geometries/MOVERSNEulerSolverIBMV4.2.cpp** – Same Locator and solver fixes as above for consistency.

Build with `make clean && make` in `MPI_Parallel_Solver` completes with zero warnings.

---

## 6. Documentation updates

- **README.md** – Directory structure, key folders, run commands (correct binary names), convert_to_vtk usage, and links to `docs/WhatsModified.md` and `docs/UPGRADE_AND_PERFORMANCE.md`.
- **MPI_Parallel_Solver/README.md** – Notes on geometry in `../geometry/` and archiving outputs to `../output/mpi_solver/`.
- **geometry/README.md** and **output/README.md** – Purpose and usage of each folder.

---

## 7. Running the solvers

- Solvers expect `GridPoints.dat` and `ImmersedBoundary.dat` in the current directory. If these were moved to `output/mpi_solver/`, copy them back (or symlink) before running, or run from a directory where they exist.
- Serial: `./MOVERSNEulerSolverIBM`
- MPI: `mpirun -np <N> ./MOVERSNEulerSolverMPI`
- VTK: `./convert_to_vtk <solution.dat> <output.vtk>`

---

## File change summary

| Area              | Files touched |
|-------------------|----------------|
| Layout            | New: `docs/`, `geometry/`, `output/`, `.gitignore`; moved/removed duplicates |
| Shaders           | New: `Shaders/FlowVisualization.metal`, `ContourOverlay.metal`, `Shaders/README.md` |
| Math/constants    | `MathFunctions.h`, `Physics.h` (new), `MOVERSNEulerSolverIBMV4.2.cpp` (both solver dirs), constants in MPI solver |
| Residual          | `MOVERSNEulerSolverIBMV4.2.cpp` (both), `Geometries/` copy |
| Build             | `Makefile`, `Config.h` |
| Warnings          | `LocatorFunctions.h`, `LocatorV4.1.cpp`, `MOVERSNEulerSolverMPI.cpp`, `MOVERSNEulerSolverIBMV4.2.cpp` (both dirs) |
| Docs              | `README.md`, `docs/UPGRADE_AND_PERFORMANCE.md`, `docs/WhatsModified.md`, `geometry/README.md`, `output/README.md`, `MPI_Parallel_Solver/README.md` |
| Navier–Stokes     | `SolverOptions.h`, `IbmBoundary.h`, `ViscousFlux.h`, `WENO5.h`, `HighOrderFlux.h`, both solver `.cpp` files |

---

## 8. Critical bug fixes (round 2 — applied from `docs/CodeReview.md`)

### Correctness bugs

- **MPI work partitioning** (`MOVERSNEulerSolverMPI.cpp`)
  Every kernel (`calculateInterfaceFluxes`, `calculateDissipativeFluxes`, `calculateTimeStep`, `updateConservativeVariables`, `updatePrimitiveVariables`, residual, `applyBoundaryConditions`, `isValidSolution`) now iterates only over the rank's owned `[local_first_i, local_last_i]` × `[2, Ny+1]` slab. Previous version had every rank computing the full domain — no real parallel speedup.
- **MPI ghost-cell exchange** rewritten to send `ghost_layers` columns of both conservative (`U1..U4`) and primitive (`Rho, u, v, P, H`) variables. Buffer sizing now reflects the 9 fields × 2 columns convention.
- **MPI immersed-wall reflection** changed from `u' = -u + 2 (u·n) n` (flipped tangential — wrong for inviscid slip wall) to `u' = u - 2 (u·n) n`, matching the serial solver.
- **Residual normalisation (RMS)** in both serial and MPI: `sqrt(Σr²)/N` → `sqrt(Σr² / N)`, which is the standard root-mean-square norm and is independent of grid size.
- **Locator polygon closure** (`LocatorV4.1.cpp`): `PointOnLine`, `Radius`, image-point, and ghost-cell-normal loops now traverse all `Nc` segments including the closing edge `(Nc, 1)` via `(k % Nc) + 1` indexing.
- **`HitOrMiss` rewritten** as a clean PNPOLY ray-cast (`LocatorFunctions.h`): half-open y-interval check + intersection x-coordinate compare. Replaces the original branch with two dead `return 0` paths.
- **Serial top/bottom farfield** copies from the *interior* row (`j±1`) instead of the *boundary* row's neighbour cell along i, eliminating the order-dependent smearing.
- **Mach number** uses velocity magnitude `|v| / a` (was `u / a`, ignoring `v`). Centralised in `updateAuxiliary()` so every site is consistent.
- **`MoversN` floating-point compare** uses `|Ul − Ur| > 1e-12` (was bare `Ul != Ur`), avoiding spurious branches when values are arithmetically equal but bit-different.
- **MPI grid read** index range fixed to match the serial pattern: `X[2..Nx+1][2..Ny+1]` (was `X[0..Nx+1][0..Ny+1]`, which left top-row metrics uninitialised and produced NaN once the masking `isnan` checks were removed).
- **MPI `Nx`/`Ny` broadcast after IBM read** so every rank sees the (smaller) cell-count dimensions, keeping ghost-exchange message lengths consistent.

### High-impact performance

- **Locator IMAGE-NEIGHBOUR-POINTS lookup** changed from O(N⁴) double-nested scan to O(1) `floor()` arithmetic per ghost cell. On a 790×1185 grid the locator now runs end-to-end in ≈10 s (previously minutes).
- **MPI `updateConservativeVariables`** drops 32 `std::isnan` checks per fluid cell per step in favour of cheap positivity floors on density and total energy. Validity is verified periodically by `isValidSolution()`.
- **OpenMP parallelism** added in the hot kernels (`calculateInterfaceFluxes`, `calculateDissipativeFluxes`, `calculateTimeStep`, `updateConservativeVariables`, `updatePrimitiveVariables`, residual reduction) of both serial and MPI solvers, gated on `_OPENMP`. Makefile auto-detects libomp on macOS (Homebrew) and `-fopenmp` on Linux. `make NO_OPENMP=1` opts out.
- **`endl → '\n'`** throughout the MPI solver (44 sites) and the serial solver's per-cell output paths, removing ~10⁶ stream flushes per output.
- **`updateAuxiliary` helper** (`Physics.h`) — replaces 7 lines of repeated E/a/H/T/Mach computation at 9+ call sites in each solver, ensuring consistency and eliminating one source of the Mach magnitude bug.
- **Output gather (`gatherFieldDouble`)** added in the MPI solver: each rank's owned slab is collected to rank 0 via `MPI_Gatherv` before writing `.dat`/`.vtk` files. Lets us gut the previous full-array `MPI_Bcast`-style approach in output paths.

### Verification

- Serial solver: 68 steps in 30 s (NO_OPENMP), 136 steps in 30 s with OpenMP (4 threads). Residuals decrease monotonically `1.20e-5 → 6.66e-6`.
- MPI solver, 1 rank: residual `1.99941e-06` at step 1, decreasing monotonically.
- MPI solver, 4 ranks: **identical** residual `1.99941e-06` at step 1 (deterministic across rank counts) — confirms partitioning is correct.
- Locator: 9.7 s end-to-end on 790×1185 (was multiple minutes previously).
- Build clean with `-Wall` on Apple Clang 20 and Homebrew GCC, both `make` and `make NO_OPENMP=1`.
- **Full-case MPI run** (790×1185, 10 000 steps, 4 ranks): completed in ~21.5 min wall time; residuals bounded, no NaNs; 20 solution/VTK snapshots written.

---

## 9. Navier–Stokes, IBM no-slip wall, and optional WENO-5 (round 3)

### Default physics (`SolverOptions.h`)

- **`USE_NAVIER_STOKES` default is `1`** — laminar Navier–Stokes is on without extra compile flags. Set to `0` for legacy **inviscid Euler + slip** IBM only.
- Transport: `LAMINAR_VISCOSITY`, `PRANDTL_NUMBER`, `VISCOUS_CFL_FRACTION`.

### Viscous flux (`ViscousFlux.h`)

- Newtonian stress τ and Fourier heat flux \(q = -k\nabla T\), \(k = \mu c_p / \mathrm{Pr}\).
- Face fluxes `VF1R`…`VF4B` added to the conservative update (with `IF` + `DF`).
- Gradients of \(u, v, T\) on cell centres via metric Jacobian \((C_x, C_y)\).
- Explicit viscous timestep cap: \(\Delta t \lesssim C \,\rho h^2 / \mu\) combined with convective CFL via `MPI_Allreduce` / local minimum.

### IBM no-slip + thermal wall (`IbmBoundary.h`)

- **Velocity (when `USE_NAVIER_STOKES`)**: after image **bilinear** interpolation, **no-slip** stationary wall: \(\mathbf{u}_\mathrm{ghost} = -\mathbf{u}_\mathrm{image}\).
- **Thermal** (choose in `SolverOptions.h`):
  - **Adiabatic (default)**: `IBM_WALL_ADIABATIC 1` → \(T_\mathrm{ghost} = T_\mathrm{image}\) (\(\partial T/\partial n = 0\)).
  - **Fixed \(T_w\)**: `IBM_WALL_FIXED_TW 1`, `WALL_TEMPERATURE_TW` → \(T_\mathrm{ghost} = 2 T_w - T_\mathrm{image}\); then \(P = \rho R T\) at the ghost cell.
- **Euler fallback** (`USE_NAVIER_STOKES 0`): inviscid **slip** wall (`ibmSlipVelocity`) unchanged.
- **Viscous faces**: `viscousFaceActive` enables flux on **fluid–fluid** and **fluid–IBM** faces (not only interior `Fl==1` pairs).
- **MPI**: IBM BC uses full bilinear (aligned with serial); ghost exchange packs **10 fields** (adds `T`) when Navier–Stokes is on.

### Optional WENO-5 (`WENO5.h`, `HighOrderFlux.h`)

- `USE_WENO5 0` by default — component-wise Jiang–Shu WENO-5 on \((\rho, u, v, P)\) at interior faces; boundary bands fall back to central flux.
- `eulerFluxDotN` evaluates \(F \cdot \hat{n}\) from reconstructed left/right states.

### New / updated files

| File | Role |
|------|------|
| `SolverOptions.h` | `USE_NAVIER_STOKES`, `USE_WENO5`, μ, Pr, IBM thermal BC |
| `IbmBoundary.h` | No-slip / slip velocity; wall temperature; `viscousFaceActive` |
| `ViscousFlux.h` | τ, heat flux, `viscousDtScale` |
| `WENO5.h` | `weno5_q_minus`, `weno5_q_plus` |
| `HighOrderFlux.h` | `eulerFluxDotN` |
| `MOVERSNEulerSolverMPI.cpp` | Viscous kernels, IBM BC, WENO interface flux, ghost `T` |
| `MOVERSNEulerSolverIBMV4.2.cpp` | Same physics in serial main loop |

### Startup message (MPI rank 0)

Example: `Physics: Navier-Stokes (mu=0.0001), IBM no-slip wall, adiabatic wall`
