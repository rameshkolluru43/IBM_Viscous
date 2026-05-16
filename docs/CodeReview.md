# Detailed Code Review and Suggested Modifications

This is a thorough scan of the MOVERS-N Euler / Immersed Boundary codebase (`MPI_Parallel_Solver/` and `Geometries/` mirrors). Issues are grouped by severity. Each item includes file:line references and a concrete recommendation.

Severity legend:

- **[BUG]** likely incorrect behavior or wrong result
- **[HOT]** correctness-OK but large performance/scaling impact
- **[STRUCT]** code quality, maintainability, robustness
- **[NIT]** small cleanup
- **[DONE]** applied — see `docs/WhatsModified.md` for the implementation

---

## Status (2026-05): items applied in the latest passes

Critical correctness bugs (1.1 – 1.8), high-impact performance items (2.1 – 2.4), and Navier–Stokes / IBM wall extensions (section 3) are implemented. See `docs/WhatsModified.md` for details.

| Item                                                          | Status   |
|---------------------------------------------------------------|----------|
| 1.1 MPI work partitioning                                     | **DONE** |
| 1.2 MPI immersed-wall reflection (Euler slip)                   | **DONE** |
| 1.3 Locator polygon closure (PointOnLine, Radius, image pts)  | **DONE** |
| 1.4 `HitOrMiss` PNPOLY rewrite                                | **DONE** |
| 1.5 Residual normalisation (RMS) — serial + MPI               | **DONE** |
| 1.6 Serial top/bottom farfield reference cell                 | **DONE** |
| 1.7 Mach number = `|v|/a` (factored into `updateAuxiliary`)   | **DONE** |
| 1.8 `MoversN` epsilon compare                                 | **DONE** |
| 2.1 Locator O(N⁴) → O(1) image neighbour lookup               | **DONE** |
| 2.2 MPI `isnan` checks out of `updateConservativeVariables`   | **DONE** |
| 2.3 OpenMP in hot kernels (gated by `_OPENMP`, Makefile)      | **DONE** |
| 2.4 `endl` → `'\n'` in solver output paths                    | **DONE** |
| MPI grid-read index alignment                                 | **DONE** |
| MPI `Nx`/`Ny` broadcast after IBM read                        | **DONE** |
| MPI per-rank output gather (`MPI_Gatherv`)                    | **DONE** |
| 3.1 Laminar Navier–Stokes viscous flux                        | **DONE** |
| 3.2 IBM no-slip wall (tied to `USE_NAVIER_STOKES`)            | **DONE** |
| 3.3 IBM adiabatic / fixed `T_w` thermal BC                    | **DONE** |
| 3.4 Viscous flux on fluid–IBM faces                           | **DONE** |
| 3.5 Optional WENO-5 inviscid interface (`USE_WENO5`)          | **DONE** |

Items still pending (longer-horizon): dynamic field allocation, 2D MPI decomposition, RK2 time integration, derived-datatype halo exchange, binary I/O, structure-of-arrays layout, characteristic WENO, turbulence models.

---

## 3. Navier–Stokes and IBM wall (implemented)

### 3.1 [DONE] Viscous flux

`ViscousFlux.h` + integration in both solvers. Default `USE_NAVIER_STOKES 1` in `SolverOptions.h`.

### 3.2 [DONE] IBM no-slip

When Navier–Stokes is enabled, `Gst==1` cells use bilinear image values then \(\mathbf{u}_\mathrm{ghost} = -\mathbf{u}_\mathrm{image}\). No separate no-slip flag — it is automatic with NS. Euler mode (`USE_NAVIER_STOKES 0`) retains slip via `ibmSlipVelocity`.

### 3.3 [DONE] Wall thermal BC

`IBM_WALL_ADIABATIC` (default) or `IBM_WALL_FIXED_TW` + `WALL_TEMPERATURE_TW` in `SolverOptions.h`.

### 3.4 [DONE] Fluid–wall viscous faces

`viscousFaceActive(Fl, Gst, …)` in `IbmBoundary.h` replaces fluid-only `viscousFaceOk`.

### 3.5 [DONE] Optional WENO-5

`USE_WENO5` in `SolverOptions.h`; `WENO5.h` + `HighOrderFlux.h`.

---

## 1. Critical correctness bugs

### 1.1 [BUG] MPI solver does not actually partition work between ranks

`MPI_Parallel_Solver/MOVERSNEulerSolverMPI.cpp`

- The residual loop checks `i >= start_i && i <= end_i` (line 443), but the actual computational kernels do **not**:
  - `calculateInterfaceFluxes()` (line 786): `for (i=2; i<=Nx; i++)`
  - `calculateDissipativeFluxes()` (line 884): same
  - `updateConservativeVariables()` (line 1020): same
  - `updatePrimitiveVariables()` (line 1123): same
  - `calculateTimeStep()` (line 989): same
- Result: every rank redundantly computes the entire global domain. The only “MPI” parts are `MPI_Allreduce` (residuals, dt) and ghost-cell exchange. So `np=4` is roughly the same wall-clock as `np=1`.

**Modification:** every kernel must loop only over its rank-owned slab `i ∈ [start_i, end_i]`, then exchange ghost layers between ranks. Pull out a single function `for_each_owned_cell(...)` that takes a body lambda, and route all kernels through it. Or pass `start_i, end_i` into each kernel.

### 1.2 [BUG] MPI immersed-boundary reflection is a different formula than serial

`MPI_Parallel_Solver/MOVERSNEulerSolverMPI.cpp`, `applyBoundaryConditions()` lines 763–772.

```cpp
Vtan = (u[im][jm] * ImNx[i][j] + v[im][jm] * ImNy[i][j]);   // actually u·n
Vnor = u[im][jm] - Vtan * ImNx[i][j];                       // tangential x-comp (mislabeled)
Vtangh = v[im][jm] - Vtan * ImNy[i][j];                     // tangential y-comp (mislabeled)
Vnorgh = -Vnor;
Vr  = Vnorgh + Vtan * ImNx[i][j];                           // = -u_x + 2 (u·n) n_x
Vrn = Vtangh + Vtan * ImNy[i][j];                           // = -u_y + 2 (u·n) n_y
```

That is `u' = -u + 2(u·n)n`, which **flips the tangential** component while keeping normal — the wrong direction for a slip wall. The correct inviscid (slip) reflection used by the serial code is `u' = u − 2(u·n)n`.

**Modification:** replace with the serial form (verified correct):

```cpp
double un  = u[im][jm] * ImNx[i][j] + v[im][jm] * ImNy[i][j]; // u · n
u[i][j]    = u[im][jm] - 2.0 * un * ImNx[i][j];
v[i][j]    = v[im][jm] - 2.0 * un * ImNy[i][j];
```

### 1.3 [BUG] Locator polygon is not closed

`MPI_Parallel_Solver/LocatorV4.1.cpp`

Three loops walk segments of the boundary polygon but skip the closing edge:

- Cell-on-line check (line 113): `for (k = 1; k < Nc; k++) PointOnLine(Cx[k], Cy[k], Cx[k+1], Cy[k+1], …)` → segments `(1,2) … (Nc-1,Nc)`. Closing edge `(Nc,1)` is missed.
- Radius (line 157): same pattern.
- Image points (line 220): same pattern.

The HitOrMiss ray-casting loop (line 132) is also subtly off:

```cpp
l = Nc - 1;
for (k = 1; k < Nc; k++) { Hits += HitOrMiss(Cx[l], Cy[l], Cx[k], Cy[k], ...); l = k; }
```

Iterates `(Nc-1,1), (1,2), …, (Nc-2,Nc-1)`. The edges `(Nc-1,Nc)` and `(Nc,1)` are missed. For a closed polygon this gives wrong inside/outside near those edges and miscomputes Mark/Rad there.

**Modification:** make all polygon walks close by indexing with modulo, e.g.

```cpp
for (int k = 1; k <= Nc; k++) {
    int kp1 = (k % Nc) + 1; // wrap Nc -> 1
    // use Cx[k], Cx[kp1]
}
```

### 1.4 [BUG] Residual normalization divides L2 sum by N (not √N)

Both serial and MPI compute `sqrt(Resn) / counter`. The standard L2 norm is `sqrt(Resn / counter)`. The current value is roughly `‖r‖₂ / N`, so it scales with grid size and is not directly comparable across resolutions.

`MOVERSNEulerSolverIBMV4.2.cpp` line 715, `MOVERSNEulerSolverMPI.cpp` lines 497–499.

**Modification:** `r_k = std::sqrt(Resn_k / counter)` (RMS) or `std::sqrt(Resn_k) / std::sqrt(counter)`.

### 1.5 [BUG] Inflow on top/bottom uses wrong reference cell

Serial `MOVERSNEulerSolverIBMV4.2.cpp` lines 388–414:

```cpp
j = 1;        for (i=1; i<=Nx+1; i++) Rho[i][j] = Rho[i+1][j]; // copies along i, on the boundary row
j = Ny+2;     for (i=1; i<=Nx+1; i++) Rho[i][j] = Rho[i-1][j]; // same
```

These copy from a cell on the same boundary row (`j=1` or `j=Ny+2`), which is itself a ghost; they should reflect from the first interior row (`j=2` or `j=Ny+1`). The current code creates a chain dependency where each cell copies the value of the next cell that has not been written yet, which is order-dependent and produces a smeared profile along the boundary.

**Modification:**

```cpp
// top farfield: copy from interior row j=2
j = 1;       for (i=1; i<=Nx+1; i++) { Rho[i][1]      = Rho[i][2];      ... }
// bottom farfield: copy from interior row j=Ny+1
j = Ny + 2;  for (i=1; i<=Nx+1; i++) { Rho[i][Ny+2] = Rho[i][Ny+1]; ... }
```

(Alternatively use freestream values as in the MPI version; pick one and use it consistently.)

### 1.6 [BUG] `Grid >> Nx >> Ny` then `while(!Grid.eof())` re-reads on EOF

Several places, e.g. `MOVERSNEulerSolverIBMV4.2.cpp` line 184. The `while(!eof())` idiom can iterate one extra time after EOF and read garbage. The fact that the inner double loop unconditionally re-reads makes this less harmful in practice, but it’s fragile.

**Modification:**

```cpp
ifstream Grid(config.gridFile);
if (!Grid) { ... }
Grid >> Nx >> Ny;
for (int j=2; j<=Ny+1; ++j)
  for (int i=2; i<=Nx+1; ++i)
    if (!(Grid >> X[i][j] >> Y[i][j]))
      throw std::runtime_error("Grid read failed at " + std::to_string(i) + "," + std::to_string(j));
```

### 1.7 [BUG] `MoversN` uses bare equality on doubles

`MathFunctions.h`:

```cpp
T S = (Ul != Ur) ? std::abs((Gl - Gr) / (Ul - Ur)) : std::abs(amin);
```

When `Ul` and `Ur` differ by < epsilon but are not bitwise equal, the divide can produce huge spurious `S`, which then gets clamped to `amax`. Not catastrophic (clamped), but it does add jitter and depends on rounding. Replace by:

```cpp
T denom = Ul - Ur;
T S = (std::fabs(denom) > T(1e-12)) ? std::fabs((Gl - Gr) / denom) : std::fabs(amin);
```

### 1.8 [BUG] `Locator::HitOrMiss` else branches both return 0

`LocatorFunctions.h` lines 153–161 — both inner branches return `0`. Likely meant `return 1` for the “line endpoints both right of px” case (or vice versa). As written the `Lin == 1` check is dead code.

**Modification:** redesign to canonical even–odd ray casting (e.g. PNPOLY):

```cpp
template <typename T>
int HitOrMiss(T x1, T y1, T x2, T y2, T px, T py) {
    bool cond = ((y1 > py) != (y2 > py)) &&
                (px < (x2 - x1) * (py - y1) / (y2 - y1) + x1);
    return cond ? 1 : 0;
}
```

This is shorter, robust, and handles closing edges correctly when combined with a wrap-around loop (see 1.3).

---

## 2. Major performance issues

### 2.1 [HOT] Locator IMAGE NEIGHBOUR POINTS is O(N⁴)

`LocatorV4.1.cpp` lines 246–263.

For each ghost cell (≈ N² cells), it loops over every `(im, jm)` cell to find which one contains the image point. For a 790×1185 grid (~9·10⁵ cells, with up to ~10⁴ ghosts), this is 10¹⁰ ops — minutes to hours.

Since `Kx, Ky` are cell centres of a Cartesian grid (uniform `dx, dy`), the lookup is O(1) per ghost cell:

```cpp
int im = (int)std::floor((ImX[i][j]   - (Xmin + 0.5*dx)) / dx) + 1;
int jm = (int)std::floor((ImY[i][j]   - (Ymin + 0.5*dy)) / dy) + 1;
ImBx[i][j] = im;  ImBy[i][j] = jm;
int gim = (int)std::floor((Geomx[i][j] - (Xmin + 0.5*dx)) / dx) + 1;
int gjm = (int)std::floor((Geomy[i][j] - (Ymin + 0.5*dy)) / dy) + 1;
ImGx[i][j] = gim; ImGy[i][j] = gjm;
```

Clamp to `[1, Nx-1]` × `[1, Ny-1]` to be safe. Expected speedup: **10⁵× or more** on the locator pre-process.

### 2.2 [HOT] MPI `MPI_Bcast` of full `MAX_GRID_X*MAX_GRID_Y` arrays

`MOVERSNEulerSolverMPI.cpp` lines ~298–312. Each broadcast moves `1500*1500*8 ≈ 18 MB` per array, and 16 arrays are broadcast — roughly 290 MB across all ranks. For small actual grids (e.g. 120×42) this is 99% wasted.

**Modification:**

- Read on rank 0, then broadcast `Nx, Ny` and only `Nx*Ny`-sized contiguous slices (or scatter them).
- Better: switch to dynamic 1D arrays sized `(Nx+pad)*(Ny+pad)` (see 3.2) and broadcast just that.

### 2.3 [HOT] NaN guards inside the hot kernel

`MOVERSNEulerSolverMPI.cpp::updateConservativeVariables()` (lines 1024–1110) does **32 `isnan` checks per fluid cell**, plus pressure/energy floors. This dominates the kernel cost and prevents vectorization. The serial version doesn’t do this.

**Modification:**

- Move the validation to a separate pass run every `validationFrequency` steps (e.g. 100 or 500).
- Inside the hot loop, write the unguarded update; rely on physical limiters at the very end (clamp `rho` and `P` only).
- If you hit NaN, log once and exit with a clear message rather than silently writing freestream.

### 2.4 [HOT] Fixed `[1500][1500]` static globals → ~1.7 GB BSS

Counting `MOVERSNEulerSolverIBMV4.2.cpp`: roughly 80+ `[1500][1500]` `double` arrays × 18 MB ≈ 1.4–1.7 GB. Loaded into the BSS at startup and uses a lot of cache footprint regardless of actual `Nx*Ny`.

**Modification:** dynamic allocation:

```cpp
struct Field2D {
    int nx, ny, stride;   // stride = ny + 2*pad
    std::vector<double> data;
    Field2D(int Nx, int Ny, int pad=2)
      : nx(Nx), ny(Ny), stride(Ny + 2*pad),
        data((Nx + 2*pad) * (Ny + 2*pad), 0.0) {}
    double& operator()(int i, int j) { return data[i*stride + j]; }
};
```

Wrap `Rho, P, u, v, …` in `Field2D` instances. Keeps the row-major `[i][j]` index pattern but sizes to actual `Nx, Ny`. Indexing as `arr(i,j)` improves cache density and lets you change grid size at runtime.

### 2.5 [HOT] Branchy `if (Fl[i][j] == 1)` in every hot kernel

Each fluid loop is dense in the `if` test, blocking auto-vectorization. Two improvements:

1. **Iterate over a precomputed cell list.** During locator/init, build `std::vector<std::pair<int,int>> fluid_cells` and iterate over it. Same correctness, much better branch prediction and vectorization.
2. **Replace `Fl` with packed integer mask** and use multiplication instead of branching for boundaries:

   ```cpp
   double m = (Fl[i][j] == 1) ? 1.0 : 0.0;   // or precomputed
   Un1[i][j] = U1[i][j] - m * ((dtg / Area[i][j]) * (...));
   ```

   This may help SIMD on Apple Silicon / x86.

### 2.6 [HOT] Verbose ASCII solution output

`Solution << ... << endl;` line by line for ~10⁶ cells. `endl` flushes every line, ~10× slower than `'\n'`. Also `iomanip` formatting is slow.

**Modification:**

- Replace `endl` with `'\n'` everywhere in solver/locator output paths.
- Use `fprintf` for formatted output, or write binary `.dat` then convert.
- Reduce `OUTPUT_FREQUENCY` defaults sensibly (every 1000–2000 steps is enough for movies).

### 2.7 [HOT] Single-stage Euler timestepping

The current scheme is forward Euler. Even RK2/RK3 (low-storage SSP) gives much better stability and can use larger CFL (factor 2 commonly), often >2× walltime improvement to convergence.

**Modification:** add a low-storage RK2 (Heun-style):

```
U^{*}     = U^n - dt L(U^n)
U^{n+1}   = 0.5 (U^n + U^{*}) - 0.5 dt L(U^{*})
```

It’s ~30 lines of code and reuses your existing flux routines.

### 2.8 [HOT] Recomputing identical sub-expressions in dissipative flux

`MOVERSNEulerSolverIBMV4.2.cpp` lines ~528–586 (and MPI mirror). Every `MoversN` call recomputes `Rho[i][j] * u[i][j] * nxr[i][j]` etc. Cache the corner state once at the top of the cell:

```cpp
const double r  = Rho[i][j], u_ = u[i][j], v_ = v[i][j], p = P[i][j], h = H[i][j], e = E[i][j];
const double rR = Rho[i+1][j], uR = u[i+1][j], …;
```

Then assemble fluxes from these locals. Already done partially in `calculateInterfaceFluxes`; extend to dissipative flux.

### 2.9 [HOT] OpenMP not used

The interior loops (flux, dissipation, update, primitive, residual) are embarrassingly parallel.

**Modification:**

- In `Makefile`, append `-fopenmp` to `CXXFLAGS` and `MPIFLAGS`.
- Annotate hot loops:

  ```cpp
  #pragma omp parallel for collapse(2) schedule(static)
  for (int j=2; j<=Ny+1; ++j)
    for (int i=2; i<=Nx;   ++i)
      if (Fl[i][j] == 1) { ... }
  ```

- For the residual loop, use reduction:

  ```cpp
  #pragma omp parallel for collapse(2) reduction(+:Resn1,Resn2,Resn3,Resn4,counter)
  ```

Expect 4–8× on a modern laptop, 16–64× on a server CPU.

---

## 3. Code structure / robustness

### 3.1 [STRUCT] Globals everywhere; near-zero unit testability

The solver state (~50 globals) is impossible to instantiate twice or mock for tests. Move state into a `SolverState` struct:

```cpp
struct SolverState {
    int Nx, Ny;
    Field2D Rho, P, u, v, ...;
    Field2D U1, U2, U3, U4, Un1, Un2, Un3, Un4;
    Field2D nxr, nyr, nxt, nyt, nxl, nyl, nxb, nyb, ...;
    Field2D Area, dt;
    SolverConfig cfg;
};
```

All kernels become free functions taking `SolverState&`. This makes testing, profiling, GPU porting, and 2D MPI decomposition far easier.

### 3.2 [STRUCT] Same fix in serial and MPI files

Many edits today need to be made in both `MPI_Parallel_Solver/` and `Geometries/`. Move everything into a single header-only library and have the two `.cpp` files just instantiate it with serial / MPI policies. Or build a common static lib.

### 3.3 [STRUCT] `std::pow(diff, 2)` still used in MPI residual

`MOVERSNEulerSolverMPI.cpp` lines 475–478 still use `pow(diff, 2)` — much slower than `diff*diff`. Apply the same fix as the serial version.

### 3.4 [STRUCT] Repeated `E, a, H, T, Mach` recomputation

The 7-line block:

```cpp
E[i][j] = (P/(rho*GAMMA_MINUS_1)) + 0.5*(u*u + v*v);
a[i][j] = sqrt(GAMMA*P/rho);
H[i][j] = (P/rho)*GAMMA_OVER_GAMMA_MINUS_1 + 0.5*(u*u + v*v);
T[i][j] = P/(rho*GAS_CONSTANT);
Mach[i][j] = u/a;
```

is duplicated 6+ times in serial and MPI. Factor into one inline:

```cpp
inline void updateAuxiliary(double rho, double u, double v, double P,
                            double& E, double& a, double& H, double& T, double& Mach) {
    const double q2 = u*u + v*v;
    E    = P/(rho*GAMMA_MINUS_1) + 0.5*q2;
    a    = std::sqrt(GAMMA*P/rho);
    H    = (P/rho)*GAMMA_OVER_GAMMA_MINUS_1 + 0.5*q2;
    T    = P/(rho*GAS_CONSTANT);
    Mach = std::sqrt(q2)/a;          // currently only u/a, which ignores v
}
```

Note also: serial code uses `Mach = u/a`, ignoring `v` → wrong on subsonic / off-axis flow. **Fix:** use velocity magnitude divided by `a`.

### 3.5 [STRUCT] MPI 1D decomposition and unused buffers

`send_buffer_top`, `send_buffer_bottom` are declared but never allocated or used; only left/right exchange exists. Either:

- delete top/bottom or implement 2D decomposition (recommended for large grids), and
- use `MPI_Type_create_subarray` to avoid manual pack/unpack.

### 3.6 [STRUCT] Output paths

All outputs are written to the current working directory. Add `--output-dir` (or `-o`) CLI argument and prepend it. Encourage runs from `MPI_Parallel_Solver/` while writing to `../output/<run_name>/`.

### 3.7 [STRUCT] CLI argument parsing

`config` is hard-coded. Accept simple args:

```text
./MOVERSNEulerSolverIBM [--cfl 0.3] [--max-iter 75000] [--out ../output/run1] \
                        [--grid GridPoints.dat] [--ibm ImmersedBoundary.dat]
```

10 lines with `getopt_long` or `argparse`.

### 3.8 [STRUCT] `validateFileOpen` always returns void semantics

The function takes a `const ifstream&` but the caller has already opened the file before calling (so the call is just a check). Replace with idiomatic C++:

```cpp
ifstream f(name);
if (!f) {
    std::cerr << "open failed: " << name << "\n";
    return -1;
}
```

### 3.9 [STRUCT] `isnan` checks scattered as `exit(1)`

Multiple places (`MOVERSNEulerSolverIBMV4.2.cpp` lines 658–677) call `exit(1)` from inside the loop. This bypasses MPI cleanup and ostreams. Better: throw, catch in `main`, finalize MPI, return non-zero.

### 3.10 [STRUCT] `Locator` corner cases

Behavior when `Nx*Ny > MAX_ARRAY_SIZE^2 = 9·10⁶` is undefined (silent overflow). Add bounds check before reading grid; `assert(Nc < MAX_NC)` and `assert(Nx < MAX_ARRAY_SIZE)`.

---

## 4. Smaller cleanups

- **[NIT]** `LocatorFunctions.h::Bilinear` parameter list is unusual (eight (x,y) pairs). Standard bilinear needs only `(x1,x2,y1,y2)` corners, four `q` values, and `(px,py)`. Simplify.
- **[NIT]** `LocatorFunctions.h::Radius`/`FootY` declare `qy` even when not needed.
- **[NIT]** Many `cout` lines remain commented-out for debug. Replace with a `verbose` flag.
- **[NIT]** Mixed `float` (X, Y, Cx, Cy, Area) and `double` (everything else). Choose one (probably `double` throughout) to avoid implicit conversions in flux assembly.
- **[NIT]** `using namespace std;` at file scope in headers and large translation units; prefer explicit `std::` qualifications.
- **[NIT]** Tabs vs spaces are mixed. Run clang-format with the project’s existing style.
- **[NIT]** `Boundary << "Zone T = Omega I = " << Nx << " J = " << Ny << endl;` — locator outputs Tecplot Zone metadata in `Boundary.dat`/`ImmersedBoundary.dat`, but the solver reads `Nx Ny` first and then plain rows. Document or unify the two formats.
- **[NIT]** `convert_to_vtk.cpp` parses `J = ` from a Zone string; brittle. Read the same Tecplot header the solver writes, with a small helper.

---

## 5. Suggested implementation order (low-risk → high-value)

1. **Fix MPI partitioning** (1.1) — single biggest correctness/perf win.
2. **Fix MPI immersed-wall reflection** (1.2).
3. **Fix Locator polygon closure and HitOrMiss** (1.3, 1.8) — solver reads Locator output, so silent geometry bugs cascade.
4. **Locator O(N⁴) → O(N²)** (2.1) — pre-processing time goes from minutes to seconds.
5. **Replace `pow(diff,2)` in MPI** (3.3) and **fix residual normalization** (1.4).
6. **Factor `updateAuxiliary` and fix Mach magnitude** (3.4).
7. **Dynamic allocation + `Field2D`** (2.4 / 3.1) — enables larger grids and SoA reorganization.
8. **OpenMP** (2.9) — single command line flag + ~20 pragmas.
9. **Cell-list iteration** (2.5).
10. **2D MPI decomposition with derived types** (3.5).
11. **RK2 timestepping** (2.7).
12. **CLI parsing, output dirs, structured config** (3.6, 3.7).

Items 1–6 should not change the discretization — they fix bugs and are validatable by comparing to the current serial output (which already runs to convergence on the bundled cases). Items 7–11 change the implementation but not the math.

---

## 6. Validation harness (recommended before/after large changes)

Add a `tests/` folder with:

- A small canonical case (e.g. `BluntBody` 120×42) wired to run for 200 iterations in <1 minute.
- A reference `Residues.dat` snapshot.
- A Python script that compares new residuals to reference within tolerance and prints diffs.

Then `make test` (or `ctest`) becomes the regression gate for every modification above.
