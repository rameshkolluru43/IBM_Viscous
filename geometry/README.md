# Geometry Definitions

This folder holds geometry definition files (`.txt`) for the Immersed Boundary solver.

- **BluntBody.txt**, **BluntBody_withAeroSpike.txt** – Blunt body and aerospike cases  
- **DBump.txt** – Double bump  
- **TriangularPrism.txt**, **Domain_TriangularPrism.txt** – Triangular prism and domain  

When running the locator or solver from `MPI_Parallel_Solver/` or `Geometries/`, point to these files, e.g.:

```sh
cd MPI_Parallel_Solver
./LocatorV4.1   # if it takes a geometry path, use ../geometry/BluntBody.txt
mpirun -np 4 ./MOVERSNEulerSolverMPI   # uses generated GridPoints.dat, ImmersedBoundary.dat in cwd
```

Generated grid/boundary files (`GridPoints.dat`, `ImmersedBoundary.dat`, etc.) are created in the solver directory when you run the locator/solver. Archived outputs are in `output/`.
