#!/usr/bin/env bash
# Supersonic aerospike (M=2): shock + viscous IBM wall — shock–boundary-layer interaction on the body.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CASE="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

if [ ! -f "$CASE/GridPoints.dat" ]; then
  echo "Missing $CASE/GridPoints.dat — copy from output/mpi_solver/"
  exit 1
fi

cp "$CASE"/GridPoints.dat "$CASE"/ImmersedBoundary.dat "$CASE"/Geometry.dat .

make NO_OPENMP=1 MOVERSNEulerSolverMPI

echo "Physics: Navier-Stokes, resolved no-slip (USE_WALL_FUNCTION=0), M=2"
echo "Grid: $(head -1 GridPoints.dat | tr '\t' ' ')"

STEPS="${1:-400}"
echo "Running up to $STEPS steps (2 MPI ranks, timeout 120s)..."

rm -f Residues_MPI.dat
timeout 120 mpirun -np 2 ./MOVERSNEulerSolverMPI || true

if [ -f Residues_MPI.dat ]; then
  N=$(wc -l < Residues_MPI.dat)
  echo "Completed $N residual lines"
  tail -3 Residues_MPI.dat
  if awk 'NR>1 && ($2=="nan" || $3=="nan") {exit 1}' Residues_MPI.dat; then
    echo "OK: no NaN in residuals"
  else
    echo "FAIL: NaN detected"
    exit 1
  fi
else
  echo "No Residues_MPI.dat produced"
  exit 1
fi
