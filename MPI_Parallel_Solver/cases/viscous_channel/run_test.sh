#!/usr/bin/env bash
# Coarse channel with IBM bottom wall: Navier-Stokes + log-law wall function (M=0.5).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CASE="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

echo "=== Build ==="
make NO_OPENMP=1 LocatorV4.1

echo "=== Locator: 80x24, domain [0,2]x[0,1] ==="
cp "$CASE/Geometry.dat" Geometry.dat
./LocatorV4.1 80 24 0 2 0 1

echo "=== Serial solvers (500 steps, M=0.5) ==="
perl -i -pe 's/timeStep < MAX_ITERATIONS \+ 1/timeStep < 501/' MOVERSNEulerSolverIBMV4.2.cpp
g++ -O3 -Wall -std=c++11 -DTEST_CHANNEL_MACH=0.5 -DUSE_WALL_FUNCTION=1 -o /tmp/solver_wf MOVERSNEulerSolverIBMV4.2.cpp
g++ -O3 -Wall -std=c++11 -DTEST_CHANNEL_MACH=0.5 -DUSE_WALL_FUNCTION=0 -o /tmp/solver_ns MOVERSNEulerSolverIBMV4.2.cpp
perl -i -pe 's/timeStep < 501/timeStep < MAX_ITERATIONS + 1/' MOVERSNEulerSolverIBMV4.2.cpp

/tmp/solver_wf
cp Aerospike_IBMSolution_500x500.dat "$CASE/solution_wallfunction.dat"
/tmp/solver_ns
cp Aerospike_IBMSolution_500x500.dat "$CASE/solution_noslip.dat"

echo "=== Residuals (wall function) ==="
tail -3 Residues.dat

echo "=== Near-wall u at x~1.0 (y, u_wf, u_ns) ==="
paste \
  <(awk 'NR>2 && $1>0.99 && $1<1.03 {print $2,$6}' "$CASE/solution_wallfunction.dat" | sort -n -u) \
  <(awk 'NR>2 && $1>0.99 && $1<1.03 {print $6}' "$CASE/solution_noslip.dat" | sort -n -u) \
  | awk 'NF>=2 && $1<0.25 {printf "  y=%-7.4f  u_wf=%-8.5f  u_ns=%-8.5f\n",$1,$2,$3}'

echo "Done."
