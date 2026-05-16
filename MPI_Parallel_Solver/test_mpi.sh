#!/bin/bash

# MPI CFD Solver Test Script
# This script compiles and runs basic tests for the MPI parallel CFD solver

echo "=========================================="
echo "MPI CFD Solver Test Script"
echo "=========================================="

# Check if MPI is available
echo "Checking MPI installation..."
if ! command -v mpirun &> /dev/null; then
    echo "ERROR: mpirun not found. Please install MPI (OpenMPI, MPICH, etc.)"
    exit 1
fi

if ! command -v mpicxx &> /dev/null; then
    echo "ERROR: mpicxx not found. Please install MPI development tools"
    exit 1
fi

echo "✓ MPI installation found"
mpirun --version | head -1

# Check for required input files
echo ""
echo "Checking for required input files..."
if [ ! -f "GridPoints.dat" ]; then
    echo "WARNING: GridPoints.dat not found"
    echo "You may need to run the locator program first to generate grid files"
fi

if [ ! -f "ImmersedBoundary.dat" ]; then
    echo "WARNING: ImmersedBoundary.dat not found"
    echo "You may need to run the locator program first to generate boundary files"
fi

# Clean previous builds
echo ""
echo "Cleaning previous builds..."
make clean

# Compile the MPI solver
echo ""
echo "Compiling MPI solver..."
if make MOVERSNEulerSolverMPI; then
    echo "✓ MPI solver compiled successfully"
else
    echo "✗ MPI solver compilation failed"
    exit 1
fi

# Compile the locator if input files are missing
if [ ! -f "GridPoints.dat" ] || [ ! -f "ImmersedBoundary.dat" ]; then
    echo ""
    echo "Compiling locator utility..."
    if make LocatorV4.1; then
        echo "✓ Locator utility compiled successfully"
        echo ""
        echo "Running locator to generate grid files..."
        if ./LocatorV4.1; then
            echo "✓ Grid files generated successfully"
        else
            echo "✗ Grid generation failed"
            exit 1
        fi
    else
        echo "✗ Locator compilation failed"
        exit 1
    fi
fi

# Test with different number of processes
echo ""
echo "Testing MPI solver with different process counts..."

# Test 1: Single process
echo ""
echo "Test 1: Running with 1 process..."
if timeout 30 mpirun -np 1 ./MOVERSNEulerSolverMPI > test_1proc.log 2>&1; then
    echo "✓ Single process test completed"
else
    echo "✗ Single process test failed or timed out"
    echo "Check test_1proc.log for details"
fi

# Test 2: Two processes
echo ""
echo "Test 2: Running with 2 processes..."
if timeout 30 mpirun -np 2 ./MOVERSNEulerSolverMPI > test_2proc.log 2>&1; then
    echo "✓ Two process test completed"
else
    echo "✗ Two process test failed or timed out"
    echo "Check test_2proc.log for details"
fi

# Test 3: Four processes (if enough cores available)
NPROC=$(nproc 2>/dev/null || echo 4)
if [ $NPROC -ge 4 ]; then
    echo ""
    echo "Test 3: Running with 4 processes..."
    if timeout 30 mpirun -np 4 ./MOVERSNEulerSolverMPI > test_4proc.log 2>&1; then
        echo "✓ Four process test completed"
    else
        echo "✗ Four process test failed or timed out"
        echo "Check test_4proc.log for details"
    fi
else
    echo ""
    echo "Skipping 4-process test (insufficient cores: $NPROC)"
fi

# Check output files
echo ""
echo "Checking output files..."
if [ -f "Residues_MPI.dat" ]; then
    echo "✓ Residuals file generated"
    echo "Last few lines of residuals:"
    tail -3 Residues_MPI.dat
else
    echo "✗ Residuals file not found"
fi

# Performance comparison (if serial version available)
echo ""
echo "Compiling serial version for comparison..."
if make MOVERSNEulerSolverIBM; then
    echo "✓ Serial solver compiled"
    
    echo ""
    echo "Performance comparison (10 iterations each):"
    echo "Serial version:"
    time timeout 10 ./MOVERSNEulerSolverIBM > serial_test.log 2>&1 || echo "Serial test completed/timed out"
    
    echo "MPI version (2 processes):"
    time timeout 10 mpirun -np 2 ./MOVERSNEulerSolverMPI > mpi_test.log 2>&1 || echo "MPI test completed/timed out"
else
    echo "Serial solver compilation failed - skipping performance comparison"
fi

# Summary
echo ""
echo "=========================================="
echo "Test Summary"
echo "=========================================="
echo "Test logs generated:"
ls -la *.log 2>/dev/null || echo "No log files found"

echo ""
echo "Output files:"
ls -la *.dat 2>/dev/null || echo "No output files found"

echo ""
echo "Executables built:"
ls -la MOVERSNEulerSolverMPI LocatorV4.1 MOVERSNEulerSolverIBM 2>/dev/null || echo "Check executable status"

echo ""
echo "=========================================="
echo "Testing completed!"
echo "=========================================="
echo ""
echo "Next steps:"
echo "1. Review test logs for any errors"
echo "2. Analyze residuals file for convergence"
echo "3. Run full simulation: mpirun -np <N> ./MOVERSNEulerSolverMPI"
echo "4. Visualize results using Tecplot, ParaView, or similar"
