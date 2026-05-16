# Makefile for MPI Parallel CFD Solver
#
# Build modes:
#   make             - default optimized build (with OpenMP if compiler supports it)
#   make NO_OPENMP=1 - disable OpenMP (single-threaded loops)
#   make debug       - debug build (-g, -O0, no OpenMP)
#   make profile     - profiling build (-pg)
#
# OpenMP is enabled by default on Linux (GCC/Clang). On macOS the system clang has no
# OpenMP runtime; install libomp (`brew install libomp`) and pass OPENMP_FLAG='-Xpreprocessor
# -fopenmp -lomp' or use `make NO_OPENMP=1`.

# Compiler
MPICXX = mpicxx
CXX = g++

# OpenMP flag — auto-detect macOS (use libomp via Homebrew) vs Linux
ifeq ($(NO_OPENMP),1)
  OPENMP_FLAG =
else
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Darwin)
    LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
    ifeq ($(LIBOMP_PREFIX),)
      # libomp not installed; fall back to no-OpenMP rather than fail.
      OPENMP_FLAG =
    else
      OPENMP_FLAG = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include -L$(LIBOMP_PREFIX)/lib -lomp
    endif
  else
    OPENMP_FLAG ?= -fopenmp
  endif
endif

# Compiler flags
CXXFLAGS = -O3 -Wall -std=c++11 -march=native -mtune=native $(OPENMP_FLAG)
MPIFLAGS = -O3 -Wall -std=c++11 -march=native -mtune=native $(OPENMP_FLAG)

# Targets
MPI_TARGET = MOVERSNEulerSolverMPI
SERIAL_TARGET = MOVERSNEulerSolverIBM
LOCATOR_TARGET = LocatorV4.1
VTK_TARGET = convert_to_vtk

# Source files
MPI_SRC = MOVERSNEulerSolverMPI.cpp
SERIAL_SRC = MOVERSNEulerSolverIBMV4.2.cpp
LOCATOR_SRC = LocatorV4.1.cpp
VTK_SRC = convert_to_vtk.cpp

# Header files
HEADERS = MathFunctions.h LocatorFunctions.h SolverFunctions.h

# Default target
all: $(MPI_TARGET) $(SERIAL_TARGET) $(LOCATOR_TARGET) $(VTK_TARGET)

# MPI parallel solver
$(MPI_TARGET): $(MPI_SRC) $(HEADERS)
	$(MPICXX) $(MPIFLAGS) -o $(MPI_TARGET) $(MPI_SRC)

# Serial IBM solver
$(SERIAL_TARGET): $(SERIAL_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(SERIAL_TARGET) $(SERIAL_SRC)

# Locator utility
$(LOCATOR_TARGET): $(LOCATOR_SRC) $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $(LOCATOR_TARGET) $(LOCATOR_SRC)

# VTK converter (no MPI)
$(VTK_TARGET): $(VTK_SRC)
	$(CXX) $(CXXFLAGS) -o $(VTK_TARGET) $(VTK_SRC)

# Clean build files
clean:
	rm -f $(MPI_TARGET) $(SERIAL_TARGET) $(LOCATOR_TARGET) $(VTK_TARGET)
	rm -f *.o
	rm -f *~

# Run MPI solver with default parameters
run-mpi: $(MPI_TARGET)
	mpirun -np 4 ./$(MPI_TARGET)

# Run MPI solver with custom number of processes
run-mpi-np: $(MPI_TARGET)
	@echo "Usage: make run-mpi-np NP=<number_of_processes>"
	@if [ -z "$(NP)" ]; then \
		echo "Error: Please specify NP=<number_of_processes>"; \
		exit 1; \
	fi
	mpirun -np $(NP) ./$(MPI_TARGET)

# Run serial solver
run-serial: $(SERIAL_TARGET)
	./$(SERIAL_TARGET)

# Run locator
run-locator: $(LOCATOR_TARGET)
	./$(LOCATOR_TARGET)

# Debug build (no OpenMP, no aggressive opt)
debug: CXXFLAGS = -O0 -g -Wall -std=c++11 -DDEBUG
debug: MPIFLAGS = -O0 -g -Wall -std=c++11 -DDEBUG
debug: $(MPI_TARGET) $(SERIAL_TARGET)

# Profile build
profile: CXXFLAGS += -pg
profile: MPIFLAGS += -pg
profile: $(MPI_TARGET)

# Help
help:
	@echo "Available targets:"
	@echo "  all         - Build all executables"
	@echo "  $(MPI_TARGET) - Build MPI parallel solver"
	@echo "  $(SERIAL_TARGET) - Build serial IBM solver"
	@echo "  $(LOCATOR_TARGET) - Build locator utility"
	@echo "  $(VTK_TARGET) - Build VTK converter"
	@echo "  clean       - Remove build files"
	@echo "  run-mpi     - Run MPI solver with 4 processes"
	@echo "  run-mpi-np  - Run MPI solver with custom number of processes (use NP=<num>)"
	@echo "  run-serial  - Run serial solver"
	@echo "  run-locator - Run locator utility"
	@echo "  debug       - Build with debug flags"
	@echo "  profile     - Build with profiling flags"
	@echo "  help        - Show this help"

.PHONY: all clean run-mpi run-mpi-np run-serial run-locator debug profile help
