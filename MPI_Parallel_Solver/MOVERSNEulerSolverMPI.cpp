/*****************************************************************************************/
/** 2D EULER CODE- MOVERS-N - MPI PARALLEL VERSION */
/** Shrinath K S */
/** Aero- IISC */
/** Parallelized with MPI domain decomposition */
/*****************************************************************************************/
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cstdlib>
#include <string>
#include <vector>
#include <limits>
#include "Config.h"
#include <mpi.h>

// Performance pragmas for compiler optimization (GCC only; Clang ignores)
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3,fast-math,unroll-loops")
#pragma GCC target("native")
#endif

#include "MathFunctions.h"
#include "Physics.h"
#include "HighOrderFlux.h"
#if USE_WENO5
#include "WENO5.h"
#endif
#include "IbmBoundary.h"
#if USE_NAVIER_STOKES
#include "ViscousFlux.h"
#endif
#include "LocatorFunctions.h"
using namespace std;

/*****************************************************************************************/
/** MPI Variables and Domain Decomposition                                                 */
/*                                                                                         */
/* Indexing convention (matches serial code):                                              */
/*   - Real fluid cells live at i = 2 .. Nx+1, j = 2 .. Ny+1                               */
/*   - Boundary/ghost cells at i = 0,1,Nx+2 and j = 0,1,Ny+2                               */
/*                                                                                         */
/* MPI 1D decomposition splits the i-direction across ranks.                               */
/*   - local_first_i, local_last_i: first/last *owned* interior cell on this rank          */
/*     (each cell is owned by exactly one rank — no double-counting)                       */
/*   - halo cells on each side hold copies of the neighbor rank's nearest interior cells   */
/*   - start_i/end_i (legacy): include `ghost_layers` halo columns on each side; kept      */
/*     for backwards compatibility with some helper functions                              */
/*****************************************************************************************/
int mpi_rank, mpi_size;
int local_Nx, local_Ny;
int start_i, end_i, start_j, end_j;
int local_first_i, local_last_i; /* owned interior i range; no halo, no overlap */
int ghost_layers = 2;

// MPI Communication buffers (top/bottom unused under 1D decomposition)
double *send_buffer_left = nullptr, *send_buffer_right = nullptr;
double *send_buffer_top  = nullptr, *send_buffer_bottom = nullptr;
double *recv_buffer_left = nullptr, *recv_buffer_right = nullptr;
double *recv_buffer_top  = nullptr, *recv_buffer_bottom = nullptr;

/*****************************************************************************************/
/** Forward Declarations */
/*****************************************************************************************/
void initializeMPI(int argc, char *argv[]);
void setupDomainDecomposition(int global_Nx, int global_Ny);
void exchangeGhostCells();
void calculateGridMetrics();
void applyBoundaryConditions();
void calculateInterfaceFluxes();
void calculateDissipativeFluxes();
#if USE_NAVIER_STOKES
void calculateViscousFluxes();
#endif
void updateConservativeVariables();
void updatePrimitiveVariables();
void calculateTimeStep();
void outputSolution(int timeStep);
void outputVTK(int timeStep);
void exportCurrentStateToVTK(const string &customName = "");
void gatherGlobalResiduals(double &global_res1, double &global_res2, double &global_res3, double &global_res4);
bool validateFileOpen(const std::ifstream &file, const std::string &filename);
bool isValidSolution();

/*****************************************************************************************/
/** Configuration Structure */
/*****************************************************************************************/
struct SolverConfiguration
{
	// Grid parameters
	int gridSizeX = 120;
	int gridSizeY = 42;

	// Flow conditions
	double machInfinity = 2.0;
	double densityInfinity = 1.4;
	double pressureInfinity = 1.0;

	// Solver parameters
	double cflNumber = CFL_NUMBER;
	int maxIterations = MAX_ITERATIONS;
	int outputFrequency = OUTPUT_FREQUENCY;
};

/*****************************************************************************************/
/** GLOBAL VARIABLES - Same as original but with MPI awareness */
/*****************************************************************************************/
static float X[MAX_GRID_X][MAX_GRID_Y], Y[MAX_GRID_X][MAX_GRID_Y], Cx[MAX_GRID_X][MAX_GRID_Y], Cy[MAX_GRID_X][MAX_GRID_Y], Area[MAX_GRID_X][MAX_GRID_Y];
static double nxr[MAX_GRID_X][MAX_GRID_Y], nxt[MAX_GRID_X][MAX_GRID_Y], nxl[MAX_GRID_X][MAX_GRID_Y], nxb[MAX_GRID_X][MAX_GRID_Y], dt[MAX_GRID_X][MAX_GRID_Y], dtg = 100, TG = 0;
static double nyr[MAX_GRID_X][MAX_GRID_Y], nyt[MAX_GRID_X][MAX_GRID_Y], nyl[MAX_GRID_X][MAX_GRID_Y], nyb[MAX_GRID_X][MAX_GRID_Y];
static double dsr[MAX_GRID_X][MAX_GRID_Y], dst[MAX_GRID_X][MAX_GRID_Y], dsl[MAX_GRID_X][MAX_GRID_Y], dsb[MAX_GRID_X][MAX_GRID_Y];
static double dxr[MAX_GRID_X][MAX_GRID_Y], dxt[MAX_GRID_X][MAX_GRID_Y], dxl[MAX_GRID_X][MAX_GRID_Y], dxb[MAX_GRID_X][MAX_GRID_Y];
static double dyr[MAX_GRID_X][MAX_GRID_Y], dyt[MAX_GRID_X][MAX_GRID_Y], dyl[MAX_GRID_X][MAX_GRID_Y], dyb[MAX_GRID_X][MAX_GRID_Y];
static double Rho[MAX_GRID_X][MAX_GRID_Y], P[MAX_GRID_X][MAX_GRID_Y], u[MAX_GRID_X][MAX_GRID_Y], v[MAX_GRID_X][MAX_GRID_Y], T[MAX_GRID_X][MAX_GRID_Y], E[MAX_GRID_X][MAX_GRID_Y], Mach[MAX_GRID_X][MAX_GRID_Y], H[MAX_GRID_X][MAX_GRID_Y];
static double U1[MAX_GRID_X][MAX_GRID_Y], U2[MAX_GRID_X][MAX_GRID_Y], U3[MAX_GRID_X][MAX_GRID_Y], U4[MAX_GRID_X][MAX_GRID_Y], a[MAX_GRID_X][MAX_GRID_Y];
static double Un1[MAX_GRID_X][MAX_GRID_Y], Un2[MAX_GRID_X][MAX_GRID_Y], Un3[MAX_GRID_X][MAX_GRID_Y], Un4[MAX_GRID_X][MAX_GRID_Y];
static double F1[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), F2[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), F3[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), F4[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused));
static double IF1R[MAX_GRID_X][MAX_GRID_Y], IF2R[MAX_GRID_X][MAX_GRID_Y], IF3R[MAX_GRID_X][MAX_GRID_Y], IF4R[MAX_GRID_X][MAX_GRID_Y];
static double IF1T[MAX_GRID_X][MAX_GRID_Y], IF2T[MAX_GRID_X][MAX_GRID_Y], IF3T[MAX_GRID_X][MAX_GRID_Y], IF4T[MAX_GRID_X][MAX_GRID_Y];
static double IF1L[MAX_GRID_X][MAX_GRID_Y], IF2L[MAX_GRID_X][MAX_GRID_Y], IF3L[MAX_GRID_X][MAX_GRID_Y], IF4L[MAX_GRID_X][MAX_GRID_Y];
static double IF1B[MAX_GRID_X][MAX_GRID_Y], IF2B[MAX_GRID_X][MAX_GRID_Y], IF3B[MAX_GRID_X][MAX_GRID_Y], IF4B[MAX_GRID_X][MAX_GRID_Y];
static double ImX[MAX_GRID_X][MAX_GRID_Y], ImY[MAX_GRID_X][MAX_GRID_Y], ImNx[MAX_GRID_X][MAX_GRID_Y], ImNy[MAX_GRID_X][MAX_GRID_Y];
static double DF1R[MAX_GRID_X][MAX_GRID_Y], DF2R[MAX_GRID_X][MAX_GRID_Y], DF3R[MAX_GRID_X][MAX_GRID_Y], DF4R[MAX_GRID_X][MAX_GRID_Y];
static double DF1T[MAX_GRID_X][MAX_GRID_Y], DF2T[MAX_GRID_X][MAX_GRID_Y], DF3T[MAX_GRID_X][MAX_GRID_Y], DF4T[MAX_GRID_X][MAX_GRID_Y];
static double DF1L[MAX_GRID_X][MAX_GRID_Y], DF2L[MAX_GRID_X][MAX_GRID_Y], DF3L[MAX_GRID_X][MAX_GRID_Y], DF4L[MAX_GRID_X][MAX_GRID_Y];
static double DF1B[MAX_GRID_X][MAX_GRID_Y], DF2B[MAX_GRID_X][MAX_GRID_Y], DF3B[MAX_GRID_X][MAX_GRID_Y], DF4B[MAX_GRID_X][MAX_GRID_Y];
#if USE_NAVIER_STOKES
static double VF1R[MAX_GRID_X][MAX_GRID_Y], VF2R[MAX_GRID_X][MAX_GRID_Y], VF3R[MAX_GRID_X][MAX_GRID_Y], VF4R[MAX_GRID_X][MAX_GRID_Y];
static double VF1T[MAX_GRID_X][MAX_GRID_Y], VF2T[MAX_GRID_X][MAX_GRID_Y], VF3T[MAX_GRID_X][MAX_GRID_Y], VF4T[MAX_GRID_X][MAX_GRID_Y];
static double VF1L[MAX_GRID_X][MAX_GRID_Y], VF2L[MAX_GRID_X][MAX_GRID_Y], VF3L[MAX_GRID_X][MAX_GRID_Y], VF4L[MAX_GRID_X][MAX_GRID_Y];
static double VF1B[MAX_GRID_X][MAX_GRID_Y], VF2B[MAX_GRID_X][MAX_GRID_Y], VF3B[MAX_GRID_X][MAX_GRID_Y], VF4B[MAX_GRID_X][MAX_GRID_Y];
#endif
static double G1L __attribute__((unused)), G2L __attribute__((unused)), G3L __attribute__((unused)), G4L __attribute__((unused)), G1R __attribute__((unused)), G2R __attribute__((unused)), G3R __attribute__((unused)), G4R __attribute__((unused));
static double a_maxR[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), a_maxT[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), a_maxL[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), a_maxB[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused));
static double a_minR[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), a_minT[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), a_minL[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused)), a_minB[MAX_GRID_X][MAX_GRID_Y] __attribute__((unused));
static double Geomx[MAX_GRID_X][MAX_GRID_Y], Geomy[MAX_GRID_X][MAX_GRID_Y], pg __attribute__((unused)), ug __attribute__((unused)), rhog __attribute__((unused)), ag __attribute__((unused)), vg __attribute__((unused)), Machg __attribute__((unused));
static double Fl[MAX_GRID_X][MAX_GRID_Y], Gst[MAX_GRID_X][MAX_GRID_Y], Rad[MAX_GRID_X][MAX_GRID_Y];
static double dx __attribute__((unused)) = 0, dy __attribute__((unused)) = 0, a1 __attribute__((unused)), b1 __attribute__((unused)), dy1 __attribute__((unused)), dy2 __attribute__((unused)), delx __attribute__((unused)), dely __attribute__((unused)), delxy __attribute__((unused));
static double Vtan __attribute__((unused)), Vnor __attribute__((unused)),
              Vtangh __attribute__((unused)), Vnorgh __attribute__((unused)),
              Vr __attribute__((unused)), Vrn __attribute__((unused)),
              C __attribute__((unused)), L = CFL_NUMBER;
static double Ct __attribute__((unused)), Cr __attribute__((unused)),
              Cl __attribute__((unused)), Cb __attribute__((unused)),
              alpha __attribute__((unused)), alpha1, alpha2, alpha3, alpha4,
              dex __attribute__((unused)), dey __attribute__((unused)),
              amax, amin;
static double c1 __attribute__((unused)) = 0, Resn1 = 0, Resn2 = 0, Resn3 = 0, Resn4 = 0;
static int i, j, timeStep, Nx = 120, Ny = 42, im, jm, ImBx[MAX_GRID_X][MAX_GRID_Y], ImBy[MAX_GRID_X][MAX_GRID_Y], ImGx[MAX_GRID_X][MAX_GRID_Y], ImGy[MAX_GRID_X][MAX_GRID_Y];

// Flow conditions
static double M_inf = (TEST_CHANNEL_MACH > 0.0) ? TEST_CHANNEL_MACH : 2.0;
static double rho_inf = 1.4;
static double p_inf = 1.0;
static double T_inf = p_inf / (rho_inf * GAS_CONSTANT);
static double a_inf = sqrt(GAMMA * GAS_CONSTANT * T_inf);
static double u_inf = M_inf * a_inf;
static double v_inf = 0.0;
static double H_inf __attribute__((unused)) = (p_inf / rho_inf) * GAMMA_OVER_GAMMA_MINUS_1 + 0.5 * (u_inf * u_inf + v_inf * v_inf);
static double E_inf = (p_inf / (rho_inf * GAMMA_MINUS_1)) + 0.5 * (u_inf * u_inf + v_inf * v_inf);

/*****************************************************************************************/
/** MPI INITIALIZATION AND DOMAIN DECOMPOSITION */
/*****************************************************************************************/
void initializeMPI(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (mpi_rank == 0)
	{
		cout << "=== MPI PARALLEL CFD SOLVER ===" << '\n';
		cout << "Number of MPI processes: " << mpi_size << '\n';
#if USE_NAVIER_STOKES
		cout << "Physics: Navier-Stokes (mu=" << LAMINAR_VISCOSITY << "), IBM wall";
#if USE_WALL_FUNCTION
		cout << " (log-law wall function)";
#else
		cout << " (resolved no-slip)";
#endif
#if IBM_WALL_FIXED_TW
		cout << ", fixed T_w=" << WALL_TEMPERATURE_TW;
#else
		cout << ", adiabatic";
#endif
		cout << '\n';
#else
		cout << "Physics: inviscid Euler, IBM slip wall" << '\n';
#endif
	}
}

void setupDomainDecomposition(int global_Nx, int global_Ny)
{
	/* Partition the interior i-range [2 .. global_Nx+1] across ranks. Each rank owns
	 * a contiguous slab; remainders are distributed to the lowest-numbered ranks.
	 * Halo layers (ghost_layers cells) are kept in addition to owned cells.       */
	const int interior_size = global_Nx; /* there are global_Nx interior cells in i */
	const int base = interior_size / mpi_size;
	const int rem  = interior_size % mpi_size;

	local_Nx = base + (mpi_rank < rem ? 1 : 0);
	local_Ny = global_Ny;

	const int offset = mpi_rank * base + std::min(mpi_rank, rem);
	local_first_i = 2 + offset;
	local_last_i  = local_first_i + local_Nx - 1;

	/* Legacy halo-inclusive bounds (some helpers still use these). */
	start_i = local_first_i;
	end_i   = local_last_i;
	if (mpi_rank > 0)
		start_i -= ghost_layers;
	if (mpi_rank < mpi_size - 1)
		end_i += ghost_layers;

	start_j = 1;
	end_j   = global_Ny;

	if (mpi_rank == 0)
	{
		cout << "Domain decomposition setup complete." << '\n';
		cout << "Per-rank owned i-cells: " << local_Nx << " (j: " << local_Ny << ")" << '\n';
	}

	/* Communication buffer: U1..U4 + Rho,u,v,P,H (+ T when Navier–Stokes). */
	const int ghost_vars = USE_NAVIER_STOKES ? 10 : 9;
	const int buffer_size = ghost_layers * (global_Ny + 3) * ghost_vars;
	if (send_buffer_left)  delete[] send_buffer_left;
	if (send_buffer_right) delete[] send_buffer_right;
	if (recv_buffer_left)  delete[] recv_buffer_left;
	if (recv_buffer_right) delete[] recv_buffer_right;
	send_buffer_left  = new double[buffer_size];
	send_buffer_right = new double[buffer_size];
	recv_buffer_left  = new double[buffer_size];
	recv_buffer_right = new double[buffer_size];
}

/* Exchange ghost-cell columns with neighboring ranks.
 *
 * We pack ghost_layers columns from the owned region into the send buffer,
 * exchange with neighbors, and unpack into the halo columns. Each j-column
 * carries 9 fields: U1..U4 (conservative) AND Rho, u, v, P, H (primitive),
 * because the flux kernels use the *primitive* variables of the neighbor. */
void exchangeGhostCells()
{
	const int J_LO = 0;
	const int J_HI = Ny + 2;        /* full j range including boundary layer */
	const int J_LEN = J_HI - J_LO + 1;
	const int VARS = USE_NAVIER_STOKES ? 10 : 9;
	const int col_size = J_LEN * VARS;

	auto pack = [&](double *buf, int i_col) {
		int idx = 0;
		for (int jj = J_LO; jj <= J_HI; ++jj) {
			buf[idx++] = U1[i_col][jj];
			buf[idx++] = U2[i_col][jj];
			buf[idx++] = U3[i_col][jj];
			buf[idx++] = U4[i_col][jj];
			buf[idx++] = Rho[i_col][jj];
			buf[idx++] = u[i_col][jj];
			buf[idx++] = v[i_col][jj];
			buf[idx++] = P[i_col][jj];
			buf[idx++] = H[i_col][jj];
			if (USE_NAVIER_STOKES)
				buf[idx++] = T[i_col][jj];
		}
	};
	auto unpack = [&](const double *buf, int i_col) {
		int idx = 0;
		for (int jj = J_LO; jj <= J_HI; ++jj) {
			U1[i_col][jj]  = buf[idx++];
			U2[i_col][jj]  = buf[idx++];
			U3[i_col][jj]  = buf[idx++];
			U4[i_col][jj]  = buf[idx++];
			Rho[i_col][jj] = buf[idx++];
			u[i_col][jj]   = buf[idx++];
			v[i_col][jj]   = buf[idx++];
			P[i_col][jj]   = buf[idx++];
			H[i_col][jj]   = buf[idx++];
			if (USE_NAVIER_STOKES)
				T[i_col][jj] = buf[idx++];
		}
	};

	MPI_Request requests[4];
	MPI_Status statuses[4];
	int req_count = 0;

	/* Send ghost_layers leftmost owned columns to left neighbor; receive into the
	 * halo columns just left of the owned region. */
	if (mpi_rank > 0) {
		for (int g = 0; g < ghost_layers; ++g)
			pack(send_buffer_left + g * col_size, local_first_i + g);
		const int n = ghost_layers * col_size;
		MPI_Isend(send_buffer_left,  n, MPI_DOUBLE, mpi_rank - 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Irecv(recv_buffer_left,  n, MPI_DOUBLE, mpi_rank - 1, 1, MPI_COMM_WORLD, &requests[req_count++]);
	}
	/* Mirror exchange with right neighbor. */
	if (mpi_rank < mpi_size - 1) {
		for (int g = 0; g < ghost_layers; ++g)
			pack(send_buffer_right + g * col_size, local_last_i - ghost_layers + 1 + g);
		const int n = ghost_layers * col_size;
		MPI_Isend(send_buffer_right, n, MPI_DOUBLE, mpi_rank + 1, 1, MPI_COMM_WORLD, &requests[req_count++]);
		MPI_Irecv(recv_buffer_right, n, MPI_DOUBLE, mpi_rank + 1, 0, MPI_COMM_WORLD, &requests[req_count++]);
	}

	MPI_Waitall(req_count, requests, statuses);

	if (mpi_rank > 0) {
		for (int g = 0; g < ghost_layers; ++g)
			unpack(recv_buffer_left + g * col_size, local_first_i - ghost_layers + g);
	}
	if (mpi_rank < mpi_size - 1) {
		for (int g = 0; g < ghost_layers; ++g)
			unpack(recv_buffer_right + g * col_size, local_last_i + 1 + g);
	}
}

void gatherGlobalResiduals(double &global_res1, double &global_res2, double &global_res3, double &global_res4)
{
	double local_residuals[4] = {Resn1, Resn2, Resn3, Resn4};
	double global_residuals[4];

	MPI_Allreduce(local_residuals, global_residuals, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	global_res1 = global_residuals[0];
	global_res2 = global_residuals[1];
	global_res3 = global_residuals[2];
	global_res4 = global_residuals[3];
}

/*****************************************************************************************/
/** MAIN FUNCTION - MPI PARALLEL VERSION */
/*****************************************************************************************/
int main(int argc, char *argv[])
{
	// Initialize MPI
	initializeMPI(argc, argv);

	// Setup domain decomposition
	setupDomainDecomposition(Nx, Ny);

	// Only rank 0 reads grid and broadcasts
	if (mpi_rank == 0)
	{
		cout << "Reading Grid from the Grid Data" << '\n';
		ifstream Grid("GridPoints.dat");
		if (!validateFileOpen(Grid, "GridPoints.dat"))
		{
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}

		Grid >> Nx >> Ny;
		cout << "Grid Size is" << Nx << "x" << Ny << '\n';

		/* The grid file uses POINT counts (Nx, Ny here are #points = #cells + 1).
		 * The solver expects grid points stored at indices i=2..Nx+1, j=2..Ny+1 so
		 * that cell (i, j) takes its lower-left corner at X[i][j].
		 * (This matches the serial solver's read pattern.) */
		for (j = 2; j <= Ny + 1; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				Grid >> X[i][j] >> Y[i][j];
			}
		}
		Grid.close();
	}

	// Broadcast grid dimensions and coordinates
	MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(X, (MAX_GRID_X * MAX_GRID_Y), MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Y, (MAX_GRID_X * MAX_GRID_Y), MPI_FLOAT, 0, MPI_COMM_WORLD);

	// Update domain decomposition with actual grid size
	setupDomainDecomposition(Nx, Ny);

	// Read immersed boundary data (rank 0) and broadcast
	if (mpi_rank == 0)
	{
		cout << "Reading Cell Centres and Immersed Boundary Data" << '\n';
		ifstream IBM("ImmersedBoundary.dat");
		if (!validateFileOpen(IBM, "ImmersedBoundary.dat"))
		{
			MPI_Abort(MPI_COMM_WORLD, 1);
			return 1;
		}

		IBM >> Nx >> Ny;
		cout << "Grid Size of Cells is" << Nx << "x" << Ny << '\n';

		for (j = 2; j <= Ny + 1; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				IBM >> Cx[i][j] >> Cy[i][j] >> Fl[i][j] >> Gst[i][j] >> Rad[i][j] >> ImNx[i][j] >> ImNy[i][j] >> ImX[i][j] >> ImY[i][j] >> ImBx[i][j] >> ImBy[i][j] >> Geomx[i][j] >> Geomy[i][j] >> ImGx[i][j] >> ImGy[i][j];
			}
		}
		IBM.close();
	}

	/* The IBM file may report a different (typically smaller by 1) Nx/Ny than the
	 * grid file — broadcast them so every rank has consistent dimensions. */
	MPI_Bcast(&Nx, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&Ny, 1, MPI_INT, 0, MPI_COMM_WORLD);
	setupDomainDecomposition(Nx, Ny);

	// Broadcast immersed boundary data
	MPI_Bcast(Cx, (MAX_GRID_X * MAX_GRID_Y), MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Cy, (MAX_GRID_X * MAX_GRID_Y), MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Area, (MAX_GRID_X * MAX_GRID_Y), MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Fl, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Gst, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Rad, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImX, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImY, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImNx, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImNy, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImBx, (MAX_GRID_X * MAX_GRID_Y), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImBy, (MAX_GRID_X * MAX_GRID_Y), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(Geomx, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(Geomy, (MAX_GRID_X * MAX_GRID_Y), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImGx, (MAX_GRID_X * MAX_GRID_Y), MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(ImGy, (MAX_GRID_X * MAX_GRID_Y), MPI_INT, 0, MPI_COMM_WORLD);

	// Calculate grid metrics
	calculateGridMetrics();

	// Initialize local domain (simplified version - full implementation would follow original pattern)
	if (mpi_rank == 0)
	{
		cout << "Domain Initialized with flow perturbation (MPI Parallel)" << '\n';
	}

	/* Initialize flow field over each rank's slab + halo. The perturbation depends
	 * only on (X, Y) which are identical across ranks, so each rank produces a
	 * consistent halo without explicit communication. We init the full j range
	 * (0..Ny+2) including boundary rows so primitives are valid in the first
	 * stencil access. */
	for (j = 0; j <= Ny + 2; j++)
	{
		for (i = std::max(0, start_i); i <= std::min(Nx + 2, end_i); i++)
		{
			if (std::isnan(X[i][j]) || std::isnan(Y[i][j]) || std::isinf(X[i][j]) || std::isinf(Y[i][j]))
			{
				if (mpi_rank == 0)
					cout << "ERROR: Invalid grid coordinates at (" << i << "," << j << "): X=" << X[i][j] << " Y=" << Y[i][j] << '\n';
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}

			Rho[i][j] = rho_inf;
			u[i][j]   = u_inf;
			v[i][j]   = v_inf;
			P[i][j]   = p_inf;
			updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
			                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
			U1[i][j] = rho_inf;
			U2[i][j] = rho_inf * u_inf;
			U3[i][j] = rho_inf * v_inf;
			U4[i][j] = rho_inf * E_inf;

			// Add small perturbation near body to trigger flow development
			const double x_pos = X[i][j];
			const double y_pos = Y[i][j];
			const double dist_from_body = std::sqrt(x_pos * x_pos + y_pos * y_pos);
			if (dist_from_body > 0.1 && dist_from_body < 2.0)
			{
				const double perturbation_factor = 1.0 + 0.01 * std::exp(-dist_from_body) * std::cos(PI * x_pos);
				P[i][j]   = p_inf * perturbation_factor;
				Rho[i][j] = rho_inf * perturbation_factor;
				updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
				                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
				U1[i][j] = Rho[i][j];
				U2[i][j] = Rho[i][j] * u[i][j];
				U3[i][j] = Rho[i][j] * v[i][j];
				U4[i][j] = Rho[i][j] * E[i][j];
			}
		}
	}

	// Initialize Un arrays to current U values for first iteration residual calculation
	for (j = 1; j <= Ny + 2; j++)
	{
		for (i = 1; i <= Nx + 1; i++)
		{
			Un1[i][j] = U1[i][j];
			Un2[i][j] = U2[i][j];
			Un3[i][j] = U3[i][j];
			Un4[i][j] = U4[i][j];
		}
	}

	// Open residual file (only rank 0)
	ofstream Residual;
	if (mpi_rank == 0)
	{
		Residual.open("Residues_MPI.dat");
	}

	// Main time-stepping loop
	for (timeStep = 1; timeStep < MAX_ITERATIONS + 1; timeStep++)
	{
		/* Order:
		 *   apply BCs -> exchange halo -> compute fluxes -> compute dt ->
		 *   update conservatives -> compute residual (Un - U_old) ->
		 *   update primitives (which also does U <- Un, advancing the time level).  */
		applyBoundaryConditions();
		exchangeGhostCells();

		calculateInterfaceFluxes();
		calculateDissipativeFluxes();
#if USE_NAVIER_STOKES
		calculateViscousFluxes();
#endif
		calculateTimeStep();
		updateConservativeVariables();

		/* Per-rank residual over owned cells (no halo, no overlap). */
		Resn1 = Resn2 = Resn3 = Resn4 = 0.0;
		int counter = 0;
		const int i_beg = std::max(2, local_first_i);
		const int i_end = std::min(Nx, local_last_i);
#if defined(_OPENMP)
		#pragma omp parallel for collapse(2) schedule(static) \
		    reduction(+:Resn1,Resn2,Resn3,Resn4,counter)
#endif
		for (int j = 2; j <= Ny + 1; j++)
		{
			for (int i = i_beg; i <= i_end; i++)
			{
				if (Fl[i][j] == 1)
				{
					counter++;
					const double d1 = Un1[i][j] - U1[i][j];
					const double d2 = Un2[i][j] - U2[i][j];
					const double d3 = Un3[i][j] - U3[i][j];
					const double d4 = Un4[i][j] - U4[i][j];
					Resn1 += d1 * d1;
					Resn2 += d2 * d2;
					Resn3 += d3 * d3;
					Resn4 += d4 * d4;
				}
			}
		}

		double global_res1, global_res2, global_res3, global_res4;
		gatherGlobalResiduals(global_res1, global_res2, global_res3, global_res4);
		int global_counter;
		MPI_Allreduce(&counter, &global_counter, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

		if (mpi_rank == 0)
		{
			const double invN = (global_counter > 0) ? 1.0 / global_counter : 1.0;
			const double r1 = std::sqrt(global_res1 * invN);
			const double r2 = std::sqrt(global_res2 * invN);
			const double r3 = std::sqrt(global_res3 * invN);
			const double r4 = std::sqrt(global_res4 * invN);
			Residual << timeStep << '\t' << r1 << '\t' << r2 << '\t' << r3 << '\t' << r4 << '\t' << TG << '\n';
			if (timeStep == 1 || timeStep % 1000 == 0)
			{
				cout << "MPI Rank 0: " << timeStep << '\t' << r1 << '\t' << r2 << '\t' << r3 << '\t' << r4 << '\t' << TG << '\n';
			}
		}

		updatePrimitiveVariables();

		// Validate solution periodically (cheap; occasional check for divergence)
		if (timeStep % VALIDATION_FREQUENCY == 0)
		{
			int local_ok = isValidSolution() ? 1 : 0;
			int global_ok;
			MPI_Allreduce(&local_ok, &global_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			if (!global_ok)
			{
				if (mpi_rank == 0)
					cout << "Solution validation failed at time step " << timeStep << '\n';
				MPI_Abort(MPI_COMM_WORLD, 1);
				return 1;
			}
		}

		if (timeStep % OUTPUT_FREQUENCY == 0)
		{
			outputSolution(timeStep);
			outputVTK(timeStep);
		}

		TG += dtg; // Time advancement
	}

	// Cleanup
	if (mpi_rank == 0)
	{
		Residual.close();
		cout << "MPI Parallel CFD simulation completed!" << '\n';
	}

	// Cleanup MPI buffers
	delete[] send_buffer_left;
	delete[] send_buffer_right;
	delete[] recv_buffer_left;
	delete[] recv_buffer_right;

	MPI_Finalize();
	return 0;
}

/*****************************************************************************************/
/** UTILITY FUNCTIONS */
/*****************************************************************************************/
bool validateFileOpen(const std::ifstream &file, const std::string &filename)
{
	if (!file.is_open())
	{
		if (mpi_rank == 0)
		{
			cerr << "Error: Could not open file " << filename << '\n';
		}
		return false;
	}
	return true;
}

/*****************************************************************************************/
/** GRID METRICS CALCULATION                                                              */
/* Computed once at startup; we keep the full Nx loop here because every rank needs the   */
/* metrics for its halo cells too, and the cost is amortized over all timesteps.          */
/*****************************************************************************************/
void calculateGridMetrics()
{
	// Calculate cell areas and side lengths
	for (i = 2; i <= Nx + 1; i++)
	{
		for (j = 2; j <= Ny + 1; j++)
		{
			Area[i][j] = 0.5 * fabs(((X[i][j + 1] - X[i + 1][j]) * (Y[i][j] - Y[i + 1][j + 1])) -
									((Y[i][j + 1] - Y[i + 1][j]) * (X[i][j] - X[i + 1][j + 1])));
		}
	}

	// Calculate side lengths and normal vectors
	for (i = 2; i <= Nx + 1; i++)
	{
		for (j = 2; j <= Ny + 1; j++)
		{
			dxr[i][j] = X[i + 1][j + 1] - X[i + 1][j];
			dxt[i][j] = X[i][j + 1] - X[i + 1][j + 1];
			dyr[i][j] = Y[i + 1][j + 1] - Y[i + 1][j];
			dyt[i][j] = Y[i][j + 1] - Y[i + 1][j + 1];
			dxl[i][j] = X[i][j] - X[i][j + 1];
			dxb[i][j] = X[i + 1][j] - X[i][j];
			dyl[i][j] = Y[i][j] - Y[i][j + 1];
			dyb[i][j] = Y[i + 1][j] - Y[i][j];

			dsr[i][j] = sqrt(dxr[i][j] * dxr[i][j] + dyr[i][j] * dyr[i][j]);
			dst[i][j] = sqrt(dxt[i][j] * dxt[i][j] + dyt[i][j] * dyt[i][j]);
			dsl[i][j] = sqrt(dxl[i][j] * dxl[i][j] + dyl[i][j] * dyl[i][j]);
			dsb[i][j] = sqrt(dxb[i][j] * dxb[i][j] + dyb[i][j] * dyb[i][j]);

			nxr[i][j] = dyr[i][j] / dsr[i][j];
			nyr[i][j] = -dxr[i][j] / dsr[i][j];
			nxt[i][j] = dyt[i][j] / dst[i][j];
			nyt[i][j] = -dxt[i][j] / dst[i][j];
			nxl[i][j] = dyl[i][j] / dsl[i][j];
			nyl[i][j] = -dxl[i][j] / dsl[i][j];
			nxb[i][j] = dyb[i][j] / dsb[i][j];
			nyb[i][j] = -dxb[i][j] / dsb[i][j];
		}
	}
}

/*****************************************************************************************/
/** BOUNDARY CONDITIONS
 *
 * BC application is rank-aware:
 *   - INLET (i=1) is owned by rank 0 only.
 *   - OUTLET (i=Nx+2) is owned by the last rank only.
 *   - Top/bottom farfield (j=1, j=Ny+2) is applied by every rank, but only over
 *     the i-cells it owns (or its halo cells that fall within the global domain).
 *   - Immersed-boundary ghost cells are applied per rank over its owned i-range.
 *
 * Navier–Stokes: IBM no-slip + adiabatic or fixed T_w (see SolverOptions.h).
 * Euler (USE_NAVIER_STOKES=0): inviscid slip wall u' = u - 2 (u . n) n.
 */
/*****************************************************************************************/
void applyBoundaryConditions()
{
	const int i_lo = std::max(1, local_first_i - ghost_layers);
	const int i_hi = std::min(Nx + 2, local_last_i + ghost_layers);

	// INLET - supersonic inflow (only rank that owns column 1's neighbor)
	if (mpi_rank == 0)
	{
		i = 1;
		for (j = 1; j <= Ny + 1; j++)
		{
			Rho[i][j] = rho_inf;
			u[i][j]   = u_inf;
			v[i][j]   = v_inf;
			P[i][j]   = p_inf;
			updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
			                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
		}
	}

	// OUTLET - supersonic outflow extrapolation (only on last rank)
	if (mpi_rank == mpi_size - 1)
	{
		i = Nx + 2;
		for (j = 1; j <= Ny + 1; j++)
		{
			Rho[i][j] = Rho[Nx + 1][j];
			u[i][j]   = u[Nx + 1][j];
			v[i][j]   = v[Nx + 1][j];
			P[i][j]   = P[Nx + 1][j];
			updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
			                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
		}
	}

	// FARFIELD - top boundary (j=1) and bottom boundary (j=Ny+2)
	for (int side = 0; side < 2; ++side)
	{
		const int jb = (side == 0) ? 1 : (Ny + 2);
		for (i = i_lo; i <= i_hi; i++)
		{
			Rho[i][jb] = rho_inf;
			u[i][jb]   = u_inf;
			v[i][jb]   = v_inf;
			P[i][jb]   = p_inf;
			updateAuxiliary(Rho[i][jb], u[i][jb], v[i][jb], P[i][jb],
			                E[i][jb], a[i][jb], H[i][jb], T[i][jb], Mach[i][jb]);
			U1[i][jb]  = Rho[i][jb];
			U2[i][jb]  = Rho[i][jb] * u[i][jb];
			U3[i][jb]  = Rho[i][jb] * v[i][jb];
			U4[i][jb]  = Rho[i][jb] * E[i][jb];
			Un1[i][jb] = U1[i][jb];
			Un2[i][jb] = U2[i][jb];
			Un3[i][jb] = U3[i][jb];
			Un4[i][jb] = U4[i][jb];
		}
	}

	// Immersed boundary ghost cells
	const int own_i_beg = std::max(2, local_first_i);
	const int own_i_end = std::min(Nx, local_last_i);
	for (j = 2; j <= Ny; j++)
	{
		for (i = own_i_beg; i <= own_i_end; i++)
		{
			if (Gst[i][j] == 1)
			{
				im = ImBx[i][j] + 1;
				jm = ImBy[i][j] + 1;
				const double Px = ImX[i][j];
				const double Py = ImY[i][j];
				const double X1 = Cx[im][jm];
				const double Y1 = Cy[im][jm];
				const double X2 = Cx[im + 1][jm];
				const double Y2 = Cy[im + 1][jm];
				const double X3 = Cx[im][jm + 1];
				const double Y3 = Cy[im][jm + 1];
				const double X4 = Cx[im + 1][jm + 1];
				const double Y4 = Cy[im + 1][jm + 1];

				Rho[i][j] = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
				                     Rho[im][jm], Rho[im + 1][jm], Rho[im][jm + 1], Rho[im + 1][jm + 1]);
				const double um = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
				                           u[im][jm], u[im + 1][jm], u[im][jm + 1], u[im + 1][jm + 1]);
				const double vm = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
				                           v[im][jm], v[im + 1][jm], v[im][jm + 1], v[im + 1][jm + 1]);
				const double Tim = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
				                            T[im][jm], T[im + 1][jm], T[im][jm + 1], T[im + 1][jm + 1]);

#if USE_NAVIER_STOKES
				ibmWallVelocity(um, vm, ImNx[i][j], ImNy[i][j], Rho[i][j], LAMINAR_VISCOSITY,
				                Cx[i][j], Cy[i][j], Geomx[i][j], Geomy[i][j], Rad[i][j],
				                u[i][j], v[i][j]);
				T[i][j] = ibmWallTemperature(Tim);
				ibmCloseWallThermodynamics(Rho[i][j], u[i][j], v[i][j], T[i][j],
				                           GAS_CONSTANT, GAMMA, GAMMA_MINUS_1,
				                           P[i][j], E[i][j], a[i][j], H[i][j], Mach[i][j]);
#else
				ibmSlipVelocity(um, vm, nx_, ny_, u[i][j], v[i][j]);
				P[i][j] = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
				                   P[im][jm], P[im + 1][jm], P[im][jm + 1], P[im + 1][jm + 1]);
				updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
				                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
#endif
			}
		}
	}
}

/*****************************************************************************************/
/** INTERFACE FLUX CALCULATION (MPI-partitioned + OpenMP) */
/*****************************************************************************************/
void calculateInterfaceFluxes()
{
	const int i_beg = std::max(2, local_first_i);
	const int i_end = std::min(Nx, local_last_i);
#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] == 1)
			{
				const double rho_ij = Rho[i][j];
				const double u_ij = u[i][j];
				const double v_ij = v[i][j];
				const double p_ij = P[i][j];
				const double h_ij = H[i][j];

				const double rho_ip1j = Rho[i + 1][j];
				const double u_ip1j = u[i + 1][j];
				const double v_ip1j = v[i + 1][j];
				const double p_ip1j = P[i + 1][j];
				const double h_ip1j = H[i + 1][j];
				const double nx_r = nxr[i][j];
				const double ny_r = nyr[i][j];

#if USE_WENO5
				const bool ok_r = (i >= 4 && i <= Nx - 2);
				if (ok_r)
				{
					const double rLm = weno5_q_minus(Rho[i - 2][j], Rho[i - 1][j], Rho[i][j], Rho[i + 1][j], Rho[i + 2][j]);
					const double rRm = weno5_q_plus(Rho[i - 1][j], Rho[i][j], Rho[i + 1][j], Rho[i + 2][j], Rho[i + 3][j]);
					const double uLm = weno5_q_minus(u[i - 2][j], u[i - 1][j], u[i][j], u[i + 1][j], u[i + 2][j]);
					const double uRm = weno5_q_plus(u[i - 1][j], u[i][j], u[i + 1][j], u[i + 2][j], u[i + 3][j]);
					const double vLm = weno5_q_minus(v[i - 2][j], v[i - 1][j], v[i][j], v[i + 1][j], v[i + 2][j]);
					const double vRm = weno5_q_plus(v[i - 1][j], v[i][j], v[i + 1][j], v[i + 2][j], v[i + 3][j]);
					const double pLm = weno5_q_minus(P[i - 2][j], P[i - 1][j], P[i][j], P[i + 1][j], P[i + 2][j]);
					const double pRm = weno5_q_plus(P[i - 1][j], P[i][j], P[i + 1][j], P[i + 2][j], P[i + 3][j]);
					double EL, ER, aL, aR, HL, HR, TL, TR, ML, MR;
					updateAuxiliary(rLm, uLm, vLm, pLm, EL, aL, HL, TL, ML);
					updateAuxiliary(rRm, uRm, vRm, pRm, ER, aR, HR, TR, MR);
					double f1a, f2a, f3a, f4a, f1b, f2b, f3b, f4b;
					eulerFluxDotN(rLm, uLm, vLm, pLm, HL, nx_r, ny_r, f1a, f2a, f3a, f4a);
					eulerFluxDotN(rRm, uRm, vRm, pRm, HR, nx_r, ny_r, f1b, f2b, f3b, f4b);
					IF1R[i][j] = 0.5 * (f1a + f1b);
					IF2R[i][j] = 0.5 * (f2a + f2b);
					IF3R[i][j] = 0.5 * (f3a + f3b);
					IF4R[i][j] = 0.5 * (f4a + f4b);
				}
				else
				{
					const double rho_u_sum_r = rho_ip1j * u_ip1j + rho_ij * u_ij;
					const double rho_v_sum_r = rho_ip1j * v_ip1j + rho_ij * v_ij;
					const double rho_uv_sum_r = rho_ip1j * u_ip1j * v_ip1j + rho_ij * u_ij * v_ij;
					IF1R[i][j] = 0.5 * (rho_u_sum_r * nx_r + rho_v_sum_r * ny_r);
					IF2R[i][j] = 0.5 * ((p_ip1j + rho_ip1j * u_ip1j * u_ip1j + p_ij + rho_ij * u_ij * u_ij) * nx_r + rho_uv_sum_r * ny_r);
					IF3R[i][j] = 0.5 * (rho_uv_sum_r * nx_r + (p_ip1j + rho_ip1j * v_ip1j * v_ip1j + p_ij + rho_ij * v_ij * v_ij) * ny_r);
					IF4R[i][j] = 0.5 * ((rho_ip1j * u_ip1j * h_ip1j + rho_ij * u_ij * h_ij) * nx_r + (rho_ip1j * v_ip1j * h_ip1j + rho_ij * v_ij * h_ij) * ny_r);
				}
#else
				const double rho_u_sum_r = rho_ip1j * u_ip1j + rho_ij * u_ij;
				const double rho_v_sum_r = rho_ip1j * v_ip1j + rho_ij * v_ij;
				const double rho_uv_sum_r = rho_ip1j * u_ip1j * v_ip1j + rho_ij * u_ij * v_ij;
				IF1R[i][j] = 0.5 * (rho_u_sum_r * nx_r + rho_v_sum_r * ny_r);
				IF2R[i][j] = 0.5 * ((p_ip1j + rho_ip1j * u_ip1j * u_ip1j + p_ij + rho_ij * u_ij * u_ij) * nx_r + rho_uv_sum_r * ny_r);
				IF3R[i][j] = 0.5 * (rho_uv_sum_r * nx_r + (p_ip1j + rho_ip1j * v_ip1j * v_ip1j + p_ij + rho_ij * v_ij * v_ij) * ny_r);
				IF4R[i][j] = 0.5 * ((rho_ip1j * u_ip1j * h_ip1j + rho_ij * u_ij * h_ij) * nx_r + (rho_ip1j * v_ip1j * h_ip1j + rho_ij * v_ij * h_ij) * ny_r);
#endif

				const double rho_ijp1 = Rho[i][j + 1];
				const double u_ijp1 = u[i][j + 1];
				const double v_ijp1 = v[i][j + 1];
				const double p_ijp1 = P[i][j + 1];
				const double h_ijp1 = H[i][j + 1];
				const double nx_t = nxt[i][j];
				const double ny_t = nyt[i][j];

#if USE_WENO5
				const bool ok_t = (j >= 4 && j <= Ny - 2);
				if (ok_t)
				{
					const double rLm = weno5_q_minus(Rho[i][j - 2], Rho[i][j - 1], Rho[i][j], Rho[i][j + 1], Rho[i][j + 2]);
					const double rRm = weno5_q_plus(Rho[i][j - 1], Rho[i][j], Rho[i][j + 1], Rho[i][j + 2], Rho[i][j + 3]);
					const double uLm = weno5_q_minus(u[i][j - 2], u[i][j - 1], u[i][j], u[i][j + 1], u[i][j + 2]);
					const double uRm = weno5_q_plus(u[i][j - 1], u[i][j], u[i][j + 1], u[i][j + 2], u[i][j + 3]);
					const double vLm = weno5_q_minus(v[i][j - 2], v[i][j - 1], v[i][j], v[i][j + 1], v[i][j + 2]);
					const double vRm = weno5_q_plus(v[i][j - 1], v[i][j], v[i][j + 1], v[i][j + 2], v[i][j + 3]);
					const double pLm = weno5_q_minus(P[i][j - 2], P[i][j - 1], P[i][j], P[i][j + 1], P[i][j + 2]);
					const double pRm = weno5_q_plus(P[i][j - 1], P[i][j], P[i][j + 1], P[i][j + 2], P[i][j + 3]);
					double EL, ER, aL, aR, HL, HR, TL, TR, ML, MR;
					updateAuxiliary(rLm, uLm, vLm, pLm, EL, aL, HL, TL, ML);
					updateAuxiliary(rRm, uRm, vRm, pRm, ER, aR, HR, TR, MR);
					double f1a, f2a, f3a, f4a, f1b, f2b, f3b, f4b;
					eulerFluxDotN(rLm, uLm, vLm, pLm, HL, nx_t, ny_t, f1a, f2a, f3a, f4a);
					eulerFluxDotN(rRm, uRm, vRm, pRm, HR, nx_t, ny_t, f1b, f2b, f3b, f4b);
					IF1T[i][j] = 0.5 * (f1a + f1b);
					IF2T[i][j] = 0.5 * (f2a + f2b);
					IF3T[i][j] = 0.5 * (f3a + f3b);
					IF4T[i][j] = 0.5 * (f4a + f4b);
				}
				else
				{
					const double rho_u_sum_t = rho_ijp1 * u_ijp1 + rho_ij * u_ij;
					const double rho_v_sum_t = rho_ijp1 * v_ijp1 + rho_ij * v_ij;
					const double rho_uv_sum_t = rho_ijp1 * u_ijp1 * v_ijp1 + rho_ij * u_ij * v_ij;
					IF1T[i][j] = 0.5 * (rho_u_sum_t * nx_t + rho_v_sum_t * ny_t);
					IF2T[i][j] = 0.5 * ((p_ijp1 + rho_ijp1 * u_ijp1 * u_ijp1 + p_ij + rho_ij * u_ij * u_ij) * nx_t + rho_uv_sum_t * ny_t);
					IF3T[i][j] = 0.5 * (rho_uv_sum_t * nx_t + (p_ijp1 + rho_ijp1 * v_ijp1 * v_ijp1 + p_ij + rho_ij * v_ij * v_ij) * ny_t);
					IF4T[i][j] = 0.5 * ((rho_ijp1 * u_ijp1 * h_ijp1 + rho_ij * u_ij * h_ij) * nx_t + (rho_ijp1 * v_ijp1 * h_ijp1 + rho_ij * v_ij * h_ij) * ny_t);
				}
#else
				const double rho_u_sum_t = rho_ijp1 * u_ijp1 + rho_ij * u_ij;
				const double rho_v_sum_t = rho_ijp1 * v_ijp1 + rho_ij * v_ij;
				const double rho_uv_sum_t = rho_ijp1 * u_ijp1 * v_ijp1 + rho_ij * u_ij * v_ij;
				IF1T[i][j] = 0.5 * (rho_u_sum_t * nx_t + rho_v_sum_t * ny_t);
				IF2T[i][j] = 0.5 * ((p_ijp1 + rho_ijp1 * u_ijp1 * u_ijp1 + p_ij + rho_ij * u_ij * u_ij) * nx_t + rho_uv_sum_t * ny_t);
				IF3T[i][j] = 0.5 * (rho_uv_sum_t * nx_t + (p_ijp1 + rho_ijp1 * v_ijp1 * v_ijp1 + p_ij + rho_ij * v_ij * v_ij) * ny_t);
				IF4T[i][j] = 0.5 * ((rho_ijp1 * u_ijp1 * h_ijp1 + rho_ij * u_ij * h_ij) * nx_t + (rho_ijp1 * v_ijp1 * h_ijp1 + rho_ij * v_ij * h_ij) * ny_t);
#endif

				const double rho_im1j = Rho[i - 1][j];
				const double u_im1j = u[i - 1][j];
				const double v_im1j = v[i - 1][j];
				const double p_im1j = P[i - 1][j];
				const double h_im1j = H[i - 1][j];
				const double nx_l = nxl[i][j];
				const double ny_l = nyl[i][j];

#if USE_WENO5
				const bool ok_l = (i >= 5 && i <= Nx - 1);
				if (ok_l)
				{
					const double rLm = weno5_q_minus(Rho[i - 3][j], Rho[i - 2][j], Rho[i - 1][j], Rho[i][j], Rho[i + 1][j]);
					const double rRm = weno5_q_plus(Rho[i - 4][j], Rho[i - 3][j], Rho[i - 2][j], Rho[i - 1][j], Rho[i][j]);
					const double uLm = weno5_q_minus(u[i - 3][j], u[i - 2][j], u[i - 1][j], u[i][j], u[i + 1][j]);
					const double uRm = weno5_q_plus(u[i - 4][j], u[i - 3][j], u[i - 2][j], u[i - 1][j], u[i][j]);
					const double vLm = weno5_q_minus(v[i - 3][j], v[i - 2][j], v[i - 1][j], v[i][j], v[i + 1][j]);
					const double vRm = weno5_q_plus(v[i - 4][j], v[i - 3][j], v[i - 2][j], v[i - 1][j], v[i][j]);
					const double pLm = weno5_q_minus(P[i - 3][j], P[i - 2][j], P[i - 1][j], P[i][j], P[i + 1][j]);
					const double pRm = weno5_q_plus(P[i - 4][j], P[i - 3][j], P[i - 2][j], P[i - 1][j], P[i][j]);
					double EL, ER, aL, aR, HL, HR, TL, TR, ML, MR;
					updateAuxiliary(rLm, uLm, vLm, pLm, EL, aL, HL, TL, ML);
					updateAuxiliary(rRm, uRm, vRm, pRm, ER, aR, HR, TR, MR);
					double f1a, f2a, f3a, f4a, f1b, f2b, f3b, f4b;
					eulerFluxDotN(rLm, uLm, vLm, pLm, HL, nx_l, ny_l, f1a, f2a, f3a, f4a);
					eulerFluxDotN(rRm, uRm, vRm, pRm, HR, nx_l, ny_l, f1b, f2b, f3b, f4b);
					IF1L[i][j] = 0.5 * (f1a + f1b);
					IF2L[i][j] = 0.5 * (f2a + f2b);
					IF3L[i][j] = 0.5 * (f3a + f3b);
					IF4L[i][j] = 0.5 * (f4a + f4b);
				}
				else
				{
					const double rho_u_sum_l = rho_im1j * u_im1j + rho_ij * u_ij;
					const double rho_v_sum_l = rho_im1j * v_im1j + rho_ij * v_ij;
					const double rho_uv_sum_l = rho_im1j * u_im1j * v_im1j + rho_ij * u_ij * v_ij;
					IF1L[i][j] = 0.5 * (rho_u_sum_l * nx_l + rho_v_sum_l * ny_l);
					IF2L[i][j] = 0.5 * ((p_im1j + rho_im1j * u_im1j * u_im1j + p_ij + rho_ij * u_ij * u_ij) * nx_l + rho_uv_sum_l * ny_l);
					IF3L[i][j] = 0.5 * (rho_uv_sum_l * nx_l + (p_im1j + rho_im1j * v_im1j * v_im1j + p_ij + rho_ij * v_ij * v_ij) * ny_l);
					IF4L[i][j] = 0.5 * ((rho_im1j * u_im1j * h_im1j + rho_ij * u_ij * h_ij) * nx_l + (rho_im1j * v_im1j * h_im1j + rho_ij * v_ij * h_ij) * ny_l);
				}
#else
				const double rho_u_sum_l = rho_im1j * u_im1j + rho_ij * u_ij;
				const double rho_v_sum_l = rho_im1j * v_im1j + rho_ij * v_ij;
				const double rho_uv_sum_l = rho_im1j * u_im1j * v_im1j + rho_ij * u_ij * v_ij;
				IF1L[i][j] = 0.5 * (rho_u_sum_l * nx_l + rho_v_sum_l * ny_l);
				IF2L[i][j] = 0.5 * ((p_im1j + rho_im1j * u_im1j * u_im1j + p_ij + rho_ij * u_ij * u_ij) * nx_l + rho_uv_sum_l * ny_l);
				IF3L[i][j] = 0.5 * (rho_uv_sum_l * nx_l + (p_im1j + rho_im1j * v_im1j * v_im1j + p_ij + rho_ij * v_ij * v_ij) * ny_l);
				IF4L[i][j] = 0.5 * ((rho_im1j * u_im1j * h_im1j + rho_ij * u_ij * h_ij) * nx_l + (rho_im1j * v_im1j * h_im1j + rho_ij * v_ij * h_ij) * ny_l);
#endif

				const double rho_ijm1 = Rho[i][j - 1];
				const double u_ijm1 = u[i][j - 1];
				const double v_ijm1 = v[i][j - 1];
				const double p_ijm1 = P[i][j - 1];
				const double h_ijm1 = H[i][j - 1];
				const double nx_b = nxb[i][j];
				const double ny_b = nyb[i][j];

#if USE_WENO5
				const bool ok_b = (j >= 5 && j <= Ny - 1);
				if (ok_b)
				{
					const double rLm = weno5_q_minus(Rho[i][j - 3], Rho[i][j - 2], Rho[i][j - 1], Rho[i][j], Rho[i][j + 1]);
					const double rRm = weno5_q_plus(Rho[i][j - 4], Rho[i][j - 3], Rho[i][j - 2], Rho[i][j - 1], Rho[i][j]);
					const double uLm = weno5_q_minus(u[i][j - 3], u[i][j - 2], u[i][j - 1], u[i][j], u[i][j + 1]);
					const double uRm = weno5_q_plus(u[i][j - 4], u[i][j - 3], u[i][j - 2], u[i][j - 1], u[i][j]);
					const double vLm = weno5_q_minus(v[i][j - 3], v[i][j - 2], v[i][j - 1], v[i][j], v[i][j + 1]);
					const double vRm = weno5_q_plus(v[i][j - 4], v[i][j - 3], v[i][j - 2], v[i][j - 1], v[i][j]);
					const double pLm = weno5_q_minus(P[i][j - 3], P[i][j - 2], P[i][j - 1], P[i][j], P[i][j + 1]);
					const double pRm = weno5_q_plus(P[i][j - 4], P[i][j - 3], P[i][j - 2], P[i][j - 1], P[i][j]);
					double EL, ER, aL, aR, HL, HR, TL, TR, ML, MR;
					updateAuxiliary(rLm, uLm, vLm, pLm, EL, aL, HL, TL, ML);
					updateAuxiliary(rRm, uRm, vRm, pRm, ER, aR, HR, TR, MR);
					double f1a, f2a, f3a, f4a, f1b, f2b, f3b, f4b;
					eulerFluxDotN(rLm, uLm, vLm, pLm, HL, nx_b, ny_b, f1a, f2a, f3a, f4a);
					eulerFluxDotN(rRm, uRm, vRm, pRm, HR, nx_b, ny_b, f1b, f2b, f3b, f4b);
					IF1B[i][j] = 0.5 * (f1a + f1b);
					IF2B[i][j] = 0.5 * (f2a + f2b);
					IF3B[i][j] = 0.5 * (f3a + f3b);
					IF4B[i][j] = 0.5 * (f4a + f4b);
				}
				else
				{
					const double rho_u_sum_b = rho_ijm1 * u_ijm1 + rho_ij * u_ij;
					const double rho_v_sum_b = rho_ijm1 * v_ijm1 + rho_ij * v_ij;
					const double rho_uv_sum_b = rho_ijm1 * u_ijm1 * v_ijm1 + rho_ij * u_ij * v_ij;
					IF1B[i][j] = 0.5 * (rho_u_sum_b * nx_b + rho_v_sum_b * ny_b);
					IF2B[i][j] = 0.5 * ((p_ijm1 + rho_ijm1 * u_ijm1 * u_ijm1 + p_ij + rho_ij * u_ij * u_ij) * nx_b + rho_uv_sum_b * ny_b);
					IF3B[i][j] = 0.5 * (rho_uv_sum_b * nx_b + (p_ijm1 + rho_ijm1 * v_ijm1 * v_ijm1 + p_ij + rho_ij * v_ij * v_ij) * ny_b);
					IF4B[i][j] = 0.5 * ((rho_ijm1 * u_ijm1 * h_ijm1 + rho_ij * u_ij * h_ij) * nx_b + (rho_ijm1 * v_ijm1 * h_ijm1 + rho_ij * v_ij * h_ij) * ny_b);
				}
#else
				const double rho_u_sum_b = rho_ijm1 * u_ijm1 + rho_ij * u_ij;
				const double rho_v_sum_b = rho_ijm1 * v_ijm1 + rho_ij * v_ij;
				const double rho_uv_sum_b = rho_ijm1 * u_ijm1 * v_ijm1 + rho_ij * u_ij * v_ij;
				IF1B[i][j] = 0.5 * (rho_u_sum_b * nx_b + rho_v_sum_b * ny_b);
				IF2B[i][j] = 0.5 * ((p_ijm1 + rho_ijm1 * u_ijm1 * u_ijm1 + p_ij + rho_ij * u_ij * u_ij) * nx_b + rho_uv_sum_b * ny_b);
				IF3B[i][j] = 0.5 * (rho_uv_sum_b * nx_b + (p_ijm1 + rho_ijm1 * v_ijm1 * v_ijm1 + p_ij + rho_ij * v_ij * v_ij) * ny_b);
				IF4B[i][j] = 0.5 * ((rho_ijm1 * u_ijm1 * h_ijm1 + rho_ij * u_ij * h_ij) * nx_b + (rho_ijm1 * v_ijm1 * h_ijm1 + rho_ij * v_ij * h_ij) * ny_b);
#endif
			}
		}
	}
}

/*****************************************************************************************/
/** DISSIPATIVE FLUX CALCULATION (MPI-partitioned + OpenMP) */
/*****************************************************************************************/
void calculateDissipativeFluxes()
{
	const int i_beg = std::max(2, local_first_i);
	const int i_end = std::min(Nx, local_last_i);
#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static) \
	    private(amax, amin, alpha1, alpha2, alpha3, alpha4)
#endif
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] == 1)
			{
				// Right interface dissipative flux
				amax = MaxEigVal(u[i + 1][j], u[i][j], v[i + 1][j], v[i][j], a[i + 1][j], a[i][j], nxr[i][j], nyr[i][j]);
				amin = MinEigVal(u[i + 1][j], u[i][j], v[i + 1][j], v[i][j], a[i + 1][j], a[i][j], nxr[i][j], nyr[i][j]);

				alpha1 = MoversN(Rho[i + 1][j] * u[i + 1][j] * nxr[i][j] + Rho[i + 1][j] * v[i + 1][j] * nyr[i][j],
								 Rho[i][j] * u[i][j] * nxr[i][j] + Rho[i][j] * v[i][j] * nyr[i][j],
								 Rho[i + 1][j], Rho[i][j], amax, amin);
				alpha2 = MoversN((Rho[i + 1][j] * u[i + 1][j] * u[i + 1][j] + P[i + 1][j]) * nxr[i][j] + (Rho[i + 1][j] * u[i + 1][j] * v[i + 1][j]) * nyr[i][j],
								 (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxr[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyr[i][j],
								 Rho[i + 1][j] * u[i + 1][j], Rho[i][j] * u[i][j], amax, amin);
				alpha3 = MoversN((Rho[i + 1][j] * v[i + 1][j] * v[i + 1][j] + P[i + 1][j]) * nyr[i][j] + (Rho[i + 1][j] * u[i + 1][j] * v[i + 1][j]) * nxr[i][j],
								 (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyr[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxr[i][j],
								 Rho[i + 1][j] * v[i + 1][j], Rho[i][j] * v[i][j], amax, amin);
				alpha4 = MoversN(Rho[i + 1][j] * u[i + 1][j] * H[i + 1][j] * nxr[i][j] + Rho[i + 1][j] * v[i + 1][j] * H[i + 1][j] * nyr[i][j],
								 Rho[i][j] * u[i][j] * H[i][j] * nxr[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyr[i][j],
								 Rho[i + 1][j] * E[i + 1][j], Rho[i][j] * E[i][j], amax, amin);

				DF1R[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i + 1][j]);
				DF2R[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i + 1][j] * u[i + 1][j]);
				DF3R[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i + 1][j] * v[i + 1][j]);
				DF4R[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i + 1][j] * E[i + 1][j]);

				// Top interface dissipative flux
				amax = MaxEigVal(u[i][j + 1], u[i][j], v[i][j + 1], v[i][j], a[i][j + 1], a[i][j], nxt[i][j], nyt[i][j]);
				amin = MinEigVal(u[i][j + 1], u[i][j], v[i][j + 1], v[i][j], a[i][j + 1], a[i][j], nxt[i][j], nyt[i][j]);

				alpha1 = MoversN(Rho[i][j + 1] * u[i][j + 1] * nxt[i][j] + Rho[i][j + 1] * v[i][j + 1] * nyt[i][j],
								 Rho[i][j] * u[i][j] * nxt[i][j] + Rho[i][j] * v[i][j] * nyt[i][j],
								 Rho[i][j + 1], Rho[i][j], amax, amin);
				alpha2 = MoversN((Rho[i][j + 1] * u[i][j + 1] * u[i][j + 1] + P[i][j + 1]) * nxt[i][j] + (Rho[i][j + 1] * u[i][j + 1] * v[i][j + 1]) * nyt[i][j],
								 (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxt[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyt[i][j],
								 Rho[i][j + 1] * u[i][j + 1], Rho[i][j] * u[i][j], amax, amin);
				alpha3 = MoversN((Rho[i][j + 1] * v[i][j + 1] * v[i][j + 1] + P[i][j + 1]) * nyt[i][j] + (Rho[i][j + 1] * u[i][j + 1] * v[i][j + 1]) * nxt[i][j],
								 (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyt[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxt[i][j],
								 Rho[i][j + 1] * v[i][j + 1], Rho[i][j] * v[i][j], amax, amin);
				alpha4 = MoversN(Rho[i][j + 1] * u[i][j + 1] * H[i][j + 1] * nxt[i][j] + Rho[i][j + 1] * v[i][j + 1] * H[i][j + 1] * nyt[i][j],
								 Rho[i][j] * u[i][j] * H[i][j] * nxt[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyt[i][j],
								 Rho[i][j + 1] * E[i][j + 1], Rho[i][j] * E[i][j], amax, amin);

				DF1T[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i][j + 1]);
				DF2T[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i][j + 1] * u[i][j + 1]);
				DF3T[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i][j + 1] * v[i][j + 1]);
				DF4T[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i][j + 1] * E[i][j + 1]);

				// Left interface dissipative flux
				amax = MaxEigVal(u[i][j], u[i - 1][j], v[i][j], v[i - 1][j], a[i][j], a[i - 1][j], nxl[i][j], nyl[i][j]);
				amin = MinEigVal(u[i][j], u[i - 1][j], v[i][j], v[i - 1][j], a[i][j], a[i - 1][j], nxl[i][j], nyl[i][j]);

				alpha1 = MoversN(Rho[i - 1][j] * u[i - 1][j] * nxl[i][j] + Rho[i - 1][j] * v[i - 1][j] * nyl[i][j],
								 Rho[i][j] * u[i][j] * nxl[i][j] + Rho[i][j] * v[i][j] * nyl[i][j],
								 Rho[i - 1][j], Rho[i][j], amax, amin);
				alpha2 = MoversN((Rho[i - 1][j] * u[i - 1][j] * u[i - 1][j] + P[i - 1][j]) * nxl[i][j] + (Rho[i - 1][j] * u[i - 1][j] * v[i - 1][j]) * nyl[i][j],
								 (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxl[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyl[i][j],
								 Rho[i - 1][j] * u[i - 1][j], Rho[i][j] * u[i][j], amax, amin);
				alpha3 = MoversN((Rho[i - 1][j] * v[i - 1][j] * v[i - 1][j] + P[i - 1][j]) * nyl[i][j] + (Rho[i - 1][j] * u[i - 1][j] * v[i - 1][j]) * nxl[i][j],
								 (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyl[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxl[i][j],
								 Rho[i - 1][j] * v[i - 1][j], Rho[i][j] * v[i][j], amax, amin);
				alpha4 = MoversN(Rho[i - 1][j] * u[i - 1][j] * H[i - 1][j] * nxl[i][j] + Rho[i - 1][j] * v[i - 1][j] * H[i - 1][j] * nyl[i][j],
								 Rho[i][j] * u[i][j] * H[i][j] * nxl[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyl[i][j],
								 Rho[i - 1][j] * E[i - 1][j], Rho[i][j] * E[i][j], amax, amin);

				DF1L[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i - 1][j]);
				DF2L[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i - 1][j] * u[i - 1][j]);
				DF3L[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i - 1][j] * v[i - 1][j]);
				DF4L[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i - 1][j] * E[i - 1][j]);

				// Bottom interface dissipative flux
				amax = MaxEigVal(u[i][j], u[i][j - 1], v[i][j], v[i][j - 1], a[i][j], a[i][j - 1], nxb[i][j], nyb[i][j]);
				amin = MinEigVal(u[i][j], u[i][j - 1], v[i][j], v[i][j - 1], a[i][j], a[i][j - 1], nxb[i][j], nyb[i][j]);

				alpha1 = MoversN(Rho[i][j - 1] * u[i][j - 1] * nxb[i][j] + Rho[i][j - 1] * v[i][j - 1] * nyb[i][j],
								 Rho[i][j] * u[i][j] * nxb[i][j] + Rho[i][j] * v[i][j] * nyb[i][j],
								 Rho[i][j - 1], Rho[i][j], amax, amin);
				alpha2 = MoversN((Rho[i][j - 1] * u[i][j - 1] * u[i][j - 1] + P[i][j - 1]) * nxb[i][j] + (Rho[i][j - 1] * u[i][j - 1] * v[i][j - 1]) * nyb[i][j],
								 (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxb[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyb[i][j],
								 Rho[i][j - 1] * u[i][j - 1], Rho[i][j] * u[i][j], amax, amin);
				alpha3 = MoversN((Rho[i][j - 1] * v[i][j - 1] * v[i][j - 1] + P[i][j - 1]) * nyb[i][j] + (Rho[i][j - 1] * u[i][j - 1] * v[i][j - 1]) * nxb[i][j],
								 (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyb[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxb[i][j],
								 Rho[i][j - 1] * v[i][j - 1], Rho[i][j] * v[i][j], amax, amin);
				alpha4 = MoversN(Rho[i][j - 1] * u[i][j - 1] * H[i][j - 1] * nxb[i][j] + Rho[i][j - 1] * v[i][j - 1] * H[i][j - 1] * nyb[i][j],
								 Rho[i][j] * u[i][j] * H[i][j] * nxb[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyb[i][j],
								 Rho[i][j - 1] * E[i][j - 1], Rho[i][j] * E[i][j], amax, amin);

				DF1B[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i][j - 1]);
				DF2B[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i][j - 1] * u[i][j - 1]);
				DF3B[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i][j - 1] * v[i][j - 1]);
				DF4B[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i][j - 1] * E[i][j - 1]);
			}
		}
	}
}

#if USE_NAVIER_STOKES
/*****************************************************************************************/
/** Laminar Newtonian viscous flux (face values; added to residual like IF/DF).           */
/*****************************************************************************************/
static inline void metric_grad_at(int i, int j,
                                double &ux, double &uy, double &vx, double &vy,
                                double &Tx, double &Ty)
{
	const double x_xi = 0.5 * (Cx[i + 1][j] - Cx[i - 1][j]);
	const double x_eta = 0.5 * (Cx[i][j + 1] - Cx[i][j - 1]);
	const double y_xi = 0.5 * (Cy[i + 1][j] - Cy[i - 1][j]);
	const double y_eta = 0.5 * (Cy[i][j + 1] - Cy[i][j - 1]);
	const double J = x_xi * y_eta - x_eta * y_xi;
	const double invJ = 1.0 / (std::fabs(J) + 1.0e-30);

	const double u_xi = 0.5 * (u[i + 1][j] - u[i - 1][j]);
	const double u_eta = 0.5 * (u[i][j + 1] - u[i][j - 1]);
	ux = invJ * (u_xi * y_eta - u_eta * y_xi);
	uy = invJ * (u_eta * x_xi - u_xi * x_eta);

	const double v_xi = 0.5 * (v[i + 1][j] - v[i - 1][j]);
	const double v_eta = 0.5 * (v[i][j + 1] - v[i][j - 1]);
	vx = invJ * (v_xi * y_eta - v_eta * y_xi);
	vy = invJ * (v_eta * x_xi - v_xi * x_eta);

	const double T_xi = 0.5 * (T[i + 1][j] - T[i - 1][j]);
	const double T_eta = 0.5 * (T[i][j + 1] - T[i][j - 1]);
	Tx = invJ * (T_xi * y_eta - T_eta * y_xi);
	Ty = invJ * (T_eta * x_xi - T_xi * x_eta);
}

void calculateViscousFluxes()
{
	const double mu = LAMINAR_VISCOSITY;
	const double kappa = thermalConductivity(mu, GAMMA, GAS_CONSTANT, GAMMA_MINUS_1, PRANDTL_NUMBER);

	const int z_i0 = std::max(2, local_first_i - 1);
	const int z_i1 = std::min(Nx + 1, local_last_i + 1);
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = z_i0; i <= z_i1; i++)
		{
			VF1R[i][j] = VF2R[i][j] = VF3R[i][j] = VF4R[i][j] = 0.0;
			VF1T[i][j] = VF2T[i][j] = VF3T[i][j] = VF4T[i][j] = 0.0;
			VF1L[i][j] = VF2L[i][j] = VF3L[i][j] = VF4L[i][j] = 0.0;
			VF1B[i][j] = VF2B[i][j] = VF3B[i][j] = VF4B[i][j] = 0.0;
		}
	}

	const int i_beg = std::max(3, local_first_i);
	const int i_end = std::min(Nx - 1, local_last_i);
#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 3; j <= Ny; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] != 1)
				continue;

			double ux1, uy1, vx1, vy1, Tx1, Ty1;
			double ux2, uy2, vx2, vy2, Tx2, Ty2;
			double f1, f2, f3, f4;

			if (viscousFaceActive(Fl[i][j], Gst[i][j], Fl[i + 1][j], Gst[i + 1][j]))
			{
				metric_grad_at(i, j, ux1, uy1, vx1, vy1, Tx1, Ty1);
				metric_grad_at(i + 1, j, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2);
				const double uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2);
				const double vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2);
				const double Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i][j] + u[i + 1][j]);
				const double vc = 0.5 * (v[i][j] + v[i + 1][j]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxr[i][j], nyr[i][j], f1, f2, f3, f4);
				VF1R[i][j] = f1;
				VF2R[i][j] = f2;
				VF3R[i][j] = f3;
				VF4R[i][j] = f4;
			}

			if (viscousFaceActive(Fl[i][j], Gst[i][j], Fl[i][j + 1], Gst[i][j + 1]))
			{
				metric_grad_at(i, j, ux1, uy1, vx1, vy1, Tx1, Ty1);
				metric_grad_at(i, j + 1, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2);
				const double uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2);
				const double vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2);
				const double Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i][j] + u[i][j + 1]);
				const double vc = 0.5 * (v[i][j] + v[i][j + 1]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxt[i][j], nyt[i][j], f1, f2, f3, f4);
				VF1T[i][j] = f1;
				VF2T[i][j] = f2;
				VF3T[i][j] = f3;
				VF4T[i][j] = f4;
			}

			if (viscousFaceActive(Fl[i - 1][j], Gst[i - 1][j], Fl[i][j], Gst[i][j]))
			{
				metric_grad_at(i - 1, j, ux1, uy1, vx1, vy1, Tx1, Ty1);
				metric_grad_at(i, j, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2);
				const double uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2);
				const double vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2);
				const double Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i - 1][j] + u[i][j]);
				const double vc = 0.5 * (v[i - 1][j] + v[i][j]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxl[i][j], nyl[i][j], f1, f2, f3, f4);
				VF1L[i][j] = f1;
				VF2L[i][j] = f2;
				VF3L[i][j] = f3;
				VF4L[i][j] = f4;
			}

			if (viscousFaceActive(Fl[i][j - 1], Gst[i][j - 1], Fl[i][j], Gst[i][j]))
			{
				metric_grad_at(i, j - 1, ux1, uy1, vx1, vy1, Tx1, Ty1);
				metric_grad_at(i, j, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2);
				const double uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2);
				const double vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2);
				const double Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i][j - 1] + u[i][j]);
				const double vc = 0.5 * (v[i][j - 1] + v[i][j]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxb[i][j], nyb[i][j], f1, f2, f3, f4);
				VF1B[i][j] = f1;
				VF2B[i][j] = f2;
				VF3B[i][j] = f3;
				VF4B[i][j] = f4;
			}
		}
	}
}
#endif // USE_NAVIER_STOKES

/*****************************************************************************************/
/** TIME STEP CALCULATION (MPI-partitioned + OpenMP)                                      */
/* Each rank scans its owned cells and computes a local minimum dt; MPI_Allreduce gives   */
/* the global minimum used by every rank.                                                 */
/*****************************************************************************************/
void calculateTimeStep()
{
	double local_dtg = 100.0;
	const int i_beg = std::max(2, local_first_i);
	const int i_end = std::min(Nx, local_last_i);
#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static) reduction(min:local_dtg)
#endif
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] == 1)
			{
				const double cr = AvgEigVal(u[i + 1][j], u[i][j], v[i + 1][j], v[i][j], a[i + 1][j], a[i][j], nxr[i][j], nyr[i][j], dsr[i][j]);
				const double ct = AvgEigVal(u[i][j + 1], u[i][j], v[i][j + 1], v[i][j], a[i][j + 1], a[i][j], nxt[i][j], nyt[i][j], dst[i][j]);
				const double cl = AvgEigVal(u[i][j], u[i - 1][j], v[i][j], v[i - 1][j], a[i][j], a[i - 1][j], nxl[i][j], nyl[i][j], dsl[i][j]);
				const double cb = AvgEigVal(u[i][j], u[i][j - 1], v[i][j], v[i][j - 1], a[i][j], a[i][j - 1], nxb[i][j], nyb[i][j], dsb[i][j]);
				double dt_ij = L * Area[i][j] / Maximum4(cr, ct, cl, cb) * 0.5;
#if USE_NAVIER_STOKES
				const double h_min = Minimum4(dsr[i][j], dst[i][j], dsl[i][j], dsb[i][j]);
				const double dt_vis = viscousDtScale(Rho[i][j], LAMINAR_VISCOSITY, h_min, VISCOUS_CFL_FRACTION);
				if (dt_vis < dt_ij)
					dt_ij = dt_vis;
#endif
				dt[i][j] = dt_ij;
				if (dt_ij < local_dtg)
					local_dtg = dt_ij;
			}
		}
	}

	double global_dtg;
	MPI_Allreduce(&local_dtg, &global_dtg, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	dtg = global_dtg;
}

/*****************************************************************************************/
/** UPDATE CONSERVATIVE VARIABLES (MPI-partitioned + OpenMP)                              */
/*                                                                                        */
/* The previous version performed 32 std::isnan() checks per fluid cell — far more        */
/* expensive than the actual update. NaN/inf produced by the flux kernels indicate a      */
/* real problem in the discretization, so we don't paper over them in the hot loop.       */
/* Numerical stability is preserved with simple positivity floors on density and total    */
/* energy. Validity is verified periodically by isValidSolution() (see main loop).        */
/*****************************************************************************************/
void updateConservativeVariables()
{
	const int i_beg = std::max(2, local_first_i);
	const int i_end = std::min(Nx, local_last_i);
	const double rho_min = 0.1 * DENSITY_INFINITY;
	const double e_min_intrinsic = 0.1 * PRESSURE_INFINITY / GAMMA_MINUS_1;

#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] == 1)
			{
				const double inv_area_dt = dtg / Area[i][j];
				const double dsr_ = dsr[i][j], dst_ = dst[i][j], dsl_ = dsl[i][j], dsb_ = dsb[i][j];
#if USE_NAVIER_STOKES
				const double f1 = (IF1R[i][j] + DF1R[i][j] + VF1R[i][j]) * dsr_ + (IF1T[i][j] + DF1T[i][j] + VF1T[i][j]) * dst_ + (IF1L[i][j] + DF1L[i][j] + VF1L[i][j]) * dsl_ + (IF1B[i][j] + DF1B[i][j] + VF1B[i][j]) * dsb_;
				const double f2 = (IF2R[i][j] + DF2R[i][j] + VF2R[i][j]) * dsr_ + (IF2T[i][j] + DF2T[i][j] + VF2T[i][j]) * dst_ + (IF2L[i][j] + DF2L[i][j] + VF2L[i][j]) * dsl_ + (IF2B[i][j] + DF2B[i][j] + VF2B[i][j]) * dsb_;
				const double f3 = (IF3R[i][j] + DF3R[i][j] + VF3R[i][j]) * dsr_ + (IF3T[i][j] + DF3T[i][j] + VF3T[i][j]) * dst_ + (IF3L[i][j] + DF3L[i][j] + VF3L[i][j]) * dsl_ + (IF3B[i][j] + DF3B[i][j] + VF3B[i][j]) * dsb_;
				const double f4 = (IF4R[i][j] + DF4R[i][j] + VF4R[i][j]) * dsr_ + (IF4T[i][j] + DF4T[i][j] + VF4T[i][j]) * dst_ + (IF4L[i][j] + DF4L[i][j] + VF4L[i][j]) * dsl_ + (IF4B[i][j] + DF4B[i][j] + VF4B[i][j]) * dsb_;
#else
				const double f1 = (IF1R[i][j] + DF1R[i][j]) * dsr_ + (IF1T[i][j] + DF1T[i][j]) * dst_ + (IF1L[i][j] + DF1L[i][j]) * dsl_ + (IF1B[i][j] + DF1B[i][j]) * dsb_;
				const double f2 = (IF2R[i][j] + DF2R[i][j]) * dsr_ + (IF2T[i][j] + DF2T[i][j]) * dst_ + (IF2L[i][j] + DF2L[i][j]) * dsl_ + (IF2B[i][j] + DF2B[i][j]) * dsb_;
				const double f3 = (IF3R[i][j] + DF3R[i][j]) * dsr_ + (IF3T[i][j] + DF3T[i][j]) * dst_ + (IF3L[i][j] + DF3L[i][j]) * dsl_ + (IF3B[i][j] + DF3B[i][j]) * dsb_;
				const double f4 = (IF4R[i][j] + DF4R[i][j]) * dsr_ + (IF4T[i][j] + DF4T[i][j]) * dst_ + (IF4L[i][j] + DF4L[i][j]) * dsl_ + (IF4B[i][j] + DF4B[i][j]) * dsb_;
#endif

				double n1 = U1[i][j] - inv_area_dt * f1;
				double n2 = U2[i][j] - inv_area_dt * f2;
				double n3 = U3[i][j] - inv_area_dt * f3;
				double n4 = U4[i][j] - inv_area_dt * f4;

				/* Positivity floors — much cheaper than isnan checks. */
				if (n1 < rho_min) n1 = rho_min;
				const double ke = 0.5 * (n2 * n2 + n3 * n3) / n1;
				const double e_min = ke + e_min_intrinsic;
				if (n4 < e_min) n4 = e_min;

				Un1[i][j] = n1;
				Un2[i][j] = n2;
				Un3[i][j] = n3;
				Un4[i][j] = n4;
			}
		}
	}
}

/*****************************************************************************************/
/** UPDATE PRIMITIVE VARIABLES (MPI-partitioned + OpenMP)                                 */
/*****************************************************************************************/
void updatePrimitiveVariables()
{
	const int i_beg = std::max(2, local_first_i);
	const int i_end = std::min(Nx, local_last_i);
	const double rho_floor = 0.1 * DENSITY_INFINITY;
	const double P_min = 0.01 * PRESSURE_INFINITY;
#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] == 1)
			{
				double rho = Un1[i][j];
				if (rho <= 0.0) { rho = rho_floor; Un1[i][j] = rho; }
				const double inv_rho = 1.0 / rho;
				const double u_ij = Un2[i][j] * inv_rho;
				const double v_ij = Un3[i][j] * inv_rho;
				const double E_ij = Un4[i][j] * inv_rho;
				const double q2   = u_ij * u_ij + v_ij * v_ij;
				double P_ij       = GAMMA_MINUS_1 * rho * (E_ij - 0.5 * q2);
				if (P_ij < P_min) P_ij = P_min;

				Rho[i][j] = rho;
				u[i][j]   = u_ij;
				v[i][j]   = v_ij;
				E[i][j]   = E_ij;
				P[i][j]   = P_ij;
				updateAuxiliary(rho, u_ij, v_ij, P_ij,
				                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);

				U1[i][j] = Un1[i][j];
				U2[i][j] = Un2[i][j];
				U3[i][j] = Un3[i][j];
				U4[i][j] = Un4[i][j];
			}
		}
	}
}

/*****************************************************************************************/
/** SOLUTION VALIDATION (per-rank, owned cells only)                                      */
/*****************************************************************************************/
bool isValidSolution()
{
	const int i_beg = std::max(2, local_first_i);
	const int i_end = std::min(Nx, local_last_i);
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = i_beg; i <= i_end; i++)
		{
			if (Fl[i][j] == 1)
			{
				if (std::isnan(P[i][j]) || std::isnan(Rho[i][j]) || std::isnan(u[i][j]) || std::isnan(v[i][j]) ||
				    P[i][j] <= 0.0 || Rho[i][j] <= 0.0)
				{
					std::cerr << "[rank " << mpi_rank << "] invalid cell (" << i << "," << j
					          << "): rho=" << Rho[i][j] << " P=" << P[i][j]
					          << " u=" << u[i][j] << " v=" << v[i][j] << '\n';
					return false;
				}
			}
		}
	}
	return true;
}

/*****************************************************************************************/
/** gatherFieldToRank0                                                                    */
/*                                                                                        */
/* Gathers a 2D scalar field from each rank's owned i-slab onto the full array on rank 0. */
/* Used before writing output files; other ranks are left with their partial views.       */
/*****************************************************************************************/
static void gatherFieldDouble(double field[MAX_GRID_X][MAX_GRID_Y])
{
	const int Jlen = Ny + 3; /* indices j = 0 .. Ny+2 */
	const int my_slab = local_Nx * Jlen;
	std::vector<double> sendbuf((size_t)my_slab);
	int idx = 0;
	for (int ii = local_first_i; ii <= local_last_i; ++ii)
		for (int jj = 0; jj < Jlen; ++jj)
			sendbuf[idx++] = field[ii][jj];

	std::vector<int> recvcounts, displs;
	std::vector<double> recvbuf;
	const int base = Nx / mpi_size;
	const int rem  = Nx % mpi_size;
	if (mpi_rank == 0) {
		recvcounts.assign(mpi_size, 0);
		displs.assign(mpi_size, 0);
		int disp = 0;
		for (int r = 0; r < mpi_size; ++r) {
			const int local_n = base + (r < rem ? 1 : 0);
			recvcounts[r] = local_n * Jlen;
			displs[r] = disp;
			disp += recvcounts[r];
		}
		recvbuf.resize((size_t)disp);
	}

	MPI_Gatherv(sendbuf.data(), my_slab, MPI_DOUBLE,
	            mpi_rank == 0 ? recvbuf.data() : nullptr,
	            mpi_rank == 0 ? recvcounts.data() : nullptr,
	            mpi_rank == 0 ? displs.data()    : nullptr,
	            MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (mpi_rank == 0) {
		int o = 0;
		for (int r = 0; r < mpi_size; ++r) {
			const int local_n = base + (r < rem ? 1 : 0);
			const int rank_first = 2 + r * base + std::min(r, rem);
			for (int k = 0; k < local_n; ++k)
				for (int jj = 0; jj < Jlen; ++jj)
					field[rank_first + k][jj] = recvbuf[o++];
		}
	}
}

static void gatherDynamicFieldsToRank0()
{
	gatherFieldDouble(Rho);
	gatherFieldDouble(P);
	gatherFieldDouble(u);
	gatherFieldDouble(v);
	gatherFieldDouble(T);
	gatherFieldDouble(Mach);
	gatherFieldDouble(U1);
	gatherFieldDouble(U4);
	gatherFieldDouble(E);
	gatherFieldDouble(H);
}

/*****************************************************************************************/
/** Enhanced VTK Output Function with MPI Data Gathering */
/*****************************************************************************************/
void outputVTK(int timeStep)
{
	gatherDynamicFieldsToRank0();
	// Gather all data to rank 0 for VTK output
	if (mpi_rank == 0)
	{
		string fileName = "Aerospike_IBM_MPI_Solution_" + to_string(Nx) + "x" + to_string(Ny) + "_t" + to_string(timeStep) + ".vtk";
		ofstream vtkFile(fileName);

		if (!vtkFile.is_open())
		{
			cout << "ERROR: Could not open VTK file " << fileName << '\n';
			return;
		}

		// VTK Header
		vtkFile << "# vtk DataFile Version 3.0" << '\n';
		vtkFile << "2D CFD Solution - Aerospike with Immersed Boundary Method - Time Step " << timeStep << '\n';
		vtkFile << "ASCII" << '\n';
		vtkFile << "DATASET STRUCTURED_GRID" << '\n';

		// Grid dimensions (exclude boundary layers as per original output)
		int grid_nx = Nx;
		int grid_ny = Ny - 2;
		vtkFile << "DIMENSIONS " << grid_nx << " " << grid_ny << " 1" << '\n';

		// Total number of points
		int total_points = grid_nx * grid_ny;
		vtkFile << "POINTS " << total_points << " double" << '\n';

		// Write grid coordinates (using cell centers)
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				vtkFile << fixed << setprecision(6) << Cx[i][j] << " " << Cy[i][j] << " 0.0" << '\n';
			}
		}

		// Write point data
		vtkFile << "POINT_DATA " << total_points << '\n';

		// Density field
		vtkFile << "SCALARS Density double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << scientific << setprecision(6) << Rho[i][j] << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n'; // Inside immersed boundary
				}
			}
		}

		// Pressure field
		vtkFile << "SCALARS Pressure double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << scientific << setprecision(6) << P[i][j] << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n';
				}
			}
		}

		// Mach number field
		vtkFile << "SCALARS Mach_Number double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << fixed << setprecision(4) << Mach[i][j] << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n';
				}
			}
		}

		// Temperature field
		vtkFile << "SCALARS Temperature double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << fixed << setprecision(4) << T[i][j] << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n';
				}
			}
		}

		// Velocity magnitude
		vtkFile << "SCALARS Velocity_Magnitude double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					double vel_mag = sqrt(u[i][j] * u[i][j] + v[i][j] * v[i][j]);
					vtkFile << fixed << setprecision(6) << vel_mag << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n';
				}
			}
		}

		// Velocity vector field
		vtkFile << "VECTORS Velocity double" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << scientific << setprecision(6) << u[i][j] << " " << v[i][j] << " 0.0" << '\n';
				}
				else
				{
					vtkFile << "0.0 0.0 0.0" << '\n';
				}
			}
		}

		// Conservative variables (useful for CFD analysis)
		vtkFile << "SCALARS Conservative_U1 double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << scientific << setprecision(6) << U1[i][j] << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n';
				}
			}
		}

		vtkFile << "SCALARS Conservative_U4 double 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
				{
					vtkFile << scientific << setprecision(6) << U4[i][j] << '\n';
				}
				else
				{
					vtkFile << "0.0" << '\n';
				}
			}
		}

		// Immersed boundary and domain indicators
		vtkFile << "SCALARS Domain_Type int 1" << '\n';
		vtkFile << "LOOKUP_TABLE default" << '\n';
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if (Fl[i][j] == 1)
				{
					vtkFile << "1" << '\n'; // Fluid region
				}
				else if (Gst[i][j] == 1)
				{
					vtkFile << "2" << '\n'; // Ghost cells
				}
				else
				{
					vtkFile << "0" << '\n'; // Solid region (inside immersed boundary)
				}
			}
		}

		// MPI rank information (for debugging domain decomposition).
		// Cell (i,j) belongs to rank r whose owned range is [first_r, last_r].
		vtkFile << "SCALARS MPI_Rank int 1\nLOOKUP_TABLE default\n";
		const int base = Nx / mpi_size;
		const int rem  = Nx % mpi_size;
		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				int cell_rank = 0;
				int acc = 2;
				for (int r = 0; r < mpi_size; ++r) {
					int local_n = base + (r < rem ? 1 : 0);
					if (i >= acc && i < acc + local_n) { cell_rank = r; break; }
					acc += local_n;
				}
				vtkFile << cell_rank << '\n';
			}
		}

		vtkFile.close();
		cout << "=== VTK OUTPUT COMPLETE ===" << '\n';
		cout << "File: " << fileName << '\n';
		cout << "Grid: " << grid_nx << " x " << grid_ny << " points" << '\n';
		cout << "Fields: Density, Pressure, Mach, Temperature, Velocity, Conservative Variables" << '\n';
		cout << "Visualization: Open in ParaView, VisIt, or any VTK-compatible viewer" << '\n';
		cout << "=== VISUALIZATION TIPS ===" << '\n';
		cout << "• Use 'Mach_Number' for contour plots" << '\n';
		cout << "• Use 'Velocity' vectors for streamlines" << '\n';
		cout << "• Use 'Domain_Type' to see immersed boundary" << '\n';
		cout << "• Use 'MPI_Rank' to visualize domain decomposition" << '\n';
	}

	// Ensure all processes wait for VTK output to complete
	MPI_Barrier(MPI_COMM_WORLD);
}

/*****************************************************************************************/
/** Quick VTK Export Function - Can be called anytime during simulation */
/*****************************************************************************************/
void exportCurrentStateToVTK(const string &customName)
{
	if (mpi_rank == 0)
	{
		string fileName;
		if (customName.empty())
		{
			fileName = "CFD_Current_State_" + to_string(Nx) + "x" + to_string(Ny) + ".vtk";
		}
		else
		{
			fileName = customName + ".vtk";
		}

		cout << "Exporting current CFD state to: " << fileName << '\n';
	}

	// Use the main VTK function with a special timestep indicator
	outputVTK(-1);
}

/*****************************************************************************************/
/** OUTPUT SOLUTION */
/*****************************************************************************************/
void outputSolution(int timeStep)
{
	gatherDynamicFieldsToRank0();
	if (mpi_rank == 0)
	{
		string fileName = "Aerospike_IBM_MPI_Solution_" + to_string(Nx) + "x" + to_string(Ny) + "_t" + to_string(timeStep) + ".dat";
		ofstream Solution(fileName);
		Solution << "TITLE = Flow Solution\nVARIABLES = X, Y, Rho, Mach, P, u, v\n";
		Solution << "Zone T = Omega I = " << Nx << " J = " << Ny - 2 << '\n';

		for (j = 3; j <= Ny; j++)
		{
			for (i = 2; i <= Nx + 1; i++)
			{
				if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
					Solution << Cx[i][j] << '\t' << Cy[i][j] << '\t' << Rho[i][j] << '\t' << Mach[i][j] << '\t' << P[i][j] << '\t' << u[i][j] << '\t' << v[i][j] << '\n';
				else
					Solution << Cx[i][j] << '\t' << Cy[i][j] << "\t0\t0\t0\t0\t0\n";
			}
		}
		Solution.close();
		cout << "Solution output written to " << fileName << '\n';
	}
}
