/*****************************************************************************************/
/** 2D EULER CODE- MOVERS-N - OPTIMIZED VERSION */
/** Shrinath K S */
/** Aero- IISC */
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

// Performance pragmas for compiler optimization (GCC only; Clang ignores)
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC optimize("O3,fast-math,unroll-loops")
#pragma GCC target("native")
#endif

// Constants - replace problematic defines
const int MAX_GRID_X = 1500;
const int MAX_GRID_Y = 1500;
const double GAMMA = 1.4;
const double GAMMA_MINUS_1 = GAMMA - 1.0;
const double GAMMA_OVER_GAMMA_MINUS_1 = GAMMA / GAMMA_MINUS_1;
const double INV_GAMMA_MINUS_1 __attribute__((unused)) = 1.0 / GAMMA_MINUS_1;
const double GAS_CONSTANT = 1.0;
const double PI __attribute__((unused)) = 3.141592653589793;
const double CFL_NUMBER = 0.3;
const int MAX_ITERATIONS = 75000;
const int OUTPUT_FREQUENCY = 2000;
#include "MathFunctions.h"
#include "Physics.h"
#include "SolverOptions.h"
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
/** Forward Declarations */
/*****************************************************************************************/
bool isValidSolution();
bool validateFileOpen(const std::ifstream &file, const std::string &filename);
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

	// File names
	std::string gridFile = "GridPoints.dat";
	std::string ibmFile = "ImmersedBoundary.dat";
	std::string outputPrefix = "Aerospike_IBM_Solution";

	// Convergence parameters
	double convergenceTolerance = 1e-6;
	int validationFrequency = 1000;
};

static SolverConfiguration config;
static float Xmin __attribute__((unused)) = 0.0, Xmax __attribute__((unused)) = 3.0, Ymin __attribute__((unused)) = 0, Ymax __attribute__((unused)) = 1.0;
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
static double Geomx[MAX_GRID_X][MAX_GRID_Y], Geomy[MAX_GRID_X][MAX_GRID_Y], pg, ug, rhog, ag, vg, Machg;
static double Fl[MAX_GRID_X][MAX_GRID_Y], Gst[MAX_GRID_X][MAX_GRID_Y], Rad[MAX_GRID_X][MAX_GRID_Y];
static double dx __attribute__((unused)) = 0, dy __attribute__((unused)) = 0, a1 __attribute__((unused)), b1 __attribute__((unused)), dy1 __attribute__((unused)), dy2 __attribute__((unused)), delx __attribute__((unused)), dely __attribute__((unused)), delxy __attribute__((unused));
static double Vtan, Vnor, Vtangh, Vnorgh, Vr __attribute__((unused)), Vrn __attribute__((unused)), C __attribute__((unused)), L = CFL_NUMBER; // L is CFL number
static double Ct, Cr, Cl, Cb, alpha __attribute__((unused)), alpha1, alpha2, alpha3, alpha4, dex __attribute__((unused)), dey __attribute__((unused)), amax, amin;
static double c1 __attribute__((unused)) = 0, Resn1 = 0, Resn2 = 0, Resn3 = 0, Resn4 = 0;
static int i, j, timeStep, Nx = 120, Ny = 42, im, jm, ImBx[MAX_GRID_X][MAX_GRID_Y], ImBy[MAX_GRID_X][MAX_GRID_Y], ImGx[MAX_GRID_X][MAX_GRID_Y], ImGy[MAX_GRID_X][MAX_GRID_Y];
static int gridI __attribute__((unused)), gridJ __attribute__((unused));
static int imageI __attribute__((unused)), imageJ __attribute__((unused));
static double densityInf __attribute__((unused)), pressureInf __attribute__((unused)), velocityXInf __attribute__((unused)), velocityYInf __attribute__((unused));
static double temperatureInf __attribute__((unused)), speedOfSoundInf __attribute__((unused)), enthalpyInf __attribute__((unused)), energyInf __attribute__((unused));

static double k1 __attribute__((unused)), k2 __attribute__((unused)), k3 __attribute__((unused)), k4 __attribute__((unused)), k5 __attribute__((unused)), k6 __attribute__((unused)), N __attribute__((unused)) = Nx * Ny;
static double X1, X2, Y1, Y2, X3, Y3, X4, Y4, um, vm, Px, Py;

/*****************************************************************************************/
/** Initial Boundary Variables - Now using configuration */
/*****************************************************************************************/
static double M_inf = (TEST_CHANNEL_MACH > 0.0) ? TEST_CHANNEL_MACH : config.machInfinity;
static double rho_inf = config.densityInfinity, p_inf = config.pressureInfinity;
static double T_inf = p_inf / (rho_inf * GAS_CONSTANT);
static double a_inf = sqrt(GAMMA * GAS_CONSTANT * T_inf);
static double u_inf = a_inf * M_inf;
static double v_inf = 0.0;
static double H_inf = (p_inf / rho_inf) * GAMMA_OVER_GAMMA_MINUS_1 + 0.5 * (u_inf * u_inf + v_inf * v_inf);
static double E_inf = (p_inf / (rho_inf * GAMMA_MINUS_1)) + 0.5 * (u_inf * u_inf + v_inf * v_inf);

/*****************************************************************************************/
/** Utility Functions - Error Checking and Validation */
/*****************************************************************************************/

/**
 * @brief Check if solution contains NaN or infinite values
 * @return true if solution is valid, false otherwise
 */
bool isValidSolution()
{
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = 2; i <= Nx; i++)
		{
			if (std::isnan(P[i][j]) || std::isinf(P[i][j]) ||
				std::isnan(u[i][j]) || std::isinf(u[i][j]) ||
				std::isnan(v[i][j]) || std::isinf(v[i][j]) ||
				std::isnan(Rho[i][j]) || std::isinf(Rho[i][j]))
			{
				std::cerr << "Invalid solution detected at (" << i << "," << j << ")" << std::endl;
				std::cerr << "P=" << P[i][j] << " u=" << u[i][j] << " v=" << v[i][j] << " rho=" << Rho[i][j] << std::endl;
				return false;
			}
		}
	}
	return true;
}

#if USE_NAVIER_STOKES
static inline void serial_metric_grad_at(int i, int j,
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

static void calculateViscousFluxes()
{
	const double mu = LAMINAR_VISCOSITY;
	const double kappa = thermalConductivity(mu, GAMMA, GAS_CONSTANT, GAMMA_MINUS_1, PRANDTL_NUMBER);
	for (int j = 2; j <= Ny + 1; j++)
	{
		for (int i = 2; i <= Nx + 1; i++)
		{
			VF1R[i][j] = VF2R[i][j] = VF3R[i][j] = VF4R[i][j] = 0.0;
			VF1T[i][j] = VF2T[i][j] = VF3T[i][j] = VF4T[i][j] = 0.0;
			VF1L[i][j] = VF2L[i][j] = VF3L[i][j] = VF4L[i][j] = 0.0;
			VF1B[i][j] = VF2B[i][j] = VF3B[i][j] = VF4B[i][j] = 0.0;
		}
	}
#if defined(_OPENMP)
	#pragma omp parallel for collapse(2) schedule(static)
#endif
	for (int j = 3; j <= Ny; j++)
	{
		for (int i = 3; i <= Nx - 1; i++)
		{
			if (Fl[i][j] != 1)
				continue;
			double ux1, uy1, vx1, vy1, Tx1, Ty1;
			double ux2, uy2, vx2, vy2, Tx2, Ty2;
			double f1, f2, f3, f4;
			if (viscousFaceActive(Fl[i][j], Gst[i][j], Fl[i + 1][j], Gst[i + 1][j]))
			{
				serial_metric_grad_at(i, j, ux1, uy1, vx1, vy1, Tx1, Ty1);
				serial_metric_grad_at(i + 1, j, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2), uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2), vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2), Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i][j] + u[i + 1][j]), vc = 0.5 * (v[i][j] + v[i + 1][j]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxr[i][j], nyr[i][j], f1, f2, f3, f4);
				VF1R[i][j] = f1;
				VF2R[i][j] = f2;
				VF3R[i][j] = f3;
				VF4R[i][j] = f4;
			}
			if (viscousFaceActive(Fl[i][j], Gst[i][j], Fl[i][j + 1], Gst[i][j + 1]))
			{
				serial_metric_grad_at(i, j, ux1, uy1, vx1, vy1, Tx1, Ty1);
				serial_metric_grad_at(i, j + 1, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2), uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2), vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2), Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i][j] + u[i][j + 1]), vc = 0.5 * (v[i][j] + v[i][j + 1]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxt[i][j], nyt[i][j], f1, f2, f3, f4);
				VF1T[i][j] = f1;
				VF2T[i][j] = f2;
				VF3T[i][j] = f3;
				VF4T[i][j] = f4;
			}
			if (viscousFaceActive(Fl[i - 1][j], Gst[i - 1][j], Fl[i][j], Gst[i][j]))
			{
				serial_metric_grad_at(i - 1, j, ux1, uy1, vx1, vy1, Tx1, Ty1);
				serial_metric_grad_at(i, j, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2), uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2), vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2), Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i - 1][j] + u[i][j]), vc = 0.5 * (v[i - 1][j] + v[i][j]);
				viscousFluxDotN(mu, kappa, ux, uy, vx, vy, Tx, Ty, uc, vc, nxl[i][j], nyl[i][j], f1, f2, f3, f4);
				VF1L[i][j] = f1;
				VF2L[i][j] = f2;
				VF3L[i][j] = f3;
				VF4L[i][j] = f4;
			}
			if (viscousFaceActive(Fl[i][j - 1], Gst[i][j - 1], Fl[i][j], Gst[i][j]))
			{
				serial_metric_grad_at(i, j - 1, ux1, uy1, vx1, vy1, Tx1, Ty1);
				serial_metric_grad_at(i, j, ux2, uy2, vx2, vy2, Tx2, Ty2);
				const double ux = 0.5 * (ux1 + ux2), uy = 0.5 * (uy1 + uy2);
				const double vx = 0.5 * (vx1 + vx2), vy = 0.5 * (vy1 + vy2);
				const double Tx = 0.5 * (Tx1 + Tx2), Ty = 0.5 * (Ty1 + Ty2);
				const double uc = 0.5 * (u[i][j - 1] + u[i][j]), vc = 0.5 * (v[i][j - 1] + v[i][j]);
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

/**
 * @brief Validate file opening operation
 * @param file Input file stream
 * @param filename Name of the file being opened
 * @return true if file opened successfully
 */
bool validateFileOpen(const std::ifstream &file, const std::string &filename)
{
	if (!file.is_open())
	{
		std::cerr << "Error: Cannot open file " << filename << std::endl;
		return false;
	}
	return true;
}

/*****************************************************************************************/
/** SOLVER BEGINS */
/*****************************************************************************************/
int main()
{
	// cout<<"M"<<M_inf<<"\t""rho"<<rho_inf<<"\t""P"<<p_inf<<"\t""a"<<a_inf<<"\t""v"<<v_inf<<"\t""u"<<u_inf<<"\t""H"<<H_inf<<"\t""E"<<E_inf<<endl;
	// cout<<"T"<<T_inf<<endl;

	/*****************************************************************************************/
	/** GRID  POINTS & CELL CENTRES */
	/*****************************************************************************************/
	ifstream Grid;
	Grid.open(config.gridFile);
	if (!validateFileOpen(Grid, config.gridFile))
	{
		return -1;
	}
	cout << "Reading Grid from the Grid Data" << endl;
	Grid >> Nx >> Ny;
	cout << "Grid Size is" << Nx << "x" << Ny << endl;
	while (!Grid.eof())
	{
		for (int j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx + 1; i++)
			{
				Grid >> X[i][j] >> Y[i][j];
				// cout<<X[i][j]<<"\t"<<Y[i][j]<<endl;
			}
		}
	}
	ifstream Cells;
	Grid.close();
	Cells.open(config.ibmFile);
	if (!validateFileOpen(Cells, config.ibmFile))
	{
		return -1;
	}
	cout << "Reading Cell Centres and Immersed Boundary Data" << endl;
	Cells >> Nx >> Ny;
	cout << "Grid Size of Cells is" << Nx << "x" << Ny << endl;
	while (!Cells.eof())
	{
		for (int j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx + 1; i++)
			{
				Cells >> Cx[i][j] >> Cy[i][j] >> Fl[i][j] >> Gst[i][j] >> Rad[i][j] >> ImNx[i][j] >> ImNy[i][j] >> ImX[i][j] >> ImY[i][j] >> ImBx[i][j] >> ImBy[i][j] >> Geomx[i][j] >> Geomy[i][j] >> ImGx[i][j] >> ImGy[i][j];
				// cout<<Cx[i][j]<<"\t"<<Cy[i][j]<<"\t"<<Fl[i][j]<<"\t"<<Gst[i][j]<<"\t"<<Rad[i][j]<<endl;
			}
		}
	}
	Cells.close();

	/* ********************************************************************************************** */
	/** AREA CALCULATION */
	/* ********************************************************************************************** */
	for (i = 2; i <= Nx + 1; i++)
	{
		for (j = 2; j <= Ny + 1; j++)
		{
			Area[i][j] = 0.5 * fabs(((X[i][j + 1] - X[i + 1][j]) * (Y[i][j] - Y[i + 1][j + 1])) - ((Y[i][j + 1] - Y[i + 1][j]) * (X[i][j] - X[i + 1][j + 1])));
		}
	}
	ofstream RampCG("RampCellCentres.dat");
	RampCG << "TITLE = Flow\nVARIABLES = X, Y, rho\n";
	RampCG << "Zone T = Omega I = " << Nx + 1 << " J = " << Ny + 1 << '\n';
	for (j = 2; j <= Ny + 1; j++)
	{
		for (i = 2; i <= Nx + 1; i++)
		{
			RampCG << Cx[i][j] << '\t' << Cy[i][j] << "\t0\n";
		}
	}
	/* ********************************************************************************************** */
	/** SideLength and Normal Vector CALCULATION */
	/* ********************************************************************************************** */
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
	/* ********************************************************************************************** */
	/** DOMAIN INITIALISATION */
	/* ********************************************************************************************** */
	for (j = 0; j <= Ny + 2; j++)
	{
		for (i = 0; i <= Nx + 2; i++)
		{
			// Basic initialization to freestream conditions
			Rho[i][j] = rho_inf;
			u[i][j] = u_inf;
			v[i][j] = v_inf;
			P[i][j] = p_inf;
			H[i][j] = H_inf;
			a[i][j] = a_inf;
			E[i][j] = E_inf;
			T[i][j] = T_inf;
			Mach[i][j] = u[i][j] / a[i][j];
			U1[i][j] = rho_inf;
			U2[i][j] = rho_inf * u_inf;
			U3[i][j] = rho_inf * v_inf;
			U4[i][j] = rho_inf * E_inf;

			// Add small perturbation near body to trigger flow development
			// Body is located around x=0, so add perturbation in that region
			double x_pos = X[i][j];
			double y_pos = Y[i][j];
			double dist_from_body = sqrt(x_pos * x_pos + y_pos * y_pos);

			// Add 1% pressure perturbation in region 0.1 < r < 2.0 from body center
			if (dist_from_body > 0.1 && dist_from_body < 2.0)
			{
				double perturbation_factor = 1.0 + 0.01 * exp(-dist_from_body) * cos(3.14159 * x_pos);
				P[i][j] = p_inf * perturbation_factor;
				Rho[i][j] = rho_inf * perturbation_factor; // Isentropic relation

				// Recalculate derived quantities
				a[i][j] = sqrt(GAMMA * P[i][j] / Rho[i][j]);
				E[i][j] = (P[i][j] / (Rho[i][j] * GAMMA_MINUS_1)) + 0.5 * (u[i][j] * u[i][j] + v[i][j] * v[i][j]);
				H[i][j] = (P[i][j] / Rho[i][j]) * GAMMA_OVER_GAMMA_MINUS_1 + 0.5 * (u[i][j] * u[i][j] + v[i][j] * v[i][j]);
				T[i][j] = P[i][j] / (Rho[i][j] * GAS_CONSTANT);
				Mach[i][j] = u[i][j] / a[i][j];

				// Update conservative variables
				U1[i][j] = Rho[i][j];
				U2[i][j] = Rho[i][j] * u[i][j];
				U3[i][j] = Rho[i][j] * v[i][j];
				U4[i][j] = Rho[i][j] * E[i][j];
			}
		}
	}
	cout << "Domain Initialised with flow perturbation" << endl;
	/* ********************************************************************************************** */
	/** SOLVER BEGINS */
	/* ********************************************************************************************** */
	ofstream Residual("Residues.dat");
	for (timeStep = 1; timeStep < MAX_ITERATIONS + 1; timeStep++)
	{
		/* ********************************************************************************************** */
		/** BOUNDARY CONDITIONS */
		/* ********************************************************************************************** */
		/** INLET */
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

		/**BOTTOM SOLID WALL **/
		for (j = 2; j <= Ny; j++)
		{
			for (i = 2; i <= Nx; i++)
			{
				if (Gst[i][j] == 1)
				{
					im = ImBx[i][j] + 1;
					jm = ImBy[i][j] + 1;
					X1 = Cx[im][jm];
					Y1 = Cy[im][jm];
					X2 = Cx[im + 1][jm];
					Y2 = Cy[im + 1][jm];
					X3 = Cx[im][jm + 1];
					Y3 = Cy[im][jm + 1];
					X4 = Cx[im + 1][jm + 1];
					Y4 = Cy[im + 1][jm + 1];
					Px = ImX[i][j];
					Py = ImY[i][j];
					// cout<<X1<<" "<<Y1<<" "<<X2<<" "<<Y2<<" "<<X3<<" "<<Y3<<" "<<X4<<" "<<Y4<<" "<<Px<<" "<<Py<<endl;
					Rho[i][j] = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
					                     Rho[im][jm], Rho[im + 1][jm], Rho[im][jm + 1], Rho[im + 1][jm + 1]);
					um = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
					              u[im][jm], u[im + 1][jm], u[im][jm + 1], u[im + 1][jm + 1]);
					vm = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
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
					ibmSlipVelocity(um, vm, ImNx[i][j], ImNy[i][j], u[i][j], v[i][j]);
					P[i][j] = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py,
					                   P[im][jm], P[im + 1][jm], P[im][jm + 1], P[im + 1][jm + 1]);
					updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
					                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
#endif
				}
			}
		}
		/* Farfield - copy from interior cell along j (not along i on the boundary row,
		 * which created an order-dependent chain). */
		j = 1;
		for (i = 1; i <= Nx + 1; i++)
		{
			Rho[i][j] = Rho[i][j + 1];
			u[i][j]   = u[i][j + 1];
			v[i][j]   = v[i][j + 1];
			P[i][j]   = P[i][j + 1];
			updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
			                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
		}
		j = Ny + 2;
		for (i = 1; i <= Nx + 1; i++)
		{
			Rho[i][j] = Rho[i][j - 1];
			u[i][j]   = u[i][j - 1];
			v[i][j]   = v[i][j - 1];
			P[i][j]   = P[i][j - 1];
			updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
			                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
		}

		/**supersonic Outflow*/
		i = Nx + 1;
		for (j = 1; j <= Ny + 1; j++)
		{

			Rho[i][j] = 2.0 * Rho[i - 1][j] - Rho[i - 2][j];
			u[i][j]   = 2.0 * u[i - 1][j]   - u[i - 2][j];
			v[i][j]   = 2.0 * v[i - 1][j]   - v[i - 2][j];
			P[i][j]   = 2.0 * P[i - 1][j]   - P[i - 2][j];
			updateAuxiliary(Rho[i][j], u[i][j], v[i][j], P[i][j],
			                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
		}
		/* ********************************************************************************************** */
		/** INTERFACE FLUX - OPTIMIZED VERSION */
		/* ********************************************************************************************** */
#if defined(_OPENMP)
		#pragma omp parallel for collapse(2) schedule(static)
#endif
		for (j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx; i++)
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
		/* ********************************************************************************************** */
		/** DISSIPATIVE FLUX */
		/* ********************************************************************************************** */
#if defined(_OPENMP)
		#pragma omp parallel for collapse(2) schedule(static) private(amax,amin,alpha1,alpha2,alpha3,alpha4)
#endif
		for (j = 2; j <= Ny + 1; j++)
		{
			for (i = 2; i <= Nx; i++)
			{
				if (Fl[i][j] == 1)
				{
					amax = MaxEigVal(u[i + 1][j], u[i][j], v[i + 1][j], v[i][j], a[i + 1][j], a[i][j], nxr[i][j], nyr[i][j]);
					amin = MinEigVal(u[i + 1][j], u[i][j], v[i + 1][j], v[i][j], a[i + 1][j], a[i][j], nxr[i][j], nyr[i][j]);

					alpha1 = MoversN(Rho[i + 1][j] * u[i + 1][j] * nxr[i][j] + Rho[i + 1][j] * v[i + 1][j] * nyr[i][j], Rho[i][j] * u[i][j] * nxr[i][j] + Rho[i][j] * v[i][j] * nyr[i][j], Rho[i + 1][j], Rho[i][j], amax, amin);
					alpha2 = MoversN((Rho[i + 1][j] * u[i + 1][j] * u[i + 1][j] + P[i + 1][j]) * nxr[i][j] + (Rho[i + 1][j] * u[i + 1][j] * v[i + 1][j]) * nyr[i][j], (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxr[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyr[i][j], Rho[i + 1][j] * u[i + 1][j], Rho[i][j] * u[i][j], amax, amin);
					alpha3 = MoversN((Rho[i + 1][j] * v[i + 1][j] * v[i + 1][j] + P[i + 1][j]) * nyr[i][j] + (Rho[i + 1][j] * u[i + 1][j] * v[i + 1][j]) * nxr[i][j], (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyr[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxr[i][j], Rho[i + 1][j] * v[i + 1][j], Rho[i][j] * v[i][j], amax, amin);
					alpha4 = MoversN(Rho[i + 1][j] * u[i + 1][j] * H[i + 1][j] * nxr[i][j] + Rho[i + 1][j] * v[i + 1][j] * H[i + 1][j] * nyr[i][j], Rho[i][j] * u[i][j] * H[i][j] * nxr[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyr[i][j], Rho[i + 1][j] * E[i + 1][j], Rho[i][j] * E[i][j], amax, amin);

					DF1R[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i + 1][j]);
					DF2R[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i + 1][j] * u[i + 1][j]);
					DF3R[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i + 1][j] * v[i + 1][j]);
					DF4R[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i + 1][j] * E[i + 1][j]);

					amax = MaxEigVal(u[i][j + 1], u[i][j], v[i][j + 1], v[i][j], a[i][j + 1], a[i][j], nxt[i][j], nyt[i][j]);
					amin = MinEigVal(u[i][j + 1], u[i][j], v[i][j + 1], v[i][j], a[i][j + 1], a[i][j], nxt[i][j], nyt[i][j]);

					alpha1 = MoversN(Rho[i][j + 1] * u[i][j + 1] * nxt[i][j] + Rho[i][j + 1] * v[i][j + 1] * nyt[i][j], Rho[i][j] * u[i][j] * nxt[i][j] + Rho[i][j] * v[i][j] * nyt[i][j], Rho[i][j + 1], Rho[i][j], amax, amin);
					alpha2 = MoversN((Rho[i][j + 1] * u[i][j + 1] * u[i][j + 1] + P[i][j + 1]) * nxt[i][j] + (Rho[i][j + 1] * u[i][j + 1] * v[i][j + 1]) * nyt[i][j], (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxt[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyt[i][j], Rho[i][j + 1] * u[i][j + 1], Rho[i][j] * u[i][j], amax, amin);
					alpha3 = MoversN((Rho[i][j + 1] * v[i][j + 1] * v[i][j + 1] + P[i][j + 1]) * nyt[i][j] + (Rho[i][j + 1] * u[i][j + 1] * v[i][j + 1]) * nxt[i][j], (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyt[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxt[i][j], Rho[i][j + 1] * v[i][j + 1], Rho[i][j] * v[i][j], amax, amin);
					alpha4 = MoversN(Rho[i][j + 1] * u[i][j + 1] * H[i][j + 1] * nxt[i][j] + Rho[i][j + 1] * v[i][j + 1] * H[i][j + 1] * nyt[i][j], Rho[i][j] * u[i][j] * H[i][j] * nxt[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyt[i][j], Rho[i][j + 1] * E[i][j + 1], Rho[i][j] * E[i][j], amax, amin);

					DF1T[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i][j + 1]);
					DF2T[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i][j + 1] * u[i][j + 1]);
					DF3T[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i][j + 1] * v[i][j + 1]);
					DF4T[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i][j + 1] * E[i][j + 1]);

					amax = MaxEigVal(u[i][j], u[i - 1][j], v[i][j], v[i - 1][j], a[i][j], a[i - 1][j], nxl[i][j], nyl[i][j]);
					amin = MinEigVal(u[i][j], u[i - 1][j], v[i][j], v[i - 1][j], a[i][j], a[i - 1][j], nxl[i][j], nyl[i][j]);

					alpha1 = MoversN(Rho[i - 1][j] * u[i - 1][j] * nxl[i][j] + Rho[i - 1][j] * v[i - 1][j] * nyl[i][j], Rho[i][j] * u[i][j] * nxl[i][j] + Rho[i][j] * v[i][j] * nyl[i][j], Rho[i - 1][j], Rho[i][j], amax, amin);
					alpha2 = MoversN((Rho[i - 1][j] * u[i - 1][j] * u[i - 1][j] + P[i - 1][j]) * nxl[i][j] + (Rho[i - 1][j] * u[i - 1][j] * v[i - 1][j]) * nyl[i][j], (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxl[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyl[i][j], Rho[i - 1][j] * u[i - 1][j], Rho[i][j] * u[i][j], amax, amin);
					alpha3 = MoversN((Rho[i - 1][j] * v[i - 1][j] * v[i - 1][j] + P[i - 1][j]) * nyl[i][j] + (Rho[i - 1][j] * u[i - 1][j] * v[i - 1][j]) * nxl[i][j], (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyl[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxl[i][j], Rho[i - 1][j] * v[i - 1][j], Rho[i][j] * v[i][j], amax, amin);
					alpha4 = MoversN(Rho[i - 1][j] * u[i - 1][j] * H[i - 1][j] * nxl[i][j] + Rho[i - 1][j] * v[i - 1][j] * H[i - 1][j] * nyl[i][j], Rho[i][j] * u[i][j] * H[i][j] * nxl[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyl[i][j], Rho[i - 1][j] * E[i - 1][j], Rho[i][j] * E[i][j], amax, amin);

					DF1L[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i - 1][j]);
					DF2L[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i - 1][j] * u[i - 1][j]);
					DF3L[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i - 1][j] * v[i - 1][j]);
					DF4L[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i - 1][j] * E[i - 1][j]);

					amax = MaxEigVal(u[i][j], u[i][j - 1], v[i][j], v[i][j - 1], a[i][j], a[i][j - 1], nxb[i][j], nyb[i][j]);
					amin = MinEigVal(u[i][j], u[i][j - 1], v[i][j], v[i][j - 1], a[i][j], a[i][j - 1], nxb[i][j], nyb[i][j]);

					alpha1 = MoversN(Rho[i][j - 1] * u[i][j - 1] * nxb[i][j] + Rho[i][j - 1] * v[i][j - 1] * nyb[i][j], Rho[i][j] * u[i][j] * nxb[i][j] + Rho[i][j] * v[i][j] * nyb[i][j], Rho[i][j - 1], Rho[i][j], amax, amin);
					alpha2 = MoversN((Rho[i][j - 1] * u[i][j - 1] * u[i][j - 1] + P[i][j - 1]) * nxb[i][j] + (Rho[i][j - 1] * u[i][j - 1] * v[i][j - 1]) * nyb[i][j], (Rho[i][j] * u[i][j] * u[i][j] + P[i][j]) * nxb[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nyb[i][j], Rho[i][j - 1] * u[i][j - 1], Rho[i][j] * u[i][j], amax, amin);
					alpha3 = MoversN((Rho[i][j - 1] * v[i][j - 1] * v[i][j - 1] + P[i][j - 1]) * nyb[i][j] + (Rho[i][j - 1] * u[i][j - 1] * v[i][j - 1]) * nxb[i][j], (Rho[i][j] * v[i][j] * v[i][j] + P[i][j]) * nyb[i][j] + (Rho[i][j] * u[i][j] * v[i][j]) * nxb[i][j], Rho[i][j - 1] * v[i][j - 1], Rho[i][j] * v[i][j], amax, amin);
					alpha4 = MoversN(Rho[i][j - 1] * u[i][j - 1] * H[i][j - 1] * nxb[i][j] + Rho[i][j - 1] * v[i][j - 1] * H[i][j - 1] * nyb[i][j], Rho[i][j] * u[i][j] * H[i][j] * nxb[i][j] + Rho[i][j] * v[i][j] * H[i][j] * nyb[i][j], Rho[i][j - 1] * E[i][j - 1], Rho[i][j] * E[i][j], amax, amin);

					DF1B[i][j] = 0.5 * fabs(alpha1) * (Rho[i][j] - Rho[i][j - 1]);
					DF2B[i][j] = 0.5 * fabs(alpha2) * (Rho[i][j] * u[i][j] - Rho[i][j - 1] * u[i][j - 1]);
					DF3B[i][j] = 0.5 * fabs(alpha3) * (Rho[i][j] * v[i][j] - Rho[i][j - 1] * v[i][j - 1]);
					DF4B[i][j] = 0.5 * fabs(alpha4) * (Rho[i][j] * E[i][j] - Rho[i][j - 1] * E[i][j - 1]);
				}
			}
		}
#if USE_NAVIER_STOKES
		calculateViscousFluxes();
#endif
		/* ********************************************************************************************** */
		/** Local Timestepping */
		/* ********************************************************************************************** */
		for (j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx; i++)
			{
				Cr = AvgEigVal(u[i + 1][j], u[i][j], v[i + 1][j], v[i][j], a[i + 1][j], a[i][j], nxr[i][j], nyr[i][j], dsr[i][j]);
				Ct = AvgEigVal(u[i][j + 1], u[i][j], v[i][j + 1], v[i][j], a[i][j + 1], a[i][j], nxt[i][j], nyt[i][j], dst[i][j]);
				Cl = AvgEigVal(u[i][j], u[i - 1][j], v[i][j], v[i - 1][j], a[i][j], a[i - 1][j], nxl[i][j], nyl[i][j], dsl[i][j]);
				Cb = AvgEigVal(u[i][j], u[i][j - 1], v[i][j], v[i][j - 1], a[i][j], a[i][j - 1], nxb[i][j], nyb[i][j], dsb[i][j]);

				double dt_ij = L * Area[i][j] / Maximum4(Cr, Ct, Cl, Cb) * 0.5;
#if USE_NAVIER_STOKES
				if (Fl[i][j] == 1)
				{
					const double h_min = Minimum4(dsr[i][j], dst[i][j], dsl[i][j], dsb[i][j]);
					const double dt_vis = viscousDtScale(Rho[i][j], LAMINAR_VISCOSITY, h_min, VISCOUS_CFL_FRACTION);
					if (dt_vis < dt_ij)
						dt_ij = dt_vis;
				}
#endif
				dt[i][j] = dt_ij;
				if (dtg >= dt[i][j])
				{
					dtg = dt[i][j];
				}
				// cout<<dt[i][j]<<"\t"<<dtg;
			}
		}
		TG = TG + dtg;
		/* ********************************************************************************************** */
		/** Global Timestepping */
		/* ********************************************************************************************** */
		/*	for(int j=2;j<=Ny+1;j++){
				for (int i=2;i<=Nx;i++){
					Cr=AvgEigVal(u[i+1][j], u[i][j], v[i+1][j], v[i][j], a[i+1][j], a[i][j], nxr[i][j], nyr[i][j], dsr[i][j]);
					Ct=AvgEigVal(u[i][j+1], u[i][j], v[i][j+1], v[i][j], a[i][j+1], a[i][j], nxt[i][j], nyt[i][j], dst[i][j]);
					Cl=AvgEigVal(u[i][j], u[i-1][j], v[i][j], v[i-1][j], a[i][j], a[i-1][j], nxl[i][j], nyl[i][j], dsl[i][j]);
					Cb=AvgEigVal(u[i][j], u[i][j-1], v[i][j], v[i][j-1], a[i][j], a[i][j-1], nxb[i][j], nyb[i][j], dsb[i][j]);

					dt[i][j]=L*Area[i][j]/Maximum4(Cr, Ct, Cl, Cb)*0.5;
					if
					//if (Resn3>=0.5){dt[i][j]=dt[i][j]*0.5;}
					cout<<dt[i][j]<<"\t";
				}
			}*/
		/* ********************************************************************************************** */
		/** UPDATE FLUX */
		/* ********************************************************************************************** */
#if defined(_OPENMP)
		#pragma omp parallel for collapse(2) schedule(static)
#endif
		for (j = 2; j <= Ny + 1; j++)
		{
			for (i = 2; i <= Nx; i++)
			{
				if (Fl[i][j] == 1)
				{
					const double inv_area_dt = dtg / Area[i][j];
					const double dsr_ = dsr[i][j], dst_ = dst[i][j], dsl_ = dsl[i][j], dsb_ = dsb[i][j];
#if USE_NAVIER_STOKES
					Un1[i][j] = U1[i][j] - inv_area_dt * ((IF1R[i][j] + DF1R[i][j] + VF1R[i][j]) * dsr_ + (IF1T[i][j] + DF1T[i][j] + VF1T[i][j]) * dst_ + (IF1L[i][j] + DF1L[i][j] + VF1L[i][j]) * dsl_ + (IF1B[i][j] + DF1B[i][j] + VF1B[i][j]) * dsb_);
					Un2[i][j] = U2[i][j] - inv_area_dt * ((IF2R[i][j] + DF2R[i][j] + VF2R[i][j]) * dsr_ + (IF2T[i][j] + DF2T[i][j] + VF2T[i][j]) * dst_ + (IF2L[i][j] + DF2L[i][j] + VF2L[i][j]) * dsl_ + (IF2B[i][j] + DF2B[i][j] + VF2B[i][j]) * dsb_);
					Un3[i][j] = U3[i][j] - inv_area_dt * ((IF3R[i][j] + DF3R[i][j] + VF3R[i][j]) * dsr_ + (IF3T[i][j] + DF3T[i][j] + VF3T[i][j]) * dst_ + (IF3L[i][j] + DF3L[i][j] + VF3L[i][j]) * dsl_ + (IF3B[i][j] + DF3B[i][j] + VF3B[i][j]) * dsb_);
					Un4[i][j] = U4[i][j] - inv_area_dt * ((IF4R[i][j] + DF4R[i][j] + VF4R[i][j]) * dsr_ + (IF4T[i][j] + DF4T[i][j] + VF4T[i][j]) * dst_ + (IF4L[i][j] + DF4L[i][j] + VF4L[i][j]) * dsl_ + (IF4B[i][j] + DF4B[i][j] + VF4B[i][j]) * dsb_);
#else
					Un1[i][j] = U1[i][j] - inv_area_dt * ((IF1R[i][j] + DF1R[i][j]) * dsr_ + (IF1T[i][j] + DF1T[i][j]) * dst_ + (IF1L[i][j] + DF1L[i][j]) * dsl_ + (IF1B[i][j] + DF1B[i][j]) * dsb_);
					Un2[i][j] = U2[i][j] - inv_area_dt * ((IF2R[i][j] + DF2R[i][j]) * dsr_ + (IF2T[i][j] + DF2T[i][j]) * dst_ + (IF2L[i][j] + DF2L[i][j]) * dsl_ + (IF2B[i][j] + DF2B[i][j]) * dsb_);
					Un3[i][j] = U3[i][j] - inv_area_dt * ((IF3R[i][j] + DF3R[i][j]) * dsr_ + (IF3T[i][j] + DF3T[i][j]) * dst_ + (IF3L[i][j] + DF3L[i][j]) * dsl_ + (IF3B[i][j] + DF3B[i][j]) * dsb_);
					Un4[i][j] = U4[i][j] - inv_area_dt * ((IF4R[i][j] + DF4R[i][j]) * dsr_ + (IF4T[i][j] + DF4T[i][j]) * dst_ + (IF4L[i][j] + DF4L[i][j]) * dsl_ + (IF4B[i][j] + DF4B[i][j]) * dsb_);
#endif
				}
			}
		}
		/* ********************************************************************************************** */
		/** UPDATE PRIMITIVE VARIABLES                                                            */
		/* NaN checks are deliberately *not* in this hot loop — they are run periodically by      */
		/* isValidSolution() (every config.validationFrequency steps).                             */
		/* ********************************************************************************************** */
#if defined(_OPENMP)
		#pragma omp parallel for collapse(2) schedule(static)
#endif
		for (int j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx; i++)
			{
				if (Fl[i][j] == 1)
				{
					const double rho = Un1[i][j];
					const double inv_rho = 1.0 / rho;
					Rho[i][j] = rho;
					u[i][j]   = Un2[i][j] * inv_rho;
					v[i][j]   = Un3[i][j] * inv_rho;
					E[i][j]   = Un4[i][j] * inv_rho;
					P[i][j]   = rho * GAMMA_MINUS_1 * (E[i][j] - 0.5 * (u[i][j] * u[i][j] + v[i][j] * v[i][j]));
					updateAuxiliary(rho, u[i][j], v[i][j], P[i][j],
					                E[i][j], a[i][j], H[i][j], T[i][j], Mach[i][j]);
				}
			}
		}

		// Validate solution periodically
		if (timeStep % config.validationFrequency == 0)
		{
			if (!isValidSolution())
			{
				std::cerr << "Solution validation failed at time step " << timeStep << std::endl;
				exit(1);
			}
		}
		/* ********************************************************************************************** */
		/** Residue Update */
		/* ********************************************************************************************** */
		Resn1 = Resn2 = Resn3 = Resn4 = 0.0;
		int counter = 0;
#if defined(_OPENMP)
		#pragma omp parallel for collapse(2) schedule(static) \
		    reduction(+:Resn1,Resn2,Resn3,Resn4,counter)
#endif
		for (int j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx; i++)
			{
				if (Fl[i][j] == 1)
				{
					++counter;
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
		const double invN = (counter > 0) ? 1.0 / counter : 1.0;
		/* RMS residuals (independent of grid size) */
		const double r1 = std::sqrt(Resn1 * invN);
		const double r2 = std::sqrt(Resn2 * invN);
		const double r3 = std::sqrt(Resn3 * invN);
		const double r4 = std::sqrt(Resn4 * invN);
		Residual << timeStep << '\t' << r1 << '\t' << r2 << '\t' << r3 << '\t' << r4 << '\t' << TG << '\n';
		if (timeStep == 1 || timeStep % 1000 == 0)
		{
			cout << timeStep << '\t' << r1 << '\t' << r2 << '\t' << r3 << '\t' << r4 << '\t' << TG << '\n';
		}
		// Resn1 =sqrt(Resn1/N);
		//  cout<<Resn1<<endl;
		// if (z ==1){  Resn2 = Resn1;  }
		// Resn3 = fabs(log10((Resn2/Resn1)));
		// cout<<z<<"\t"<<Resn3<<endl;
		/* ********************************************************************************************** */
		/** UPDATE CONSERVED VARIABLES */
		/* ********************************************************************************************** */
		for (int j = 2; j <= Ny + 1; j++)
		{
			for (int i = 2; i <= Nx; i++)
			{
				if (Fl[i][j] == 1)
				{
					U1[i][j] = Un1[i][j];
					U2[i][j] = Un2[i][j];
					U3[i][j] = Un3[i][j];
					U4[i][j] = Un4[i][j];
				}
			}
		}

		/* ********************************************************************************************** */
		/** INTERMITTENT OUTPUT */
		/* ********************************************************************************************** */
		if (timeStep % OUTPUT_FREQUENCY == 0)
		{
			string Fname(config.outputPrefix + "_" + std::to_string(Nx) + "x" + std::to_string(Ny) + "_t" + std::to_string(timeStep) + ".dat");
			ofstream Solution(Fname);
			Solution << "TITLE = Flow Solution\nVARIABLES = X, Y, Rho, Mach, P, u, v\n";
			Solution << "Zone T = Omega I = " << Nx << " J = " << Ny - 2 << '\n';

			for (j = 3; j <= Ny; j++)
			{
				for (i = 2; i <= Nx + 1; i++)
				{
					if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
					{
						Solution << Cx[i][j] << '\t' << Cy[i][j] << '\t' << Rho[i][j] << '\t' << Mach[i][j] << '\t' << P[i][j] << '\t' << u[i][j] << '\t' << v[i][j] << '\n';
					}
					else
					{
						Solution << Cx[i][j] << '\t' << Cy[i][j] << "\t0\t0\t0\t0\t0\n";
					}
				}
			}
		}
		/***********************************************************************************************/
	}
	/*****************************************************************************************/
	/**SOLVER ENDS*/
	/* ********************************************************************************************** */
	/** OUTPUT */
	/* ********************************************************************************************** */
	ofstream Solution("Aerospike_IBMSolution_500x500.dat");
	Solution << "TITLE = Flow Solution\nVARIABLES = X, Y, Rho, Mach, P, u, v\n";
	Solution << "Zone T = Omega I = " << Nx << " J = " << Ny - 2 << '\n';

	for (j = 3; j <= Ny; j++)
	{
		for (i = 2; i <= Nx + 1; i++)
		{
			if ((Fl[i][j] == 1) || (Gst[i][j] == 1))
			{
				Solution << Cx[i][j] << '\t' << Cy[i][j] << '\t' << Rho[i][j] << '\t' << Mach[i][j] << '\t' << P[i][j] << '\t' << u[i][j] << '\t' << v[i][j] << '\n';
			}
			else
			{
				Solution << Cx[i][j] << '\t' << Cy[i][j] << "\t0\t0\t0\t0\t0\n";
			}
		}
	}
	ofstream Boundary("BoundarySolution.dat");
	Boundary << "TITLE = Flow Solution\nVARIABLES = X, Y, Rho, Mach, P, u, v\n";
	for (j = 2; j <= Ny; j++)
	{
		for (i = 2; i <= Nx; i++)
		{
			if (Gst[i][j] == 1)
			{
				im = ImGx[i][j] + 1;
				jm = ImGy[i][j] + 1;
				X1 = Cx[im][jm];
				Y1 = Cy[im][jm];
				X2 = Cx[im + 1][jm];
				Y2 = Cy[im + 1][jm];
				X3 = Cx[im][jm + 1];
				Y3 = Cy[im][jm + 1];
				X4 = Cx[im + 1][jm + 1];
				Y4 = Cy[im + 1][jm + 1];
				Px = Geomx[i][j];
				Py = Geomy[i][j];
				// cout<<X1<<" "<<Y1<<" "<<X2<<" "<<Y2<<" "<<X3<<" "<<Y3<<" "<<X4<<" "<<Y4<<" "<<Px<<" "<<Py<<endl;
				pg = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, Px, Py, P[im][jm], P[im + 1][jm], P[im][jm + 1], P[im + 1][jm + 1]);
				rhog = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, ImX[i][j], ImY[i][j], Rho[im][jm], Rho[im + 1][jm], Rho[im][jm + 1], Rho[im + 1][jm + 1]);
				ug = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, ImX[i][j], ImY[i][j], u[im][jm], u[im + 1][jm], u[im][jm + 1], u[im + 1][jm + 1]);
				vg = Bilinear(X1, Y1, X2, Y2, X3, Y3, X4, Y4, ImX[i][j], ImY[i][j], v[im][jm], v[im + 1][jm], v[im][jm + 1], v[im + 1][jm + 1]);
				ag = sqrt((GAMMA * pg) / rhog);
				Machg = std::sqrt(ug * ug + vg * vg) / ag;
				Boundary << Geomx[i][j] << '\t' << Geomy[i][j] << '\t' << rhog << '\t' << Machg << '\t' << pg << '\t' << ug << '\t' << vg << '\n';
			}
		}
	}
	Boundary.close();
	Solution.close();
	/***********************************************************************************************/
	return 0;
	/***********************************************************************************************/
}
