#ifndef VISCOUSFLUX_H
#define VISCOUSFLUX_H

#include <cmath>

/*****************************************************************************************/
/** Newtonian laminar viscous flux for 2D compressible Navier–Stokes (primitive vars).   */
/**                                                                                       */
/** Face flux uses stress tensor τ and Fourier heat flux q = -k∇T dotted with face      */
/** outward normal **n**, matched to the same semi-discrete form as the Euler fluxes:     */
/**   F2_visc = -(τ·n)_x ,  F3_visc = -(τ·n)_y ,                                          */
/**   F4_visc = -(u (τ·n)_x + v (τ·n)_y) + k (n·∇T)   with q = -k∇T ⇒ energy flux term.   */
/** Mass viscous flux = 0.                                                                */
/**                                                                                       */
/** IBM fluid–wall faces use viscousFaceActive() from IbmBoundary.h (include after NS on). */
/*****************************************************************************************/

inline double thermalConductivity(double mu, double gamma, double gas_R, double gmo, double Pr)
{
	const double cp = gamma * gas_R / gmo;
	return mu * cp / Pr;
}

inline void newtonianTau(double mu, double ux, double uy, double vx, double vy,
                         double &tau_xx, double &tau_xy, double &tau_yy)
{
	const double div = ux + vy;
	const double t23 = (2.0 / 3.0) * mu * div;
	tau_xx = 2.0 * mu * ux - t23;
	tau_yy = 2.0 * mu * vy - t23;
	tau_xy = mu * (uy + vx);
}

/** Viscous flux (F1..F4) dotted with face normal (nx,ny); mass flux F1=0. */
inline void viscousFluxDotN(double mu, double kappa,
                            double ux, double uy, double vx, double vy,
                            double Tx, double Ty,
                            double u, double v,
                            double nx, double ny,
                            double &F1, double &F2, double &F3, double &F4)
{
	double tau_xx, tau_xy, tau_yy;
	newtonianTau(mu, ux, uy, vx, vy, tau_xx, tau_xy, tau_yy);
	const double fnx = tau_xx * nx + tau_xy * ny;
	const double fny = tau_xy * nx + tau_yy * ny;
	const double qn = kappa * (Tx * nx + Ty * ny);
	F1 = 0.0;
	F2 = -fnx;
	F3 = -fny;
	F4 = -(u * fnx + v * fny) + qn;
}

/** Explicit viscous stability scale: dt_vis ~ coeff * rho * h^2 / mu with h = edge length. */
inline double viscousDtScale(double rho, double mu, double h_edge, double coeff)
{
	if (mu <= 0.0 || rho <= 0.0 || h_edge <= 0.0)
		return 1.0e300;
	const double h2 = h_edge * h_edge;
	return coeff * rho * h2 / mu;
}

#endif // VISCOUSFLUX_H
