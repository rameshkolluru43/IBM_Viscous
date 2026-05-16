#ifndef IBM_BOUNDARY_H
#define IBM_BOUNDARY_H

/*****************************************************************************************/
/** Immersed-boundary wall models (used when USE_NAVIER_STOKES is enabled).              */
/**                                                                                       */
/** Velocity: no-slip u_g = -u_im, or wall-function shear matching (USE_WALL_FUNCTION).    */
/** Thermal (pick one in SolverOptions.h):                                                */
/**   IBM_WALL_ADIABATIC  — ∂T/∂n = 0  ⇒  T_g = T_im                                      */
/**   IBM_WALL_FIXED_TW   — T = T_w at wall  ⇒  T_g = 2 T_w - T_im                        */
/*****************************************************************************************/

#include "SolverOptions.h"

#include <algorithm>
#include <cmath>

#if USE_NAVIER_STOKES

/** Face is active for viscous flux: fluid–fluid or fluid–IBM. */
inline bool viscousFaceActive(int Fl0, int Gst0, int Fl1, int Gst1)
{
	const bool fluid0 = (Fl0 == 1);
	const bool fluid1 = (Fl1 == 1);
	const bool ibm0 = (Gst0 == 1);
	const bool ibm1 = (Gst1 == 1);
	return (fluid0 && fluid1) || (fluid0 && ibm1) || (ibm0 && fluid1);
}

/** Stationary no-slip: negate image velocity (antisymmetric about the wall). */
inline void ibmNoSlipVelocity(double um, double vm, double &u_out, double &v_out)
{
	u_out = -um;
	v_out = -vm;
}

/** Wall-normal distance from cell centre (cx,cy) to wall foot (wx, wy). */
inline double ibmWallDistance(double cx, double cy, double wx, double wy)
{
	const double dx = cx - wx;
	const double dy = cy - wy;
	return std::sqrt(dx * dx + dy * dy);
}

/** Solve u_τ from u_∥ and d_w using blended laminar / log law (explicit Newton). */
inline double wallFunctionTauW(double rho, double u_parallel, double d_w, double mu)
{
	if (rho <= 0.0 || mu <= 0.0 || d_w <= 0.0 || u_parallel < 1.0e-14)
		return 0.0;

	const double nu = mu / rho;
	const double kappa = WALL_KAPPA;
	const double E = WALL_E;
	double u_tau = std::sqrt(std::max(mu * u_parallel / (rho * d_w), 1.0e-12));

	const double u_tau_lam = std::sqrt(std::max(mu * u_parallel / (rho * d_w), 0.0));
	const double u_tau_cap = std::max(0.15 * u_parallel, u_tau_lam);

	for (int it = 0; it < WALL_FUNCTION_NEWTON_ITERS; ++it)
	{
		const double yplus = std::max(u_tau * d_w / nu, 1.0e-8);
		double u_plus;
		double du_plus_du_tau;
		if (yplus < 11.0)
		{
			u_plus = yplus;
			du_plus_du_tau = d_w / nu;
		}
		else
		{
			u_plus = (1.0 / kappa) * std::log(yplus) + E;
			du_plus_du_tau = (1.0 / kappa) * (d_w / nu) / yplus;
		}
		const double f = u_parallel - u_tau * u_plus;
		const double df = -(u_plus + u_tau * du_plus_du_tau);
		if (std::fabs(df) < 1.0e-30)
			break;
		u_tau -= f / df;
		u_tau = std::min(std::max(u_tau, 1.0e-10), u_tau_cap);
		if (std::fabs(f) < 1.0e-8 * std::max(u_parallel, 1.0))
			break;
	}
	return rho * u_tau * u_tau;
}

/** Wall-function ghost velocity: u_n mirrored; u_t set from τ_w and (d_w + d_g). */
inline void ibmWallFunctionVelocity(double um, double vm, double nx, double ny,
                                    double rho, double mu, double d_w, double d_g,
                                    double &u_out, double &v_out)
{
	const double un = um * nx + vm * ny;
	double utx = um - un * nx;
	double uty = vm - un * ny;
	const double u_par = std::sqrt(utx * utx + uty * uty);
	const double tau_w = wallFunctionTauW(rho, u_par, d_w, mu);
	const double d_eff = std::max(d_w + d_g, 1.0e-12);
	const double tau_lam = mu * u_par / d_eff;
	const double tau_use = std::min(tau_w, 2.0 * tau_lam);
	double u_t_wall = u_par - tau_use * d_eff / mu;
	u_t_wall = std::max(0.0, std::min(u_t_wall, u_par));

	if (u_par > 1.0e-14)
	{
		const double scale = -u_t_wall / u_par;
		utx *= scale;
		uty *= scale;
	}
	else
	{
		utx = 0.0;
		uty = 0.0;
	}

	const double un_ghost = -un;
	u_out = un_ghost * nx + utx;
	v_out = un_ghost * ny + uty;

	if (!std::isfinite(u_out) || !std::isfinite(v_out) ||
	    (u_out * u_out + v_out * v_out) > 25.0 * (um * um + vm * vm + 1.0e-12))
	{
		ibmNoSlipVelocity(um, vm, u_out, v_out);
	}
}

/** Apply IBM velocity BC (no-slip or wall function). rad_dist > 0 uses locator wall distance. */
inline void ibmWallVelocity(double um, double vm, double nx, double ny,
                            double rho, double mu,
                            double cx_g, double cy_g, double wx, double wy,
                            double rad_dist,
                            double &u_out, double &v_out)
{
#if USE_WALL_FUNCTION
	double d_w = (rad_dist > 1.0e-12) ? rad_dist : ibmWallDistance(cx_g, cy_g, wx, wy);
	d_w = std::max(d_w, 1.0e-12);
	const double d_g = std::max(0.5 * d_w, 1.0e-12);
	ibmWallFunctionVelocity(um, vm, nx, ny, rho, mu, d_w, d_g, u_out, v_out);
#else
	(void)rad_dist;
	(void)nx;
	(void)ny;
	(void)rho;
	(void)mu;
	(void)cx_g;
	(void)cy_g;
	(void)wx;
	(void)wy;
	ibmNoSlipVelocity(um, vm, u_out, v_out);
#endif
}

/** Thermal ghost value from image temperature T_im. */
inline double ibmWallTemperature(double T_im)
{
#if IBM_WALL_FIXED_TW
	return 2.0 * WALL_TEMPERATURE_TW - T_im;
#else
	return T_im;
#endif
}

/** After rho, u, v, T are set at the IBM cell, close the wall state with ideal gas. */
inline void ibmCloseWallThermodynamics(double rho, double u_, double v_, double T_,
                                       double gas_R, double gamma, double gmo,
                                       double &P_, double &E_, double &a_, double &H_,
                                       double &Mach_)
{
	P_ = rho * gas_R * T_;
	const double q2 = u_ * u_ + v_ * v_;
	E_ = P_ / (rho * gmo) + 0.5 * q2;
	a_ = std::sqrt(gamma * P_ / rho);
	H_ = (P_ / rho) * (gamma / gmo) + 0.5 * q2;
	Mach_ = std::sqrt(q2) / a_;
}

#else /* Euler: inviscid slip wall in normal/tangential form */

inline bool viscousFaceActive(int /*Fl0*/, int /*Gst0*/, int /*Fl1*/, int /*Gst1*/)
{
	return false;
}

inline void ibmSlipVelocity(double um, double vm, double nx, double ny,
                            double &u_out, double &v_out)
{
	const double Vtan = ny * um - nx * vm;
	const double Vnor = nx * um + ny * vm;
	const double Vtangh = Vtan;
	const double Vnorgh = -Vnor;
	u_out = Vtangh * ny + Vnorgh * nx;
	v_out = -Vtangh * nx + Vnorgh * ny;
}

#endif

#endif // IBM_BOUNDARY_H
