#ifndef IBM_BOUNDARY_H
#define IBM_BOUNDARY_H

/*****************************************************************************************/
/** Immersed-boundary wall models (used when USE_NAVIER_STOKES is enabled).              */
/**                                                                                       */
/** Velocity: stationary no-slip  u_g = -u_im  (after image interpolation).               */
/** Thermal (pick one in SolverOptions.h):                                                */
/**   IBM_WALL_ADIABATIC  — ∂T/∂n = 0  ⇒  T_g = T_im                                      */
/**   IBM_WALL_FIXED_TW   — T = T_w at wall  ⇒  T_g = 2 T_w - T_im                        */
/*****************************************************************************************/

#include "SolverOptions.h"

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
