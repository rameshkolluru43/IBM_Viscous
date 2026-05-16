#ifndef PHYSICS_H
#define PHYSICS_H

#include <cmath>

/*****************************************************************************************/
/** Common physics helpers shared between the serial and MPI solvers.                     */
/*                                                                                        */
/* Auxiliary-variable update from primitive (rho, u, v, P).                               */
/*   E    = total specific energy                                                         */
/*   a    = speed of sound                                                                */
/*   H    = specific enthalpy                                                             */
/*   T    = temperature                                                                   */
/*   Mach = |v|/a (uses velocity magnitude, not just u, so it is correct in subsonic /    */
/*                 off-axis flow)                                                         */
/*                                                                                        */
/* Constants GAMMA, GAS_CONSTANT, GAMMA_MINUS_1, GAMMA_OVER_GAMMA_MINUS_1 must be defined */
/* before including this header.                                                          */
/*****************************************************************************************/

inline void updateAuxiliary(double rho, double u_, double v_, double P_,
                            double &E_, double &a_, double &H_, double &T_, double &Mach_)
{
	const double q2 = u_ * u_ + v_ * v_;
	E_    = P_ / (rho * GAMMA_MINUS_1) + 0.5 * q2;
	a_    = std::sqrt(GAMMA * P_ / rho);
	H_    = (P_ / rho) * GAMMA_OVER_GAMMA_MINUS_1 + 0.5 * q2;
	T_    = P_ / (rho * GAS_CONSTANT);
	Mach_ = std::sqrt(q2) / a_;
}

#endif // PHYSICS_H
