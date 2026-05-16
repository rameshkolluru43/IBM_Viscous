#ifndef HIGHORDERFLUX_H
#define HIGHORDERFLUX_H

/*****************************************************************************************/
/** Inviscid Euler flux F·n̂ for 2D conservative variables (primitive form).               */
/** Used by central and WENO-based interface closures.                                   */
/*****************************************************************************************/

inline void eulerFluxDotN(double rho, double u, double v, double P, double H,
                          double nx, double ny,
                          double &f1, double &f2, double &f3, double &f4)
{
	const double un = u * nx + v * ny;
	f1 = rho * un;
	f2 = rho * u * un + P * nx;
	f3 = rho * v * un + P * ny;
	f4 = rho * H * un;
}

#endif // HIGHORDERFLUX_H
