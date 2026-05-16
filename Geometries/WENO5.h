#ifndef WENO5_H
#define WENO5_H

#include <cmath>

/*****************************************************************************************/
/** Fifth-order WENO (Jiang–Shu) scalar reconstruction at a 1D cell interface.            */
/**                                                                                       */
/** For interface between k and k+1:                                                     */
/**   q_minus  uses stencil {k-2, k-1, k, k+1, k+2}  (left-biased reconstruction)        */
/**   q_plus   uses stencil {k-1, k, k+1, k+2, k+3}  (right-biased reconstruction)       */
/**                                                                                       */
/** WENO_EPSILON should be set before include (default 1e-6).                              */
/*****************************************************************************************/

#ifndef WENO_EPSILON
#define WENO_EPSILON 1.0e-6
#endif

inline double weno5_q_minus(double v0, double v1, double v2, double v3, double v4)
{
	const double beta0 = 13.0 / 12.0 * (v0 - 2.0 * v1 + v2) * (v0 - 2.0 * v1 + v2)
	                     + 0.25 * (v0 - 4.0 * v1 + 3.0 * v2) * (v0 - 4.0 * v1 + 3.0 * v2);
	const double beta1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3)
	                     + 0.25 * (v1 - v3) * (v1 - v3);
	const double beta2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4)
	                     + 0.25 * (3.0 * v2 - 4.0 * v3 + v4) * (3.0 * v2 - 4.0 * v3 + v4);

	const double d0 = 0.1;
	const double d1 = 0.6;
	const double d2 = 0.3;

	const double eps = WENO_EPSILON;
	const double a0 = d0 / ((eps + beta0) * (eps + beta0));
	const double a1 = d1 / ((eps + beta1) * (eps + beta1));
	const double a2 = d2 / ((eps + beta2) * (eps + beta2));
	const double wsum = a0 + a1 + a2;
	const double w0 = a0 / wsum;
	const double w1 = a1 / wsum;
	const double w2 = a2 / wsum;

	const double p0 = (1.0 / 3.0) * v0 - (7.0 / 6.0) * v1 + (11.0 / 6.0) * v2;
	const double p1 = (-1.0 / 6.0) * v1 + (5.0 / 6.0) * v2 + (1.0 / 3.0) * v3;
	const double p2 = (1.0 / 3.0) * v2 + (5.0 / 6.0) * v3 - (1.0 / 6.0) * v4;

	return w0 * p0 + w1 * p1 + w2 * p2;
}

inline double weno5_q_plus(double v0, double v1, double v2, double v3, double v4)
{
	const double beta0 = 13.0 / 12.0 * (v0 - 2.0 * v1 + v2) * (v0 - 2.0 * v1 + v2)
	                     + 0.25 * (v0 - 4.0 * v1 + 3.0 * v2) * (v0 - 4.0 * v1 + 3.0 * v2);
	const double beta1 = 13.0 / 12.0 * (v1 - 2.0 * v2 + v3) * (v1 - 2.0 * v2 + v3)
	                     + 0.25 * (v1 - v3) * (v1 - v3);
	const double beta2 = 13.0 / 12.0 * (v2 - 2.0 * v3 + v4) * (v2 - 2.0 * v3 + v4)
	                     + 0.25 * (3.0 * v2 - 4.0 * v3 + v4) * (3.0 * v2 - 4.0 * v3 + v4);

	const double d0 = 0.3;
	const double d1 = 0.6;
	const double d2 = 0.1;

	const double eps = WENO_EPSILON;
	const double a0 = d0 / ((eps + beta0) * (eps + beta0));
	const double a1 = d1 / ((eps + beta1) * (eps + beta1));
	const double a2 = d2 / ((eps + beta2) * (eps + beta2));
	const double wsum = a0 + a1 + a2;
	const double w0 = a0 / wsum;
	const double w1 = a1 / wsum;
	const double w2 = a2 / wsum;

	const double p0 = (-1.0 / 6.0) * v0 + (5.0 / 6.0) * v1 + (1.0 / 3.0) * v2;
	const double p1 = (1.0 / 3.0) * v1 + (5.0 / 6.0) * v2 - (1.0 / 6.0) * v3;
	const double p2 = (-1.0 / 6.0) * v2 + (5.0 / 6.0) * v3 + (1.0 / 3.0) * v4;

	return w0 * p0 + w1 * p1 + w2 * p2;
}

#endif // WENO5_H
