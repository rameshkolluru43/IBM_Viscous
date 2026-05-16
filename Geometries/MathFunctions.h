#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

#include <algorithm>
#include <cmath>

template <typename T>
inline T Maximum2(T a1, T a2) { return std::max(a1, a2); }

template <typename T>
inline T Minimum2(T a1, T a2) { return std::min(a1, a2); }

template <typename T>
inline T Maximum3(T a1, T a2, T a3) { return std::max({a1, a2, a3}); }

template <typename T>
inline T Minimum3(T a1, T a2, T a3) { return std::min({a1, a2, a3}); }

template <typename T>
inline T Maximum4(T a1, T a2, T a3, T a4) { return std::max({a1, a2, a3, a4}); }

template <typename T>
inline T Minimum4(T a1, T a2, T a3, T a4) { return std::min({a1, a2, a3, a4}); }

template <typename T>
inline T Maximum5(T a1, T a2, T a3, T a4, T a5) { return std::max({a1, a2, a3, a4, a5}); }

template <typename T>
inline T Minimum5(T a1, T a2, T a3, T a4, T a5) { return std::min({a1, a2, a3, a4, a5}); }

template <typename T>
inline T Max6(T a1, T a2, T a3, T a4, T a5, T a6) { return std::max({a1, a2, a3, a4, a5, a6}); }

template <typename T>
inline T Min6(T a1, T a2, T a3, T a4, T a5, T a6) { return std::min({a1, a2, a3, a4, a5, a6}); }

template <typename T>
inline T MaxEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny)
{
	T L1r = std::abs(ur * nx + vr * ny + ar);
	T L1l = std::abs(ul * nx + vl * ny + al);
	T L2r = std::abs(ur * nx + vr * ny);
	T L2l = std::abs(ul * nx + vl * ny);
	T L3r = std::abs(ur * nx + vr * ny - ar);
	T L3l = std::abs(ul * nx + vl * ny - al);
	return Maximum2(Maximum3(L1r, L2r, L3r), Maximum3(L1l, L2l, L3l));
}

template <typename T>
inline T MinEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny)
{
	T L1r = std::abs(ur * nx + vr * ny + ar);
	T L1l = std::abs(ul * nx + vl * ny + al);
	T L2r = std::abs(ur * nx + vr * ny);
	T L2l = std::abs(ul * nx + vl * ny);
	T L3r = std::abs(ur * nx + vr * ny - ar);
	T L3l = std::abs(ul * nx + vl * ny - al);
	return Maximum2(Minimum3(L1r, L2r, L3r), Minimum3(L1l, L2l, L3l));
}

template <typename T>
inline T AvgEigVal(T ur, T ul, T vr, T vl, T ar, T al, T nx, T ny, T ds)
{
	return ((std::abs(ur * nx + vr * ny) + ar) + (std::abs(ul * nx + vl * ny) + al)) * ds * T(0.5);
}

template <typename T>
inline T MoversN(T Gl, T Gr, T Ul, T Ur, T amax, T amin)
{
	const T denom = Ul - Ur;
	const T eps = T(1e-12);
	T S = (std::fabs(denom) > eps) ? std::fabs((Gl - Gr) / denom) : std::fabs(amin);
	if (S > amax) return std::fabs(amax);
	if (S < amin) return std::fabs(amin);
	return std::fabs(S);
}

#endif // MATHFUNCTIONS_H
