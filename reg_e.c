#include "ulp.h"

#include <math.h>

// The domains of floating point operations are separated into "additive" operations,
// which use register group F and "multiplicative" operations, which use register group E.
//
// [1.7e-77, INFINITY)

static inline double upper_half(double x) {
	double secator = 134217729.0;
	double p = x * secator;
	return p + (x - p);
}

static inline double emulated_fma(double a, double b, double c) {
	double aup = upper_half(a);
	double alo = a - aup;
	double bup = upper_half(b);
	double blo = b - bup;

	double high = aup * bup; 
	double mid = aup * blo + alo * bup;
	double low = alo * blo;
	double ab = high + mid;
	double resab = (high - ab) + mid + low;

	double fma = ab + c;
	return resab + fma;
}

static inline double mul_residue(double a, double b, double c) {
	return emulated_fma(a, b, -c);
}

static inline double mul_residue_fma(double a, double b, double c) {
	return fma(a, b, -c);
}

// ties to even
double fmul_0(double a, double b) {
	return a * b;
}

// toward negative infinity
double fmul_1(double a, double b) {
	double c = a * b;
	if (!isfinite(c) && isfinite(a) && isfinite(b) && b != 0.0) {
		return nextafter_1(c);
	}
	int expa, expb;
	double a_scaled = frexp(a, &expa);
	double b_scaled = frexp(b, &expb);
	double c_scaled = ldexp(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);

	if (res < 0.0) {
		return nextafter_1(c);
	}

	return c;
}

double fmul_2(double a, double b) {
	double c = a * b;
	if (!isfinite(c) && isfinite(a) && isfinite(b) && b != 0.0) {
		return nextafter_2(c);
	}
	int expa, expb;
	double a_scaled = frexp(a, &expa);
	double b_scaled = frexp(b, &expb);
	double c_scaled = ldexp(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);

	if (res > 0.0) {
		return nextafter_2(c);
	}

	return c;
}

double fmul_3(double a, double b) {
	double c = a * b;
	if (!isfinite(c) && isfinite(a) && isfinite(b) && b != 0.0) {
		return nextafter_3(c);
	}
	int expa, expb;
	double a_scaled = frexp(a, &expa);
	double b_scaled = frexp(b, &expb);
	double c_scaled = ldexp(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);

	if (res < 0.0) {
		return nextafter_3(c);
	}

	return c;
}

// toward negative infinity
double fmul_fma_1(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	if (res < 0.0) {
		return nextafter_1(c);
	}

	return c;
}

double fmul_fma_2(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	if (res > 0.0) {
		return nextafter_2(c);
	}

	return c;
}

double fmul_fma_3(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	if (res < 0.0) {
		return nextafter_3(c);
	}

	return c;
}
