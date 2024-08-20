#include "ulp.h"

// The domains of floating point operations are separated into "additive" operations,
// which use register group F and "multiplicative" operations, which use register group E.
//
// [-3.0e+14, 3.0e+14]

static inline double sum_residue(double a, double b, double c) {
	double delta_a = a - (c - b);
	double delta_b = b - (c - a);
	double res = delta_a + delta_b;
	return res;
}

// ties to even
double fadd_0(double a, double b) {
	return a + b;
}

// toward negative infinity
double fadd_1(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);

	if (res < 0.0) {
		return nextafter_1_reg_f(c);
	}

	return c;
}

// toward positive infinity
double fadd_2(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);
	
	if (res > 0.0) {
		return nextafter_2_reg_f(c);
	}

	return c;
}

// toward zero
double fadd_3(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);

	if ((res > 0.0 && c < 0.0) || (res < 0.0 && c > 0.0)) {
		return nextafter_3_reg_f(c);
	}

	return c;
}

// ties to even
double fsub_0(double a, double b) {
	return a - b;
}

// toward negative infinity
double fsub_1(double a, double b) {
	return fadd_1(a, -b);
}

// toward positive infinity
double fsub_2(double a, double b) {
	return fadd_2(a, -b);
}

// toward zero
double fsub_3(double a, double b) {
	return fadd_3(a, -b);
}
