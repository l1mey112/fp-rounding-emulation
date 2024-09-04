#include "ulp.h"
#include <stdint.h>
#include <stdio.h>

// The domains of floating point operations are separated into "additive" operations,
// which use register group F and "multiplicative" operations, which use register group E.
//
// [-3.0e+14, 3.0e+14]

int64_t fadd_1_visit, fadd_2_visit, fadd_3_visit;
int64_t fadd_1_taken, fadd_2_taken, fadd_3_taken;

__attribute__((destructor))
void destructor_f() {
	printf("fadd_1: %g%% taken (%lu/%lu)\n", (double)fadd_1_taken / (double)fadd_1_visit * 100.0, fadd_1_taken, fadd_1_visit);
	printf("fadd_2: %g%% taken (%lu/%lu)\n", (double)fadd_2_taken / (double)fadd_2_visit * 100.0, fadd_2_taken, fadd_2_visit);
	printf("fadd_3: %g%% taken (%lu/%lu)\n", (double)fadd_3_taken / (double)fadd_3_visit * 100.0, fadd_3_taken, fadd_3_visit);
}

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

	fadd_1_visit++;
	if (res < 0.0) {
		fadd_1_taken++;
		return nextafter_1_reg_f(c);
	}

	return c;
}

// toward positive infinity
double fadd_2(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);
	
	fadd_2_visit++;
	if (res > 0.0) {
		fadd_2_taken++;
		return nextafter_2_reg_f(c);
	}

	return c;
}

// toward zero
double fadd_3(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);

	fadd_3_visit++;
	if ((res > 0.0 && c < 0.0) || (res < 0.0 && c > 0.0)) {
		fadd_3_taken++;
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
