#include "ulp.h"

#include <stdio.h>
#include <stdint.h>
#include <math.h>

// The domains of floating point operations are separated into "additive" operations,
// which use register group F and "multiplicative" operations, which use register group E.
//
// [1.7e-77, INFINITY)

int64_t fmul_1_visit, fmul_2_visit, fmul_3_visit, fmul_1_taken, fmul_2_taken, fmul_3_taken;
int64_t fmul_fma_1_visit, fmul_fma_2_visit, fmul_fma_3_visit, fmul_fma_1_taken, fmul_fma_2_taken, fmul_fma_3_taken;
int64_t fdiv_1_visit, fdiv_2_visit, fdiv_3_visit, fdiv_1_taken, fdiv_2_taken, fdiv_3_taken;
int64_t fdiv_fma_1_visit, fdiv_fma_2_visit, fdiv_fma_3_visit, fdiv_fma_1_taken, fdiv_fma_2_taken, fdiv_fma_3_taken;
int64_t fsqrt_1_visit, fsqrt_2_visit, fsqrt_3_visit, fsqrt_1_taken, fsqrt_2_taken, fsqrt_3_taken;
int64_t fsqrt_fma_1_visit, fsqrt_fma_2_visit, fsqrt_fma_3_visit, fsqrt_fma_1_taken, fsqrt_fma_2_taken, fsqrt_fma_3_taken;

int64_t fmul_1_isinf, fmul_2_isinf, fmul_3_isinf;

__attribute__((destructor))
void destructor_e() {
	printf("fmul_1: %g%% taken (%lu/%lu)\n", (double)fmul_1_taken / (double)fmul_1_visit * 100.0, fmul_1_taken, fmul_1_visit);
	printf("fmul_2: %g%% taken (%lu/%lu)\n", (double)fmul_2_taken / (double)fmul_2_visit * 100.0, fmul_2_taken, fmul_2_visit);
	printf("fmul_3: %g%% taken (%lu/%lu)\n", (double)fmul_3_taken / (double)fmul_3_visit * 100.0, fmul_3_taken, fmul_3_visit);
	printf("fmul_fma_1: %g%% taken (%lu/%lu)\n", (double)fmul_fma_1_taken / (double)fmul_fma_1_visit * 100.0, fmul_fma_1_taken, fmul_fma_1_visit);
	printf("fmul_fma_2: %g%% taken (%lu/%lu)\n", (double)fmul_fma_2_taken / (double)fmul_fma_2_visit * 100.0, fmul_fma_2_taken, fmul_fma_2_visit);
	printf("fmul_fma_3: %g%% taken (%lu/%lu)\n", (double)fmul_fma_3_taken / (double)fmul_fma_3_visit * 100.0, fmul_fma_3_taken, fmul_fma_3_visit);
	printf("fdiv_1: %g%% taken (%lu/%lu)\n", (double)fdiv_1_taken / (double)fdiv_1_visit * 100.0, fdiv_1_taken, fdiv_1_visit);
	printf("fdiv_2: %g%% taken (%lu/%lu)\n", (double)fdiv_2_taken / (double)fdiv_2_visit * 100.0, fdiv_2_taken, fdiv_2_visit);
	printf("fdiv_3: %g%% taken (%lu/%lu)\n", (double)fdiv_3_taken / (double)fdiv_3_visit * 100.0, fdiv_3_taken, fdiv_3_visit);
	printf("fdiv_fma_1: %g%% taken (%lu/%lu)\n", (double)fdiv_fma_1_taken / (double)fdiv_fma_1_visit * 100.0, fdiv_fma_1_taken, fdiv_fma_1_visit);
	printf("fdiv_fma_2: %g%% taken (%lu/%lu)\n", (double)fdiv_fma_2_taken / (double)fdiv_fma_2_visit * 100.0, fdiv_fma_2_taken, fdiv_fma_2_visit);
	printf("fdiv_fma_3: %g%% taken (%lu/%lu)\n", (double)fdiv_fma_3_taken / (double)fdiv_fma_3_visit * 100.0, fdiv_fma_3_taken, fdiv_fma_3_visit);
	printf("fsqrt_1: %g%% taken (%lu/%lu)\n", (double)fsqrt_1_taken / (double)fsqrt_1_visit * 100.0, fsqrt_1_taken, fsqrt_1_visit);
	printf("fsqrt_2: %g%% taken (%lu/%lu)\n", (double)fsqrt_2_taken / (double)fsqrt_2_visit * 100.0, fsqrt_2_taken, fsqrt_2_visit);
	printf("fsqrt_3: %g%% taken (%lu/%lu)\n", (double)fsqrt_3_taken / (double)fsqrt_3_visit * 100.0, fsqrt_3_taken, fsqrt_3_visit);
	printf("fsqrt_fma_1: %g%% taken (%lu/%lu)\n", (double)fsqrt_fma_1_taken / (double)fsqrt_fma_1_visit * 100.0, fsqrt_fma_1_taken, fsqrt_fma_1_visit);
	printf("fsqrt_fma_2: %g%% taken (%lu/%lu)\n", (double)fsqrt_fma_2_taken / (double)fsqrt_fma_2_visit * 100.0, fsqrt_fma_2_taken, fsqrt_fma_2_visit);
	printf("fsqrt_fma_3: %g%% taken (%lu/%lu)\n", (double)fsqrt_fma_3_taken / (double)fsqrt_fma_3_visit * 100.0, fsqrt_fma_3_taken, fsqrt_fma_3_visit);

	printf("fmul_1: %g%% isinf (%lu/%lu)\n", (double)fmul_1_isinf / (double)fmul_1_visit * 100.0, fmul_1_isinf, fmul_1_visit);
	printf("fmul_2: %g%% isinf (%lu/%lu)\n", (double)fmul_2_isinf / (double)fmul_2_visit * 100.0, fmul_2_isinf, fmul_2_visit);
	printf("fmul_3: %g%% isinf (%lu/%lu)\n", (double)fmul_3_isinf / (double)fmul_3_visit * 100.0, fmul_3_isinf, fmul_3_visit);
}

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
	fmul_1_visit++;
	double c = a * b;
	if (!isfinite(c)) {
		fmul_1_isinf++;
		// fin * fin = _inf_; round down to the nearest representable number
		if (isfinite(a) && isfinite(b)) {
			return nextafter_1_reg_e(c);
		}
		// inf * inf = inf
		return c;
	}
	int expa, expb;
	double a_scaled = frexp_reg_e_nozero_noinf(a, &expa);
	double b_scaled = frexp_reg_e_nozero_noinf(b, &expb);
	double c_scaled = ldexp_reg_e_nozero_noinf(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);

	if (res < 0.0) {
		fmul_1_taken++;
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
double fmul_2(double a, double b) {
	fmul_2_visit++;
	double c = a * b;
	if (!isfinite(c)) {
		fmul_2_isinf++;
		// fin * fin = inf
		// inf * inf = inf
		return c;
	}
	int expa, expb;
	double a_scaled = frexp_reg_e_nozero_noinf(a, &expa);
	double b_scaled = frexp_reg_e_nozero_noinf(b, &expb);
	double c_scaled = ldexp_reg_e_nozero_noinf(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);

	if (res > 0.0) {
		fmul_2_taken++;
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
double fmul_3(double a, double b) {
	fmul_3_visit++;
	double c = a * b;
	if (!isfinite(c)) {
		fmul_3_isinf++;
		// fin * fin = _inf_; round down to the nearest representable number
		if (isfinite(a) && isfinite(b)) {
			return nextafter_3_reg_e(c);
		}
		// inf * inf = inf
		return c;
	}
	int expa, expb;
	double a_scaled = frexp_reg_e_nozero_noinf(a, &expa);
	double b_scaled = frexp_reg_e_nozero_noinf(b, &expb);
	double c_scaled = ldexp_reg_e_nozero_noinf(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);

	if (res < 0.0) {
		fmul_3_taken++;
		return nextafter_3_reg_e(c);
	}

	return c;
}

// toward negative infinity
double fmul_fma_1(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	fmul_fma_1_visit++;
	if (res < 0.0) {
		fmul_fma_1_taken++;
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
double fmul_fma_2(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	fmul_fma_2_visit++;
	if (res > 0.0) {
		fmul_fma_2_taken++;
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
double fmul_fma_3(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	fmul_fma_3_visit++;
	if (res < 0.0) {
		fmul_fma_3_taken++;
		return nextafter_3_reg_e(c);
	}

	return c;
}

// ties to even
double fdiv_0(double a, double b) {
	return a / b;
}

// toward negative infinity
// b can never be inf; inf / inf = nan
// b can be 0.0; fin / 0.0 = inf; fdiv instructions read from memory, possible zero
double fdiv_1(double a, double b) {
	double c = a / b;

	// FAIL: [fdiv fprc(1)] truth(inf, 0x1.efb83b51bc60bp+810) == inf != op(...) == 0x1.fffffffffffffp+1023
	if (!isfinite(c) && isfinite(a) && isfinite(b)) {
		return nextafter_1_reg_e(c);
	}
	
	int expa, expb;
	double a_scaled = frexp_reg_e_nozero(a, &expa);
	double b_scaled = frexp_reg_e_nozero(b, &expb);
	double c_scaled = ldexp_reg_e_nozero_noinf(c, -expa + expb);
	double res = mul_residue(c_scaled, b_scaled, a_scaled);

	fdiv_1_visit++;
	if (-res < 0.0) {
		fdiv_1_taken++;
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
// b can never be inf; inf / inf = nan
// b can be 0.0; fin / 0.0 = inf
double fdiv_2(double a, double b) {
	double c = a / b;

	if (!isfinite(c) && isfinite(a) && isfinite(b)) {
		return nextafter_2_reg_e(c);
	}
	
	int expa, expb;
	double a_scaled = frexp_reg_e_nozero(a, &expa);
	double b_scaled = frexp_reg_e_nozero(b, &expb);
	double c_scaled = ldexp_reg_e_nozero_noinf(c, -expa + expb);
	double res = mul_residue(c_scaled, b_scaled, a_scaled);

	fdiv_2_visit++;
	if (-res > 0.0) {
		fdiv_2_taken++;
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
// b can never be inf; inf / inf = nan
// b can be 0.0; fin / 0.0 = inf
double fdiv_3(double a, double b) {
	double c = a / b;

	if (!isfinite(c) && isfinite(a) && isfinite(b)) {
		return nextafter_3_reg_e(c);
	}
	
	int expa, expb;
	double a_scaled = frexp_reg_e_nozero(a, &expa);
	double b_scaled = frexp_reg_e_nozero(b, &expb);
	double c_scaled = ldexp_reg_e_nozero_noinf(c, -expa + expb);
	double res = mul_residue(c_scaled, b_scaled, a_scaled);

	fdiv_3_visit++;
	if (-res < 0.0) {
		fdiv_3_taken++;
		return nextafter_3_reg_e(c);
	}

	return c;
}

// toward negative infinity
// b can never be inf; inf / inf = nan
// b can be 0.0; fin / 0.0 = inf; fdiv instructions read from memory, possible zero
double fdiv_fma_1(double a, double b) {
	double c = a / b;
	double res = mul_residue_fma(c, b, a);

	fdiv_fma_1_visit++;
	if (-res < 0.0) {
		fdiv_fma_1_taken++;
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward negative infinity
// b can never be inf; inf / inf = nan
// b can be 0.0; fin / 0.0 = inf; fdiv instructions read from memory, possible zero
double fdiv_fma_2(double a, double b) {
	double c = a / b;
	double res = mul_residue_fma(c, b, a);

	fdiv_fma_2_visit++;
	if (-res > 0.0) {
		fdiv_fma_2_taken++;
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
// b can never be inf; inf / inf = nan
// b can be 0.0; fin / 0.0 = inf; fdiv instructions read from memory, possible zero
double fdiv_fma_3(double a, double b) {
	double c = a / b;
	double res = mul_residue_fma(c, b, a);

	fdiv_fma_3_visit++;
	if (-res < 0.0) {
		fdiv_fma_3_taken++;
		return nextafter_3_reg_e(c);
	}

	return c;
}

// ties to even
double fsqrt_0(double a) {
	return sqrt(a);
}

// toward negative infinity
double fsqrt_1(double a) {
	if (!isfinite(a)) {
		return a;
	}

	double c = sqrt(a);
	int expc;
	double c_scaled = frexp_reg_e_nozero_noinf(c, &expc);
	double a_scaled = ldexp_reg_e_nozero_noinf(a, -2 * expc);
	double res = mul_residue(c_scaled, c_scaled, a_scaled);

	fsqrt_1_visit++;
	if (-res < 0.0) {
		fsqrt_1_taken++;
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
double fsqrt_2(double a) {
	if (!isfinite(a)) {
		return a;
	}
	
	double c = sqrt(a);
	int expc;
	double c_scaled = frexp_reg_e_nozero_noinf(c, &expc);
	double a_scaled = ldexp_reg_e_nozero_noinf(a, -2 * expc);
	double res = mul_residue(c_scaled, c_scaled, a_scaled);

	fsqrt_2_visit++;
	if (-res > 0.0) {
		fsqrt_2_taken++;
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
double fsqrt_3(double a) {
	if (!isfinite(a)) {
		return a;
	}
	
	double c = sqrt(a);
	int expc;
	double c_scaled = frexp_reg_e_nozero_noinf(c, &expc);
	double a_scaled = ldexp_reg_e_nozero_noinf(a, -2 * expc);
	double res = mul_residue(c_scaled, c_scaled, a_scaled);

	fsqrt_3_visit++;
	if (-res < 0.0) {
		fsqrt_3_taken++;
		return nextafter_3_reg_e(c);
	}

	return c;
}

// toward negative infinity
double fsqrt_fma_1(double a) {
	if (!isfinite(a)) {
		return a;
	}

	double c = sqrt(a);
	double res = mul_residue_fma(c, c, a);

	fsqrt_fma_1_visit++;
	if (-res < 0.0) {
		fsqrt_fma_1_taken++;
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
double fsqrt_fma_2(double a) {
	if (!isfinite(a)) {
		return a;
	}

	double c = sqrt(a);
	double res = mul_residue_fma(c, c, a);

	fsqrt_fma_2_visit++;
	if (-res > 0.0) {
		fsqrt_fma_2_taken++;
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
double fsqrt_fma_3(double a) {
	if (!isfinite(a)) {
		return a;
	}

	double c = sqrt(a);
	double res = mul_residue_fma(c, c, a);

	fsqrt_fma_3_visit++;
	if (-res < 0.0) {
		fsqrt_fma_3_taken++;
		return nextafter_3_reg_e(c);
	}

	return c;
}
