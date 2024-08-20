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
	if (!isfinite(c)) {
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
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
double fmul_2(double a, double b) {
	double c = a * b;
	if (!isfinite(c)) {
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
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
double fmul_3(double a, double b) {
	double c = a * b;
	if (!isfinite(c)) {
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
		return nextafter_3_reg_e(c);
	}

	return c;
}

// toward negative infinity
double fmul_fma_1(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	if (res < 0.0) {
		return nextafter_1_reg_e(c);
	}

	return c;
}

// toward positive infinity
double fmul_fma_2(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	if (res > 0.0) {
		return nextafter_2_reg_e(c);
	}

	return c;
}

// toward zero
double fmul_fma_3(double a, double b) {
	double c = a * b;
	double res = mul_residue_fma(a, b, c);

	if (res < 0.0) {
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

	if (-res < 0.0) {
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

	if (-res > 0.0) {
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

	if (-res < 0.0) {
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

	if (-res < 0.0) {
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

	if (-res > 0.0) {
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

	if (-res < 0.0) {
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

	if (-res < 0.0) {
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

	if (-res > 0.0) {
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

	if (-res < 0.0) {
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

	if (-res < 0.0) {
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

	if (-res > 0.0) {
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

	if (-res < 0.0) {
		return nextafter_3_reg_e(c);
	}

	return c;
}
