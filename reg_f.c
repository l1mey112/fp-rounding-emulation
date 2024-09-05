#include "ulp.h"
#include <stdint.h>
#include <stdio.h>

// The domains of floating point operations are separated into "additive" operations,
// which use register group F and "multiplicative" operations, which use register group E.
//
// [-3.0e+14, 3.0e+14]

int64_t fadd_1_visit, fadd_2_visit, fadd_3_visit;
int64_t fadd_1_taken, fadd_2_taken, fadd_3_taken;

__attribute__((destructor)) void destructor_f() {
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

#define SIMDE_ENABLE_NATIVE_ALIASES
#include "simde/simde/wasm/simd128.h"

static inline v128_t sum_residue_simd(v128_t a, v128_t b, v128_t c) {
	v128_t delta_a = wasm_f64x2_sub(a, wasm_f64x2_sub(c, b));
	v128_t delta_b = wasm_f64x2_sub(b, wasm_f64x2_sub(c, a));
	v128_t res = wasm_f64x2_add(delta_a, delta_b);
	return res;
}

v128_t fadd_1_ssse(v128_t a, v128_t b) {
	v128_t c = wasm_f64x2_add(a, b);
	v128_t res = sum_residue_simd(a, b, c);

	// round toward negative infinity
	v128_t to_round = wasm_f64x2_lt(res, wasm_f64x2_const(0.0, 0.0));
	
	v128_t is_zero = wasm_f64x2_eq(c, wasm_f64x2_const(0.0, 0.0));

	v128_t k = wasm_i64x2_add(c, wasm_v128_or(wasm_f64x2_gt(c, wasm_f64x2_const(0.0, 0.0)), wasm_i64x2_const(1, 1)));

	k = wasm_v128_bitselect(wasm_i32x4_const(0x80000000, 1, 0x80000000, 1), k, is_zero); // is_zero ? +-minsubnormal : k

	// to_round ? k : c
	v128_t d = wasm_v128_bitselect(k, c, to_round);

	return d;
}

double fadd_1(double a, double b) {
	v128_t a_ssse = wasm_f64x2_splat(a);
	v128_t b_ssse = wasm_f64x2_splat(b);
	v128_t result = fadd_1_ssse(a_ssse, b_ssse);
	double c = wasm_f64x2_extract_lane(result, 0);
	return c;
}

// toward negative infinity
/* double fadd_1(double a, double b) {
    double c = a + b;
    double res = sum_residue(a, b, c);

    fadd_1_visit++;
    if (res < 0.0) {
        fadd_1_taken++;
        return nextafter_1_reg_f(c);
    }

    return c;
} */

/* v128_t fadd_2_ssse(v128_t a, v128_t b) {
	double a_0 = wasm_f64x2_extract_lane(a, 0);
	double a_1 = wasm_f64x2_extract_lane(a, 1);
	double b_0 = wasm_f64x2_extract_lane(b, 0);
	double b_1 = wasm_f64x2_extract_lane(b, 1);

	double c_0 = a_0 + b_0;
	double c_1 = a_1 + b_1;

	double res_0 = sum_residue(a_0, b_0, c_0);
	double res_1 = sum_residue(a_1, b_1, c_1);

	if (res_0 > 0.0) {
		c_0 = nextafter_2_reg_f(c_0);
	}

	if (res_1 > 0.0) {
		c_1 = nextafter_2_reg_f(c_1);
	}

	return wasm_f64x2_make(c_0, c_1);
} */

v128_t fadd_2_ssse(v128_t a, v128_t b) {
	v128_t c = wasm_f64x2_add(a, b);
	v128_t res = sum_residue_simd(a, b, c);

	// round toward positive infinity
	v128_t to_round = wasm_f64x2_gt(res, wasm_f64x2_const(0.0, 0.0));
	
	v128_t is_zero = wasm_f64x2_eq(c, wasm_f64x2_const(0.0, 0.0));

	v128_t k = wasm_i64x2_add(c, wasm_v128_or(wasm_f64x2_lt(c, wasm_f64x2_const(0.0, 0.0)), wasm_i64x2_const(1, 1)));

	k = wasm_v128_bitselect(wasm_i32x4_const(0x00000000, 1, 0x00000000, 1), k, is_zero); // is_zero ? +-minsubnormal : k

	// to_round ? k : c
	v128_t d = wasm_v128_bitselect(k, c, to_round);

	return d;
}

double fadd_2(double a, double b) {
	v128_t a_ssse = wasm_f64x2_splat(a);
	v128_t b_ssse = wasm_f64x2_splat(b);
	v128_t result = fadd_2_ssse(a_ssse, b_ssse);
	double c = wasm_f64x2_extract_lane(result, 0);
	return c;
}

// toward positive infinity
/* double fadd_2(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);

	if (res > 0.0) {
		return nextafter_2_reg_f(c);
	}

	return c;
} */

v128_t fadd_3_ssse(v128_t a, v128_t b) {
	v128_t c = wasm_f64x2_add(a, b);
	v128_t res = sum_residue_simd(a, b, c);

	// round toward zero
	// if (res > 0.0 && c < 0.0) || (res < 0.0 && c > 0.0)
	v128_t to_round = wasm_v128_or(
		wasm_v128_and(wasm_f64x2_gt(res, wasm_f64x2_const(0.0, 0.0)), wasm_f64x2_lt(c, wasm_f64x2_const(0.0, 0.0))),
		wasm_v128_and(wasm_f64x2_lt(res, wasm_f64x2_const(0.0, 0.0)), wasm_f64x2_gt(c, wasm_f64x2_const(0.0, 0.0)))
	);

	v128_t is_zero = wasm_f64x2_eq(c, wasm_f64x2_const(0.0, 0.0));

	v128_t k = wasm_i64x2_sub(c, wasm_i64x2_const(1, 1));

	k = wasm_v128_bitselect(wasm_f64x2_const(0.0, 0.0), k, is_zero); // is_zero ? 0 : k

	// to_round ? k : c
	v128_t d = wasm_v128_bitselect(k, c, to_round);

	return d;
}

double fadd_3(double a, double b) {
	v128_t a_ssse = wasm_f64x2_splat(a);
	v128_t b_ssse = wasm_f64x2_splat(b);
	v128_t result = fadd_3_ssse(a_ssse, b_ssse);
	double c = wasm_f64x2_extract_lane(result, 0);
	return c;
}

// toward zero
/* double fadd_3(double a, double b) {
	double c = a + b;
	double res = sum_residue(a, b, c);

	fadd_3_visit++;
	if ((res > 0.0 && c < 0.0) || (res < 0.0 && c > 0.0)) {
		fadd_3_taken++;
		return nextafter_3_reg_f(c);
	}

	return c;
} */

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
