#include <stdint.h>
#include <fenv.h>
#include <stdbool.h>
#include <stdio.h>
#include <math.h>

// FE_TONEAREST  - round ties to even
// FE_DOWNWARD   - round toward negative infinity
// FE_UPWARD     - round toward positive infinity
// FE_TOWARDZERO - round toward zero

static uint64_t state = 1;

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

double u64_double(uint64_t x) {
	union {
		uint64_t u;
		double d;
	} u;
	u.u = x;
	return u.d;
}

uint64_t double_u64(double x) {
	union {
		uint64_t u;
		double d;
	} u;
	u.d = x;
	return u.u;
}

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng) {
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

double pcgd_random_r(pcg32_random_t* rng) {
	uint32_t low32 = pcg32_random_r(rng);
	uint32_t high32 = pcg32_random_r(rng);
	return u64_double((uint64_t)low32 | ((uint64_t)high32 << 32));
}

// [0, 1)
double pcgd01_random_r(pcg32_random_t* rng) {
	return (double)pcg32_random_r(rng) / (double)UINT32_MAX;
}

double pc_uniform(double a, double b, pcg32_random_t *state) {
	return a + (b - a) * pcgd01_random_r(state) / UINT32_MAX;
}

void compare(const char *desc, int fe, int iterations, double lo, double hi, double reference(double, double), double op(double, double)) {
	pcg32_random_t state = {
		// .state = time(NULL),
	};

	// ±3.0e+14 ∪ [1.7E-77, ∞) => (-3.0e+14, ∞)

	int complete = 0;
	for (int i = 0; i < iterations; i++) {
		double a;
		double b;
		double c;

		do {
			fesetround(FE_TONEAREST);
			a = pcgd_random_r(&state);
			b = pcgd_random_r(&state);
			fesetround(fe);
			c = reference(a, b);
		} while (isnan(a) || fpclassify(a) == FP_SUBNORMAL || isnan(b) || fpclassify(b) == FP_SUBNORMAL || isnan(c) || fpclassify(c) == FP_SUBNORMAL
			|| a <= lo || a >= hi || b <= lo || b >= hi || c <= lo || c >= hi
		);

		fesetround(FE_TONEAREST);
		double d = op(a, b);

		// check NaNs and if the results are different
		if (!(isnan(c) || isnan(d)) && c == d) {
			complete++;
		} else {
			printf("FAIL: [%s] truth(%a, %a) == %a != op(...) == %a\n", desc, a, b, c, d);
		}
	}

	printf("%s: %d/%d\n", desc, complete, iterations);
}

void compare_unary(const char *desc, int fe, int iterations, double lo, double hi, double reference(double), double op(double)) {
	pcg32_random_t state = {
		// .state = time(NULL),
	};

	int complete = 0;
	for (int i = 0; i < iterations; i++) {
		double a;
		double c;

		do {
			fesetround(FE_TONEAREST);
			a = pcgd_random_r(&state);
			fesetround(fe);
			c = reference(a);
		} while (isnan(a) || fpclassify(a) == FP_SUBNORMAL || isnan(c) || fpclassify(c) == FP_SUBNORMAL
			|| a <= lo || a >= hi || c <= lo || c >= hi
		);

		fesetround(FE_TONEAREST);
		double d = op(a);

		// check NaNs and if the results are different
		if (!(isnan(c) || isnan(d)) && c == d) {
			complete++;
		} else {
			printf("FAIL: [%s] truth(...) == %a != op(...) == %a\n", desc, c, d);
		}
	}

	printf("%s: %d/%d\n", desc, complete, iterations);
}


#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x

double nextafter_better(double x, double y) {
	int	hx,hy,ix,iy;
	unsigned lx,ly;

	hx = __HI(x);		/* high word of x */
	lx = __LO(x);		/* low  word of x */
	hy = __HI(y);		/* high word of y */
	ly = __LO(y);		/* low  word of y */
	ix = hx&0x7fffffff;		/* |x| */
	iy = hy&0x7fffffff;		/* |y| */

	if(hx>=0) {				/* x > 0 */
	    if(hx>hy||((hx==hy)&&(lx>ly))) {	/* x > y, x -= ulp */
			if(lx==0) hx -= 1;
			lx -= 1;
	    } else {				/* x < y, x += ulp */
			lx += 1;
			if(lx==0) hx += 1;
	    }
	} else {				/* x < 0 */
	    if(hy>=0||hx>hy||((hx==hy)&&(lx>ly))){/* x < y, x -= ulp */
			if(lx==0) hx -= 1;
			lx -= 1;
	    } else {				/* x > y, x += ulp */
			lx += 1;
			if(lx==0) hx += 1;
	    }
	}
	hy = hx&0x7ff00000;
	if(hy>=0x7ff00000) return x+x;	/* overflow  */
	if(hy<0x00100000) {		/* underflow */
	    y = x*x;
	    if(y!=x) {		/* raise underflow flag */
			__HI(y) = hx; __LO(y) = lx;
			return y;
	    }
	}
	__HI(x) = hx; __LO(x) = lx;
	return x;
}

double round_toward(double c, double res, int fe) {
	double direction;
	switch (fe) {
		case FE_DOWNWARD: direction = -INFINITY; break;
		case FE_UPWARD: direction = INFINITY; break;
		case FE_TOWARDZERO: direction = 0.0; break;
	}

	if (res > 0.0 && c < direction) return nextafter_better(c, direction);
	if (res < 0.0 && c > direction) return nextafter_better(c, direction);
	return c;
}

double upper_half(double x) {
	double secator = 134217729.0;
	double p = x * secator;
	return p + (x - p);
}

double emulated_nfma(double a, double b, double c) {
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

double mul_residue(double a, double b, double c) {
	return emulated_nfma(a, b, -c);
}

double sum_residue(double a, double b, double c) {
	double delta_a = a - (c - b);
	double delta_b = b - (c - a);
	double res = delta_a + delta_b;
	return res;
}

double add_reference(double a, double b) {
	return a + b;
}

double add_impl(double a, double b, int fe) {
	double c = a + b;
	if (!isfinite(c) && isfinite(a) && isfinite(b) && b != 0.0) {
		return round_toward(c, -c, fe);
	}
	double res = sum_residue(a, b, c);
	return round_toward(c, res, fe);
}

double add_fe_upward(double a, double b) {
	return add_impl(a, b, FE_UPWARD);
}

double add_fe_downward(double a, double b) {
	return add_impl(a, b, FE_DOWNWARD);
}

double add_fe_towardzero(double a, double b) {
	return add_impl(a, b, FE_TOWARDZERO);
}

double sub_reference(double a, double b) {
	return a - b;
}

double sub_fe_upward(double a, double b) {
	return add_fe_upward(a, -b);
}

double sub_fe_downward(double a, double b) {
	return add_fe_downward(a, -b);
}

double sub_fe_towardzero(double a, double b) {
	return add_fe_towardzero(a, -b);
}

double mul_reference(double a, double b) {
	return a * b;
}

double mul_impl(double a, double b, int fe) {
	double c = a * b;
	if (!isfinite(c) && isfinite(a) && isfinite(b) && b != 0.0) {
		return round_toward(c, -c, fe);
	}
	int expa, expb;
	double a_scaled = frexp(a, &expa);
	double b_scaled = frexp(b, &expb);
	double c_scaled = ldexp(c, -expa - expb);
	double res = mul_residue(a_scaled, b_scaled, c_scaled);
	return round_toward(c, res, fe);
}

double mul_fe_upward(double a, double b) {
	return mul_impl(a, b, FE_UPWARD);
}

double mul_fe_downward(double a, double b) {
	return mul_impl(a, b, FE_DOWNWARD);
}

double mul_fe_towardzero(double a, double b) {
	return mul_impl(a, b, FE_TOWARDZERO);
}

double div_reference(double a, double b) {
	return a / b;
}

double div_impl(double a, double b, int fe) {
	double c = a / b;
	if (!isfinite(c) && isfinite(a) && isfinite(b) && b != 0.0) {
		return round_toward(c, -c, fe);
	}
	int expa, expb;
	double a_scaled = frexp(a, &expa);
	double b_scaled = frexp(b, &expb);
	double c_scaled = ldexp(c, -expa + expb);
	double res = mul_residue(c_scaled, b_scaled, a_scaled);
	return round_toward(c, -res, fe);
}

double div_fe_upward(double a, double b) {
	return div_impl(a, b, FE_UPWARD);
}

double div_fe_downward(double a, double b) {
	return div_impl(a, b, FE_DOWNWARD);
}

double div_fe_towardzero(double a, double b) {
	return div_impl(a, b, FE_TOWARDZERO);
}

double sqrt_reference(double a) {
	return sqrt(a);
}

double sqrt_impl(double a, int fe) {
	double c = sqrt(a);
	int expc;
	double c_scaled = frexp(c, &expc);
	double a_scaled = ldexp(a, -2 * expc);
	double res = mul_residue(c_scaled, c_scaled, a_scaled);
	return round_toward(c, -res, fe);
}

double sqrt_fe_upward(double a) {
	return sqrt_impl(a, FE_UPWARD);
}

double sqrt_fe_downward(double a) {
	return sqrt_impl(a, FE_DOWNWARD);
}

double sqrt_fe_towardzero(double a) {
	return sqrt_impl(a, FE_TOWARDZERO);
}

int main() {
	int iterations = 10000000;

	compare("add(FE_TOWARDZERO)", FE_TOWARDZERO, iterations, -3.0e+14, 3.0e+14, add_reference, add_fe_towardzero);
	compare("add(FE_UPWARD)", FE_UPWARD, iterations, -3.0e+14, 3.0e+14, add_reference, add_fe_upward);
	compare("add(FE_DOWNWARD)", FE_DOWNWARD, iterations, -3.0e+14, 3.0e+14, add_reference, add_fe_downward);

	compare("sub(FE_TOWARDZERO)", FE_TOWARDZERO, iterations, -3.0e+14, 3.0e+14, sub_reference, sub_fe_towardzero);
	compare("sub(FE_UPWARD)", FE_UPWARD, iterations, -3.0e+14, 3.0e+14, sub_reference, sub_fe_upward);
	compare("sub(FE_DOWNWARD)", FE_DOWNWARD, iterations, -3.0e+14, 3.0e+14, sub_reference, sub_fe_downward);

	compare("mul(FE_TOWARDZERO)", FE_TOWARDZERO, iterations, 1.7E-77, INFINITY, mul_reference, mul_fe_towardzero);
	compare("mul(FE_UPWARD)", FE_UPWARD, iterations, 1.7E-77, INFINITY, mul_reference, mul_fe_upward);
	compare("mul(FE_DOWNWARD)", FE_DOWNWARD, iterations, 1.7E-77, INFINITY, mul_reference, mul_fe_downward);

	compare("div(FE_TOWARDZERO)", FE_TOWARDZERO, iterations, 1.7E-77, INFINITY, div_reference, div_fe_towardzero);
	compare("div(FE_UPWARD)", FE_UPWARD, iterations, 1.7E-77, INFINITY, div_reference, div_fe_upward);
	compare("div(FE_DOWNWARD)", FE_DOWNWARD, iterations, 1.7E-77, INFINITY, div_reference, div_fe_downward);

	compare_unary("sqrt(FE_TOWARDZERO)", FE_TOWARDZERO, iterations, 1.7E-77, INFINITY, sqrt_reference, sqrt_fe_towardzero);
	compare_unary("sqrt(FE_UPWARD)", FE_UPWARD, iterations, 1.7E-77, INFINITY, sqrt_reference, sqrt_fe_upward);
	compare_unary("sqrt(FE_DOWNWARD)", FE_DOWNWARD, iterations, 1.7E-77, INFINITY, sqrt_reference, sqrt_fe_downward);
}
