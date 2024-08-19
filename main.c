#include <assert.h>
#include <fenv.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <time.h>

#include "reg.h"

typedef struct driver_t driver_t;

struct driver_t {
	const char *desc;
	double lo;
	double hi;
	double (*fprc[4])(
		double, double); // probably safe to call one arg function with two args
};

// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct {
	uint64_t state;
	uint64_t inc;
} pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t *rng) {
	uint64_t oldstate = rng->state;
	// Advance internal state
	rng->state = oldstate * 6364136223846793005ULL + (rng->inc | 1);
	// Calculate output function (XSH RR), uses old state for max ILP
	uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
	uint32_t rot = oldstate >> 59u;
	return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

uint64_t pcg64_random_r(pcg32_random_t *rng) {
	uint64_t low32 = pcg32_random_r(rng);
	uint64_t high32 = pcg32_random_r(rng);
	return low32 | (high32 << 32);
}

// END PCG32

// http://marc-b-reynolds.github.io/math/2020/06/16/UniformFloat.html
// (0, 1]
double pcgd_random(pcg32_random_t *rng) {
	double uni = (double)(pcg64_random_r(rng) >> (64 - 53));
	return fma(uni, 0x1p-53, 0x1p-53);
}

double pcgd_uniform(double a, double b, pcg32_random_t *state) {
	double st = pcgd_random(state);

	if (st == 1.0) {
		return b;
	}

	if (isinf(b)) {
		// set to largest finite number
		b = 0x1.fffffffffffffp1023;
	}

	return a + (b - a) * st;
}

typedef struct running_avg_t running_avg_t;

struct running_avg_t {
	double avg;
	unsigned count;
};

static void running_avg_update(running_avg_t *avg, double value) {
	if (value < 0.0) {
		return;
	}
	avg->count++;
	double a = 1.0 / avg->count;
	double b = 1.0 - a;
	avg->avg = a * value + b * avg->avg;
}

void bench(driver_t *driver, unsigned samples) {
	pcg32_random_t state = {
		// .state = time(NULL),
	};

	const int fprc_env[] = {
		FE_TONEAREST,
		FE_DOWNWARD,
		FE_UPWARD,
		FE_TOWARDZERO,
	};

	running_avg_t running_avg[4] = {};
	unsigned complete[4] = {};

#define TIMEIT(v, f)                                                                \
	do {                                                                            \
		struct timespec start, end;                                                 \
		clock_gettime(CLOCK_MONOTONIC_RAW, &start);                                 \
		asm volatile("" ::: "memory");                                              \
		f;                                                                          \
		asm volatile("" ::: "memory");                                              \
		clock_gettime(CLOCK_MONOTONIC_RAW, &end);                                   \
		running_avg_update(&running_avg[v], (double)(end.tv_nsec - start.tv_nsec)); \
	} while (0)

	for (unsigned i = 0; i < samples; i++) {
		double a, b, c, d;

		a = pcgd_uniform(driver->lo, driver->hi, &state);
		b = pcgd_uniform(driver->lo, driver->hi, &state);

		for (int v = 1; v < 4; v++) {
			fesetround(fprc_env[v]);
			TIMEIT(0, { c = driver->fprc[0](a, b); });

			assert(!(isnan(a) || fpclassify(a) == FP_SUBNORMAL || isnan(b) || fpclassify(b) == FP_SUBNORMAL || isnan(c) || fpclassify(c) == FP_SUBNORMAL
				|| a <= driver->lo || a >= driver->hi || b <= driver->lo || b >= driver->hi));

			fesetround(FE_TONEAREST);
			TIMEIT(v, { d = driver->fprc[v](a, b); });

			if (!(isnan(c) || isnan(d)) && c == d) {
				complete[v]++;
			} else {
				// printf("FAIL: [%s] truth(%a, %a) == %a != op(...) == %a\n", driver->desc, a, b, c, d);
			}
		}
	}

#undef TIMEIT

	printf("%s:\n", driver->desc);
	for (int v = 0; v < 4; v++) {
		if (v == 0) {
			printf("  fprc(%d): %uns\n", v, (unsigned)running_avg[v].avg);
		} else {
			printf("  fprc(%d): %uns [%d/%d]\n", v, (unsigned)running_avg[v].avg, complete[v], samples);
		}
	}
}

static driver_t driver[] = {
	{"fadd", -3.0e+14, 3.0e+14, fadd_0, fadd_1, fadd_2, fadd_3},
	{"fsub", -3.0e+14, 3.0e+14, fsub_0, fsub_1, fsub_2, fsub_3},
	{"fmul", 1.7e-77, INFINITY, fmul_0, fmul_1, fmul_2, fmul_3},
	{"fmul (fma)", 1.7e-77, INFINITY, fmul_0, fmul_fma_1, fmul_fma_2, fmul_fma_3},
};

int main() {
	for (int i = 0; i < sizeof(driver) / sizeof(driver[0]); i++) {
		bench(&driver[i], 1000000);
	}
}
