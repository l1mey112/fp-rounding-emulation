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

double pcgd_random_r(pcg32_random_t *rng) {
	union {
		uint64_t u;
		double d;
	} u;
	u.u = pcg64_random_r(rng);
	return u.d;
}

// END PCG32

bool pcgd_truth(double lo, double hi, double v) {
	return !(isnan(v) || v < lo || v > hi || fpclassify(v) == FP_SUBNORMAL);
}

double pcgd_random_some(double lo, double hi, pcg32_random_t *rng) {
	double v;

	do {
		v = pcgd_random_r(rng);
	} while (!pcgd_truth(lo, hi, v));

	return v;
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

// https://www.ccsl.carleton.ca/%7Ejamuir/rdtscpm1.pdf
// https://en.wikipedia.org/wiki/Time_Stamp_Counter
// https://github.com/xuwd1/rdtsc-notes

static volatile inline uint64_t rdtsc() {
	uint32_t lo, hi;
	asm volatile(
		"xorl %%eax, %%eax\n\t"
		"cpuid\n\t"
		"rdtsc\n\t"
		"lfence\n\t"
		: "=a"(lo), "=d"(hi)
		:
		: "%ebx", "%ecx");
	return (uint64_t)hi << 32 | lo;
}

static volatile uint32_t rdtsc_overhead() {
	uint64_t start, end;
	start = rdtsc();
	end = rdtsc();
	return end - start;
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

#define TIMEIT(v, f)                                                    \
	do {                                                                \
		uint32_t rdtsc_base = rdtsc_overhead();                         \
		uint32_t start = rdtsc();                                       \
		f;                                                              \
		uint32_t end = rdtsc();                                         \
		int32_t ncycles = (int32_t)(end - start) - rdtsc_base;          \
		running_avg_update(&running_avg[v], ncycles < 0 ? 0 : ncycles); \
	} while (0)

	// boundary conditions
	long bound_idx = 0;
	double bound[] = {
		0.0,
		-0.0,
		1.0,
		-1.0,
		INFINITY,
		-INFINITY,
	};

	long bound_size = sizeof(bound) / sizeof(bound[0]);

	for (unsigned i = 0; i < samples; i++) {
		double a, b, c, d;

	retry:
		a = pcgd_random_some(driver->lo, driver->hi, &state);
		b = pcgd_random_some(driver->lo, driver->hi, &state);
		if (bound_idx < bound_size && pcgd_truth(driver->lo, driver->hi, bound[bound_idx])) {
			a = bound[bound_idx];
		} else if (bound_idx < bound_size * 2 && pcgd_truth(driver->lo, driver->hi, bound[bound_idx - bound_size])) {
			b = bound[bound_idx - bound_size];
		} else if (bound_idx < bound_size * 3 && pcgd_truth(driver->lo, driver->hi, bound[bound_idx - bound_size * 2])) {
			a = bound[bound_idx - bound_size * 2];
			b = bound[bound_idx - bound_size * 2];
		}
		bound_idx++;

		for (int v = 1; v < 4; v++) {
			fesetround(fprc_env[v]);
			TIMEIT(0, { c = driver->fprc[0](a, b); });

			if (!pcgd_truth(driver->lo, driver->hi, c)) {
				goto retry;
			}

			fesetround(FE_TONEAREST);
			TIMEIT(v, { d = driver->fprc[v](a, b); });

			if (!(isnan(c) || isnan(d)) && c == d) {
				complete[v]++;
			} else {
				// printf("FAIL: [%s fprc(%d)] truth(%a, %a) == %a != op(...) == %a\n", driver->desc, v, a, b, c, d);
			}
		}
	}

#undef TIMEIT

	printf("%s:\n", driver->desc);
	for (int v = 0; v < 4; v++) {
		if (v == 0) {
			printf("  fprc(%d): %uclks\n", v, (unsigned)running_avg[v].avg);
		} else {
			printf("  fprc(%d): %uclks [%d/%d] x%u overhead\n", v, (unsigned)running_avg[v].avg, complete[v], samples, (unsigned)round(running_avg[v].avg / running_avg[0].avg));
		}
	}
}

static driver_t driver[] = {
	{"fadd", -3.0e+14, 3.0e+14, fadd_0, fadd_1, fadd_2, fadd_3},
	{"fsub", -3.0e+14, 3.0e+14, fsub_0, fsub_1, fsub_2, fsub_3},
	{"fmul", 1.7e-77, INFINITY, fmul_0, fmul_1, fmul_2, fmul_3},
	{"fmul (fma)", 1.7e-77, INFINITY, fmul_0, fmul_fma_1, fmul_fma_2, fmul_fma_3},
	{"fdiv", 1.7e-77, INFINITY, fdiv_0, fdiv_1, fdiv_2, fdiv_3},
	{"fdiv (fma)", 1.7e-77, INFINITY, fdiv_0, fdiv_fma_1, fdiv_fma_2, fdiv_fma_3},
	{"fsqrt", 1.7e-77, INFINITY, fsqrt_0, fsqrt_1, fsqrt_2, fsqrt_3},
	{"fsqrt (fma)", 1.7e-77, INFINITY, fsqrt_0, fsqrt_fma_1, fsqrt_fma_2, fsqrt_fma_3},
};

int main() {
	for (int i = 0; i < sizeof(driver) / sizeof(driver[0]); i++) {
		bench(&driver[i], 100000);
	}
}
