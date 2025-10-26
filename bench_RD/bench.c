#include "bench.h"
#include <stdio.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>

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

/* bool pcgd_truth(double lo, double hi, double v) {
	return !(isnan(v) || v < lo || v > hi || fpclassify(v) == FP_SUBNORMAL);
} */

double pcgd_random_normal_in_range(double lo, double hi, bool is_F, pcg32_random_t *rng) {
	if (is_F) {
		// Use the lowest bit of a random number for a 50/50 choice.
		uint32_t choice = pcg32_random_r(rng) & 1;
		if (choice == 1) {
			// Flip the interval to [-hi, -lo]
			double tmp = lo;
			lo = -hi;
			hi = -tmp;
		}
	}
	
	const uint64_t MANTISSA_MASK = 0x000FFFFFFFFFFFFFULL;
	const uint64_t ONE_POINT_ZERO_BITS = 0x3FF0000000000000ULL;

	union {
		uint64_t u;
		double d;
	} u;

	uint64_t random_bits = pcg64_random_r(rng);

	u.u = ONE_POINT_ZERO_BITS | (random_bits & MANTISSA_MASK);
	return lo + (u.d - 1.0) * (hi - lo);
}

typedef struct running_avg_t running_avg_t;

/* struct running_avg_t {
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
} */

#include "timer.h"

typedef v128_t op(v128_t, v128_t);

/* artificial use of all of memory */
# define BENCH_CLOBBER() asm volatile("":::"memory")
/* artificial dependency of x on all of memory and all of memory on x */
# define BENCH_VOLATILE(x) asm volatile("" : "+g"(x) : "g"(x) : "memory")
# define BENCH_VOLATILE_REG(x) asm volatile("" : "+r"(x) : "r"(x) : "memory")
# define BENCH_VOLATILE_MEM(x) asm volatile("" : "+m"(x) : "m"(x) : "memory")

typedef struct {
	v128_t a;
	v128_t b;
} v128_pair_t;

__attribute__((noinline))
void bench(const char *name, unsigned samples, v128_pair_t sample_pairs[samples], double lo, double hi, bool is_F, op fn) {
	// fprc = 1 (FE_DOWNWARD)
	uint64_t timing_overhead = measure_overhead();

	#define UNROLL 4

	uint64_t temp = start_timer();
	for (unsigned i = 0; i < samples; i++) {
		v128_pair_t pair = sample_pairs[i];
		v128_t c;

		c = fn(pair.a, pair.b);
		BENCH_VOLATILE_MEM(c);
		c = fn(pair.a, pair.b);
		BENCH_VOLATILE_MEM(c);
		c = fn(pair.a, pair.b);
		BENCH_VOLATILE_MEM(c);
		c = fn(pair.a, pair.b);
		BENCH_VOLATILE_MEM(c);
	}
	temp = end_timer() - temp;

	uint64_t clks = (temp / samples);
	if (clks < timing_overhead) {
		clks = 0;
	} else {
		clks -= timing_overhead;
	}
	clks /= UNROLL;

	printf("  fprc(1): \"%s\" %luclks\n", name, clks);

#undef UNROLL
}

__attribute__((noinline))
v128_t hard_fadd_1(v128_t dest, v128_t src) {
	return wasm_f64x2_add(dest, src);
}

#define SAMPLES 100000000

v128_pair_t sample_pairs[SAMPLES];

int main(void) {
	// F (1e+4, 1e+14) and (-1e+14, -1e+4)
	// E {0} and (1e-40, inf)

	pcg32_random_t state = {
		.state = time(NULL),
	};

	for (unsigned i = 0; i < SAMPLES; i++) {
		double a0, a1, b0, b1;

		a0 = pcgd_random_normal_in_range(1e+4, 1e+14, true, &state);
		a1 = pcgd_random_normal_in_range(1e+4, 1e+14, true, &state);
		b0 = pcgd_random_normal_in_range(1e+4, 1e+14, true, &state);
		b1 = pcgd_random_normal_in_range(1e+4, 1e+14, true, &state);

		sample_pairs[i].a = wasm_f64x2_make(a0, a1);
		sample_pairs[i].b = wasm_f64x2_make(b0, b1);
	}

	bench("semi_fadd_1", SAMPLES, sample_pairs, 1e+4, 1e+14, true, semi_fadd_1);
	bench("soft_fadd_1", SAMPLES, sample_pairs, 1e+4, 1e+14, true, soft_fadd_1);
	
	fesetround(FE_DOWNWARD);
	bench("hard_fadd_1", SAMPLES, sample_pairs, 1e+4, 1e+14, true, hard_fadd_1);
	fesetround(FE_TONEAREST);
}