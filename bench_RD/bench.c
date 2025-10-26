#include "bench.h"
#include <math.h>
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
double bench(const char *name, unsigned samples, v128_pair_t sample_pairs[samples], op fn, double ref) {
	// fprc = 1 (FE_DOWNWARD)
	uint64_t timing_overhead = measure_overhead();

	v128_t c;
	#define UNROLL 16

	uint64_t temp = start_timer();
	for (unsigned i = 0; i < samples / UNROLL; i += UNROLL) {
		#define LOAD_PAIR(n) v128_pair_t pair##n = sample_pairs[i + n];
		#define THING(n) c = fn(pair##n.a, pair##n.b); BENCH_VOLATILE_MEM(c); BENCH_CLOBBER();
		
		LOAD_PAIR(0)  LOAD_PAIR(1)  LOAD_PAIR(2)  LOAD_PAIR(3)
		LOAD_PAIR(4)  LOAD_PAIR(5)  LOAD_PAIR(6)  LOAD_PAIR(7)
		LOAD_PAIR(8)  LOAD_PAIR(9)  LOAD_PAIR(10) LOAD_PAIR(11)
		LOAD_PAIR(12) LOAD_PAIR(13) LOAD_PAIR(14) LOAD_PAIR(15)

		THING(0)  THING(1)  THING(2)  THING(3)
		THING(4)  THING(5)  THING(6)  THING(7)
		THING(8)  THING(9)  THING(10) THING(11)
		THING(12) THING(13) THING(14) THING(15)
	}
	temp = end_timer() - temp;

	double ticks = (temp - timing_overhead);
	ticks /= (double)samples;
	
	BENCH_VOLATILE_MEM(c);

	if (isnan(ref)) {
		printf("  fprc(1): %8.2f ticks \"%s\"\n", ticks, name);	
	} else {
		double ref_overhead = ticks / ref;
		printf("  fprc(1): %8.2f ticks \"%s\" (x%.2f overhead)\n", ticks, name, ref_overhead);
	}

#undef UNROLL
	return ticks;
}

__attribute__((noinline))
v128_t hard_fadd_1(v128_t dest, v128_t src) {
	return wasm_f64x2_add(dest, src);
}

__attribute__((noinline))
v128_t hard_fmul_1(v128_t dest, v128_t src) {
	return wasm_f64x2_mul(dest, src);
}

#define SAMPLES 160000

v128_pair_t sample_pairs[SAMPLES];

void bench_fill(double lo, double hi, bool is_F) {
	pcg32_random_t state = {
		.state = time(NULL),
	};
	
	for (unsigned i = 0; i < SAMPLES; i++) {
		double a0, a1, b0, b1;

		a0 = pcgd_random_normal_in_range(lo, hi, is_F, &state);
		a1 = pcgd_random_normal_in_range(lo, hi, is_F, &state);
		b0 = pcgd_random_normal_in_range(lo, hi, is_F, &state);
		b1 = pcgd_random_normal_in_range(lo, hi, is_F, &state);

		sample_pairs[i].a = wasm_f64x2_make(a0, a1);
		sample_pairs[i].b = wasm_f64x2_make(b0, b1);
	}
}

int main(void) {
	// F (1e+4, 1e+14) and (-1e+14, -1e+4)
	// E {0} and (1e-40, inf)

	bench_fill(1e+4, 1e+14, true);

	double add_ref, mul_ref;

	fesetround(FE_DOWNWARD);
	add_ref = bench("hard_fadd_1", SAMPLES, sample_pairs, hard_fadd_1, NAN);
	fesetround(FE_TONEAREST);
	bench("semi_fadd_1", SAMPLES, sample_pairs, semi_fadd_1, add_ref);
	bench("soft_fadd_1", SAMPLES, sample_pairs, soft_fadd_1, add_ref);

	bench_fill(1e-40, 1e+50, false);
	printf("\n");

	fesetround(FE_DOWNWARD);
	mul_ref = bench("hard_fmul_1", SAMPLES, sample_pairs, hard_fmul_1, NAN);
	fesetround(FE_TONEAREST);
	bench("semi_fmul_1", SAMPLES, sample_pairs, semi_fmul_1, mul_ref);
	bench("semi_fmul_fma_1", SAMPLES, sample_pairs, semi_fmul_fma_1, mul_ref);
	bench("soft_fmul_1", SAMPLES, sample_pairs, soft_fmul_1, mul_ref);
}
