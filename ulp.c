#include <assert.h>
#include <math.h>

#include "ulp.h"

#define __HI(x) *(1 + (int *)&x)
#define __LO(x) *(int *)&x
#define __HIp(x) *(1 + (int *)x)
#define __LOp(x) *(int *)x

// assume always for all numbers, no subnormals and no NaNs

// toward negative infinity
double nextafter_1_reg_f(double x) {
	int hx, hy, ix;
	unsigned lx, ly;

	hx = __HI(x);         /* high word of x */
	lx = __LO(x);         /* low  word of x */
	ix = hx & 0x7fffffff; /* |x| */

	// y = -inf

	if ((ix | lx) == 0) { /* x == 0 */
		// __HI(y) & 0x80000000 /* extract sign, rest zero */
		__HI(x) = 0x80000000; /* return +-minsubnormal */
		__LO(x) = 1;
		return x;
	}

	if (hx >= 0) { /* x > 0 */
		/* x > y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	} else { /* x < 0 */
		/* x > y, x += ulp */
		lx += 1;
		if (lx == 0)
			hx += 1;
	}

	__HI(x) = hx;
	__LO(x) = lx;
	return x;
}

// toward positive infinity
double nextafter_2_reg_f(double x) {
	int hx, hy, ix;
	unsigned lx, ly;

	hx = __HI(x);         /* high word of x */
	lx = __LO(x);         /* low  word of x */
	ix = hx & 0x7fffffff; /* |x| */

	// y = inf

	if ((ix | lx) == 0) { /* x == 0 */
		// __HI(y) & 0x80000000 /* extract sign, rest zero */
		__HI(x) = 0x00000000; /* return +-minsubnormal */
		__LO(x) = 1;
		return x;
	}

	if (hx >= 0) { /* x > 0 */
		/* x < y, x += ulp */
		lx += 1;
		if (lx == 0)
			hx += 1;
	} else { /* x < 0 */
		/* x < y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	}

	__HI(x) = hx;
	__LO(x) = lx;
	return x;
}

// toward zero
double nextafter_3_reg_f(double x) {
	int hx, hy, ix;
	unsigned lx, ly;

	hx = __HI(x);         /* high word of x */
	lx = __LO(x);         /* low  word of x */
	ix = hx & 0x7fffffff; /* |x| */

	// y = 0

	if ((ix | lx) == 0) { /* x == 0 */
		return x; /* x=y, return x */
	}

	if (hx >= 0) { /* x > 0 */
		/* x > y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	} else { /* x < 0 */
		/* x < y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	}

	__HI(x) = hx;
	__LO(x) = lx;
	return x;
}

// toward negative infinity
// x cannot be zero or negative infinity
double nextafter_1_reg_e(double x) {
	int hx, hy, ix;
	unsigned lx, ly;

	hx = __HI(x);         /* high word of x */
	lx = __LO(x);         /* low  word of x */
	ix = hx & 0x7fffffff; /* |x| */

	// y = -inf

	if (hx >= 0) { /* x > 0 */
		/* x > y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	} else { /* x < 0 */
		/* x > y, x += ulp */
		lx += 1;
		if (lx == 0)
			hx += 1;
	}

	__HI(x) = hx;
	__LO(x) = lx;
	return x;
}

// toward positive infinity
// x cannot be zero
double nextafter_2_reg_e(double x) {
	int hx, hy, ix;
	unsigned lx, ly;

	hx = __HI(x);         /* high word of x */
	lx = __LO(x);         /* low  word of x */
	ix = hx & 0x7fffffff; /* |x| */

	// y = inf

	if (x == INFINITY)
		return x; /* x=y, return x */

	if (hx >= 0) { /* x > 0 */
		/* x < y, x += ulp */
		lx += 1;
		if (lx == 0)
			hx += 1;
	} else { /* x < 0 */
		/* x < y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	}

	__HI(x) = hx;
	__LO(x) = lx;
	return x;
}

// toward zero
// x cannot be zero
double nextafter_3_reg_e(double x) {
	int hx, hy, ix;
	unsigned lx, ly;

	hx = __HI(x);         /* high word of x */
	lx = __LO(x);         /* low  word of x */
	ix = hx & 0x7fffffff; /* |x| */

	// y = 0

	if (hx >= 0) { /* x > 0 */
		/* x > y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	} else { /* x < 0 */
		/* x < y, x -= ulp */
		if (lx == 0)
			hx -= 1;
		lx -= 1;
	}

	__HI(x) = hx;
	__LO(x) = lx;
	return x;
}

/*
 * for non-zero x
 *	x = frexp(arg,&exp);
 * return a double fp quantity x such that 0.5 <= |x| <1.0
 * and the corresponding binary exponent "exp". That is
 *	arg = x*2^exp.
 * If arg is inf, 0.0, or NaN, then frexp(arg,&exp) returns arg
 * with *exp=0.
 */
static const double
	two54 = 1.80143985094819840000e+16,  /* 0x43500000, 0x00000000 */
	twom54 = 5.55111512312578270212e-17, /* 0x3C900000, 0x00000000 */
	huge = 1.0e+300,
	tiny = 1.0e-300;

// assume no NaN, subnormal, inf, or zero
double frexp_reg_e_nozero_noinf(double x, int *eptr) {
	int hx, ix, lx;
	hx = __HI(x);
	ix = 0x7fffffff & hx;
	lx = __LO(x);
	*eptr = 0;
	*eptr += (ix >> 20) - 1022;
	hx = (hx & 0x800fffff) | 0x3fe00000;
	__HI(x) = hx;
	return x;
}

// assume no NaN, subnormal, or zero
double frexp_reg_e_nozero(double x, int *eptr) {
	int hx, ix, lx;
	hx = __HI(x);
	ix = 0x7fffffff & hx;
	lx = __LO(x);
	*eptr = 0;
	if(ix>=0x7ff00000) return x; /* inf,nan */
	*eptr += (ix >> 20) - 1022;
	hx = (hx & 0x800fffff) | 0x3fe00000;
	__HI(x) = hx;
	return x;
}

// assume no NaN, subnormal, inf, or zero
double ldexp_reg_e_nozero_noinf(double x, int n) {
	int hx, k;
	hx = __HI(x);
	k = (hx & 0x7ff00000) >> 20; /* extract exponent */

	k += n; /* add the exponent */

	/* The result must be a normal number, so k must stay within the range [1, 0x7fe]. */
	__HI(x) = (hx & 0x800fffff) | (k << 20); /* set the new exponent */

	return x;
}
