#include <math.h>

#include "ulp.h"

#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x

// assume always for all numbers, no subnormals and no NaNs

// toward negative infinity
double nextafter_1(double x) {
	int	hx,hy,ix;
	unsigned lx,ly;

	hx = __HI(x);           /* high word of x */
	lx = __LO(x);           /* low  word of x */
	ix = hx&0x7fffffff;     /* |x| */
	
	// y = -inf

	if (x == -INFINITY) return x; /* x=y, return x */

	if((ix|lx)==0) { /* x == 0 */
		// __HI(y) & 0x80000000 /* extract sign, rest zero */
		__HI(x) = 0x80000000; /* return +-minsubnormal */
		__LO(x) = 1;
		return x;
	}

	if (hx>=0) { /* x > 0 */
		/* x > y, x -= ulp */
		if (lx==0) hx -= 1;
		lx -= 1;
	} else { /* x < 0 */
		/* x > y, x += ulp */
		lx += 1;
		if (lx==0) hx += 1;
	}

	__HI(x) = hx; __LO(x) = lx;
	return x;
}

// toward positive infinity
double nextafter_2(double x) {
	int	hx,hy,ix;
	unsigned lx,ly;

	hx = __HI(x);           /* high word of x */
	lx = __LO(x);           /* low  word of x */
	ix = hx&0x7fffffff;     /* |x| */
	
	// y = inf

	if (x == INFINITY) return x; /* x=y, return x */

	if((ix|lx)==0) { /* x == 0 */
		// __HI(y) & 0x80000000 /* extract sign, rest zero */
		__HI(x) = 0x00000000; /* return +-minsubnormal */
		__LO(x) = 1;
		return x;
	}

	if (hx>=0) { /* x > 0 */
		/* x < y, x += ulp */
		lx += 1;
		if (lx==0) hx += 1;
	} else { /* x < 0 */
		/* x < y, x -= ulp */
		if (lx==0) hx -= 1;
		lx -= 1;
	}

	__HI(x) = hx; __LO(x) = lx;
	return x;
}

// toward zero
double nextafter_3(double x) {
	int	hx,hy,ix;
	unsigned lx,ly;

	hx = __HI(x);           /* high word of x */
	lx = __LO(x);           /* low  word of x */
	ix = hx&0x7fffffff;     /* |x| */
	
	// y = 0

	if (x == 0.0) return x; /* x=y, return x */

	if((ix|lx)==0) { /* x == 0 */
		// __HI(y) & 0x80000000 /* extract sign, rest zero */
		__HI(x) = 0x00000000; /* return +-minsubnormal */
		__LO(x) = 1;
		return x;
	}

	if (hx>=0) { /* x > 0 */
		/* x > y, x -= ulp */
		if (lx==0) hx -= 1;
		lx -= 1;
	} else { /* x < 0 */
		/* x < y, x -= ulp */
		if (lx==0) hx -= 1;
		lx -= 1;
	}

	__HI(x) = hx; __LO(x) = lx;
	return x;
}
