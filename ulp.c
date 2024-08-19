#define __HI(x) *(1+(int*)&x)
#define __LO(x) *(int*)&x
#define __HIp(x) *(1+(int*)x)
#define __LOp(x) *(int*)x

// assume always for all numbers, no subnormals and no NaNs

// toward negative infinity
// x cannot be inf, cannot return inf
double nextafter_1_noinf(double x) {
	int	hx,hy,ix;
	unsigned lx,ly;

	hx = __HI(x);           /* high word of x */
	lx = __LO(x);           /* low  word of x */
	ix = hx&0x7fffffff;     /* |x| */
	
	// y = -inf

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
// x cannot be inf, cannot return inf
double nextafter_2_noinf(double x) {
	int	hx,hy,ix;
	unsigned lx,ly;

	hx = __HI(x);           /* high word of x */
	lx = __LO(x);           /* low  word of x */
	ix = hx&0x7fffffff;     /* |x| */
	
	// y = inf

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
// x cannot be inf, cannot return inf
double nextafter_3_noinf(double x) {
	int	hx,hy,ix;
	unsigned lx,ly;

	hx = __HI(x);           /* high word of x */
	lx = __LO(x);           /* low  word of x */
	ix = hx&0x7fffffff;     /* |x| */
	
	// y = 0

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
