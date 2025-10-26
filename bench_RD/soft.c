#include <stdbool.h>
#include <stdint.h>
#include "bench.h"

typedef struct { uint64_t v; } float64_t;
union ui64_f64 { uint64_t ui; float64_t f; };

#define signF64UI( a ) ((bool) ((uint64_t) (a)>>63))
#define expF64UI( a ) ((int_fast16_t) ((a)>>52) & 0x7FF)
#define fracF64UI( a ) ((a) & UINT64_C( 0x000FFFFFFFFFFFFF ))
#define packToF64UI( sign, exp, sig ) ((uint64_t) (((uint_fast64_t) (sign)<<63) + ((uint_fast64_t) (exp)<<52) + (sig)))

#define isNaNF64UI( a ) (((~(a) & UINT64_C( 0x7FF0000000000000 )) == 0) && ((a) & UINT64_C( 0x000FFFFFFFFFFFFF )))

static inline uint64_t softfloat_shiftRightJam64( uint64_t a, uint_fast32_t dist )
{
    return
        (dist < 63) ? a>>dist | ((uint64_t) (a<<(-dist & 63)) != 0) : (a != 0);
}

enum {
    softfloat_round_near_even   = 0,
    softfloat_round_minMag      = 1,
    softfloat_round_min         = 2,
    softfloat_round_max         = 3,
    softfloat_round_near_maxMag = 4,
    softfloat_round_odd         = 6
};
//uint_fast8_t softfloat_roundingMode = softfloat_round_min;
#define softfloat_roundingMode softfloat_round_min
/*

softfloat_round_near_even 	round to nearest, with ties to even
softfloat_round_near_maxMag   	round to nearest, with ties to maximum magnitude (away from zero)
softfloat_round_minMag 	round to minimum magnitude (toward zero)
softfloat_round_min 	round to minimum (down)
softfloat_round_max 	round to maximum (up)
softfloat_round_odd 	round to odd (jamming), if supported by the SoftFloat port

*/

static float64_t
 softfloat_roundPackToF64( bool sign, int_fast16_t exp, uint_fast64_t sig )
{
    uint_fast16_t roundIncrement, roundBits;
    bool isTiny;
    uint_fast64_t uiZ;
    union ui64_f64 uZ;

    roundIncrement = sign ? 0x3F : 0;
    roundBits = sig & 0x3FF;
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    // overflow to infinity and underflow to subnormals
    /* if ( 0x7FD <= (uint16_t) exp ) {
        if ( exp < 0 ) {
            sig = softfloat_shiftRightJam64( sig, -exp );
            exp = 0;
            roundBits = sig & 0x3FF;
        } else if (
            (0x7FD < exp)
                || (UINT64_C( 0x8000000000000000 ) <= sig + roundIncrement)
        ) {
            uiZ = packToF64UI( sign, 0x7FF, 0 ) - ! roundIncrement;
            goto uiZ;
        }
    } */
    sig = (sig + roundIncrement)>>10;
    sig &= ~(uint_fast64_t) (! (roundBits ^ 0x200));

    // can never be called from softfloat_addMagsF64 as 0 + 0 == 0 is the only value that equals zero
    // see (1*)
    //if ( ! sig ) exp = 0;

    uiZ = packToF64UI( sign, exp, sig );
    uZ.ui = uiZ;
    return uZ.f;
}

// performs |A| + |B|
// (1*) |A| + |B| != 0 for nonzero A, B. so under these circumstances |A| + |B| == 0 <=> |A| == 0, |B| == 0
static float64_t
 softfloat_addMagsF64( uint_fast64_t uiA, uint_fast64_t uiB, bool signZ )
{
    int_fast16_t expA;
    uint_fast64_t sigA;
    int_fast16_t expB;
    uint_fast64_t sigB;
    int_fast16_t expDiff;
    uint_fast64_t uiZ;
    int_fast16_t expZ;
    uint_fast64_t sigZ;
    union ui64_f64 uZ;

    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    expA = expF64UI( uiA );
    sigA = fracF64UI( uiA );
    expB = expF64UI( uiB );
    sigB = fracF64UI( uiB );
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
    expDiff = expA - expB;
    if ( ! expDiff ) {
        /*--------------------------------------------------------------------
        *--------------------------------------------------------------------*/
        // expA == expB == 0, the numbers must be zero, no subnormals
        if ( ! expA ) {
            uiZ = uiA + sigB;
            goto uiZ;
        }
        // infinities don't happen
        /* if ( expA == 0x7FF ) {
            uiZ = uiA;
            goto uiZ;
        } */
        expZ = expA;
        sigZ = UINT64_C( 0x0020000000000000 ) + sigA + sigB;
        sigZ <<= 9;
    } else {
        sigA <<= 9;
        sigB <<= 9;

        if ( expDiff < 0 ) {
            // infinities don't happen
            /* if ( expB == 0x7FF ) {
                uiZ = packToF64UI( signZ, 0x7FF, 0 );
                goto uiZ;
            } */
            expZ = expB;
            if ( expA ) {
                sigA += UINT64_C( 0x2000000000000000 );
            } else {
                // if expA == 0, then sigA == 0. no subnormals
                //sigA <<= 1;
            }
            sigA = softfloat_shiftRightJam64( sigA, -expDiff );
        } else {
            // infinities don't happen
            /* if ( expA == 0x7FF ) {
                uiZ = uiA;
                goto uiZ;
            } */
            expZ = expA;
            if ( expB ) {
                sigB += UINT64_C( 0x2000000000000000 );
            } else {
                // if expB == 0, then sigB == 0. no subnormals
                //sigB <<= 1;
            }
            sigB = softfloat_shiftRightJam64( sigB, expDiff );
        }
        sigZ = UINT64_C( 0x2000000000000000 ) + sigA + sigB;
        if ( sigZ < UINT64_C( 0x4000000000000000 ) ) {
            --expZ;
            sigZ <<= 1;
        }
    }
    return softfloat_roundPackToF64( signZ, expZ, sigZ );
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
 uiZ:
    uZ.ui = uiZ;
    return uZ.f;

}

static inline uint_fast8_t softfloat_countLeadingZeros64( uint64_t a ) {
    return __builtin_clz(a);
}

static float64_t
 softfloat_normRoundPackToF64( bool sign, int_fast16_t exp, uint_fast64_t sig )
{
    int_fast8_t shiftDist;
    union ui64_f64 uZ;

    shiftDist = softfloat_countLeadingZeros64( sig ) - 1;
    exp -= shiftDist;
    // no infinities
    if ( (10 <= shiftDist) /* && ((unsigned int) exp < 0x7FD) */ ) {
        // sig is nonzero
        uZ.ui = packToF64UI( sign, exp /* sig ? exp : 0 */, sig<<(shiftDist - 10) );
        return uZ.f;
    } else {
        return softfloat_roundPackToF64( sign, exp, sig<<shiftDist );
    }

}

static float64_t
 softfloat_subMagsF64( uint_fast64_t uiA, uint_fast64_t uiB, bool signZ )
{
    int_fast16_t expA;
    uint_fast64_t sigA;
    int_fast16_t expB;
    uint_fast64_t sigB;
    int_fast16_t expDiff;
    uint_fast64_t uiZ;
    int_fast64_t sigDiff;
    int_fast8_t shiftDist;
    int_fast16_t expZ;
    uint_fast64_t sigZ;
    union ui64_f64 uZ;

    expA = expF64UI( uiA );
    sigA = fracF64UI( uiA );
    expB = expF64UI( uiB );
    sigB = fracF64UI( uiB );

    expDiff = expA - expB;

    if ( ! expDiff ) {
        // no NaN no infnities
        /* if ( expA == 0x7FF ) {
            if ( sigA | sigB ) goto propagateNaN;
            softfloat_raiseFlags( softfloat_flag_invalid );
            uiZ = defaultNaNF64UI;
            goto uiZ;
        } */
        sigDiff = sigA - sigB;

        // |A| - |B| == 0 case, skips softfloat_normRoundPackToF64 (1*)
        if ( ! sigDiff ) {
            // 0 is produced here
            uiZ =
                packToF64UI(
                    (softfloat_roundingMode == softfloat_round_min), 0, 0 );
            goto uiZ;
        }
        if ( expA ) --expA;
        if ( sigDiff < 0 ) {
            signZ = ! signZ;
            sigDiff = -sigDiff;
        }
        shiftDist = softfloat_countLeadingZeros64( sigDiff ) - 11;
        expZ = expA - shiftDist;
        // no subnormals
        /* if ( expZ < 0 ) {
            shiftDist = expA;
            expZ = 0;
        } */
        uiZ = packToF64UI( signZ, expZ, sigDiff<<shiftDist );
        goto uiZ;
    } else {
        sigA <<= 10;
        sigB <<= 10;
        if ( expDiff < 0 ) {
            signZ = ! signZ;
            // no infinities
            /* if ( expB == 0x7FF ) {
                if ( sigB ) goto propagateNaN;
                uiZ = packToF64UI( signZ, 0x7FF, 0 );
                goto uiZ;
            } */
            sigA += expA ? UINT64_C( 0x4000000000000000 ) : sigA;
            sigA = softfloat_shiftRightJam64( sigA, -expDiff );
            sigB |= UINT64_C( 0x4000000000000000 );
            expZ = expB;
            sigZ = sigB - sigA;
        } else {
            // no infinities
            /* if ( expA == 0x7FF ) {
                if ( sigA ) goto propagateNaN;
                uiZ = uiA;
                goto uiZ;
            } */
            sigB += expB ? UINT64_C( 0x4000000000000000 ) : sigB;
            sigB = softfloat_shiftRightJam64( sigB, expDiff );
            sigA |= UINT64_C( 0x4000000000000000 );
            expZ = expA;
            sigZ = sigA - sigB;
        }
        return softfloat_normRoundPackToF64( signZ, expZ - 1, sigZ );
    }
    /*------------------------------------------------------------------------
    *------------------------------------------------------------------------*/
 uiZ:
    uZ.ui = uiZ;
    return uZ.f;

}


static float64_t f64_add(float64_t a, float64_t b) {
    union ui64_f64 uA;
    uint_fast64_t uiA;
    bool signA;
    union ui64_f64 uB;
    uint_fast64_t uiB;
    bool signB;

    uA.f = a;
    uiA = uA.ui;
    signA = signF64UI( uiA );
    uB.f = b;
    uiB = uB.ui;
    signB = signF64UI( uiB );
    if ( signA == signB ) {
        return softfloat_addMagsF64( uiA, uiB, signA );
    } else {
        return softfloat_subMagsF64( uiA, uiB, signA );
    }
}

__attribute__((noinline))
v128_t soft_fadd_1(v128_t dest, v128_t src) {
    uint64_t dest_0 = wasm_i64x2_extract_lane(dest, 0);
    uint64_t src_0 = wasm_i64x2_extract_lane(src, 0);

    uint64_t dest_1 = wasm_i64x2_extract_lane(dest, 1);
    uint64_t src_1 = wasm_i64x2_extract_lane(src, 1);

    uint64_t res_0 = f64_add((float64_t){dest_0}, (float64_t){src_0}).v;
    uint64_t res_1 = f64_add((float64_t){dest_1}, (float64_t){src_1}).v;

    return wasm_i64x2_make(res_0, res_1);
}
