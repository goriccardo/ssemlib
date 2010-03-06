/* Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>  *
 *                                                          *
 * This file is part of ssemlib and it's released under the *
 * BSD license, see LICENSE                                 */

#include "ssemath.h"

/* Mask definitions */
PI64_CONST(mantissa_mask, 0x7ff0000000000000LL);
PI64_CONST(inv_mantissa_mask, ~0x7ff0000000000000LL);
PI32_CONST(bias4i, 0x3fe);
PD_CONST(half_mask, 0.5);


/** Return the significand and the exponent both in double precision **/
__m128d frexp_sse(__m128d x, __m128d *e)
{
/* Integer exponent */
__m128i ei;

/* Save the exponent */
ei = _mm_castpd_si128(x);
ei = _mm_and_si128(ei, *(__m128i*)pi64_mantissa_mask);
ei = _mm_srli_epi64(ei, 52);
ei = _mm_shuffle_epi32(ei,216);
ei = _mm_sub_epi32(ei, *(__m128i*)pi32_bias4i);

*e = _mm_cvtepi32_pd(ei);

/* Save the significand */
x = _mm_and_pd(x, *(__m128d*)pi64_inv_mantissa_mask);
x = _mm_or_pd(x, *(__m128d*)pd_half_mask);

return x;
}
