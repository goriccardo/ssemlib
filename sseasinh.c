/* Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>  *
 *                                                          *
 * This file is part of ssemlib and it's released under the *
 * BSD license, see LICENSE                                 */

#include "ssemlib.h"

PD_CONST(one, 1.);

/* Mask */
PI64_CONST(sign_mask, 0x8000000000000000LL);
PI64_CONST(no_sign_mask, ~0x8000000000000000LL);

__m128d asinh_sse(__m128d x)
{
/* Save and remove the sign */
__m128d smask = _mm_and_pd(x, *(__m128d*)pi64_sign_mask);
x = _mm_and_pd(x, *(__m128d*)pi64_no_sign_mask);

/* Slower but more general */
__m128d z = _mm_mul_pd(x, x);
z = _mm_add_pd(z, *(__m128d*)pd_one);
z = _mm_sqrt_pd(z);
z = _mm_add_pd(z, x);
z = log_sse(z);

/* Restore sign */
z = _mm_or_pd(z, smask);

return z;
}

