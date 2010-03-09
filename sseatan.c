/* Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>  *
 *                                                          *
 * This file is part of ssemlib and it's released under the *
 * BSD license, see LICENSE                                 */
 
#include "ssemlib.h"

PD_CONST(one, 1.);

/* Mask */
PI64_CONST(sign_mask, 0x8000000000000000LL);
PI64_CONST(no_sign_mask, ~0x8000000000000000LL);
 
/* CEPHES constants */
PD_CONST(cephes_atan_p0, -8.750608600031904122785E-1);
PD_CONST(cephes_atan_p1, -1.615753718733365076637E1);
PD_CONST(cephes_atan_p2, -7.500855792314704667340E1);
PD_CONST(cephes_atan_p3, -1.228866684490136173410E2);
PD_CONST(cephes_atan_p4, -6.485021904942025371773E1);
PD_CONST(cephes_atan_q1, 2.485846490142306297962E1);
PD_CONST(cephes_atan_q2, 1.650270098316988542046E2);
PD_CONST(cephes_atan_q3, 4.328810604912902668951E2);
PD_CONST(cephes_atan_q4, 4.853903996359136964868E2);
PD_CONST(cephes_atan_q5, 1.945506571482613964425E2);

/* Trigonometric constants */
PD_CONST(T3P8, 2.41421356237309504880);     /* tan( 3*pi/8 ) */
PD_CONST(PI, 3.14159265358979323846);       /* pi */
PD_CONST(PIO2, 1.57079632679489661923);     /* pi/2 */
PD_CONST(PIO4, 7.85398163397448309616E-1);  /* pi/4 */
PD_CONST(Z66, 0.66);

__m128d atan_sse(__m128d x)
{
/* TODO, handle infinity */

/* Save and remove the sign */
__m128d smask = _mm_and_pd(x, *(__m128d*)pi64_sign_mask);
x = _mm_and_pd(x, *(__m128d*)pi64_no_sign_mask);

__m128d one = *(__m128d*)pd_one;

/* Range reduction:
    if( x > T3P8 ) {
     y = PIO2;
     flag = 1;
     x = -( 1.0/x );
    } else if( x <= 0.66 ) {
     y = 0.0;
    } else {
     y = PIO4;
     flag = 2;
     x = (x-1.0)/(x+1.0);
    } */
__m128d range_1 = _mm_cmple_pd(x, *(__m128d*)pd_T3P8);
__m128d y = _mm_andnot_pd(range_1, *(__m128d*)pd_PIO2);
__m128d num = _mm_and_pd(range_1, x);
__m128d den = _mm_and_pd(range_1, one);

__m128d range_2 = _mm_cmpgt_pd(x, *(__m128d*)pd_Z66);

range_1 = _mm_and_pd(range_2, range_1);

__m128d tmp = _mm_and_pd(range_1, *(__m128d*)pd_PIO4);
y = _mm_add_pd(y, tmp);

tmp = _mm_and_pd(range_2, one);
num = _mm_sub_pd(num, tmp);
tmp = _mm_and_pd(range_2, x);
den = _mm_add_pd(den, tmp);

x = _mm_div_pd(num, den);

/* Numeratore P(z) */
__m128d z = _mm_mul_pd(x, x);
__m128d p = _mm_mul_pd(z, *(__m128d*)pd_cephes_atan_p0);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_atan_p1);
p = _mm_mul_pd(p, z);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_atan_p2);
p = _mm_mul_pd(p, z);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_atan_p3);
p = _mm_mul_pd(p, z);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_atan_p4);

/* Denominatore Q(z) */
__m128d q = _mm_add_pd(z, *(__m128d*)pd_cephes_atan_q1);
q = _mm_mul_pd(q, z);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_atan_q2);
q = _mm_mul_pd(q, z);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_atan_q3);
q = _mm_mul_pd(q, z);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_atan_q4);
q = _mm_mul_pd(q, z);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_atan_q5);

/* Result */
p = _mm_div_pd(p, q);
z = _mm_mul_pd(p, z);
z = _mm_mul_pd(z, x);
z = _mm_add_pd(z, x);
y = _mm_add_pd(y, z);

/* TODO: Investigate about morebits */

/* Restore sign */
y = _mm_or_pd(y, smask);

return y;
}

__m128d atan2_sse(__m128d x, __m128d y)
{
__m128d mask = _mm_cmplt_pd(x,_mm_setzero_pd());  // x < 0
__m128d w = _mm_and_pd(mask, *(__m128d*)pd_PI);
mask = _mm_cmplt_pd(y, _mm_setzero_pd());
mask = _mm_and_pd(mask, *(__m128d*)pi64_sign_mask);
w = _mm_xor_pd(w, mask);

__m128d z = _mm_div_pd(x, y);
z = atan_sse(z);
z = _mm_add_pd(w, z);
return z;
}

