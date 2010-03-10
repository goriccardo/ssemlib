/* Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>  *
 *                                                          *
 * This file is part of ssemlib and it's released under the *
 * BSD license, see LICENSE                                 */

#include "ssemlib.h"

PD_CONST(one, 1.);
PD_CONST(half, 0.5);
PI64_CONST(min_norm, 0x0010000000000000);

/* CEPHES constants */
PD_CONST(cephes_SQRTHF, 0.70710678118654752440);
PD_CONST(cephes_log_p0, 1.01875663804580931796E-4);
PD_CONST(cephes_log_p1, 4.97494994976747001425E-1);
PD_CONST(cephes_log_p2, 4.70579119878881725854E0);
PD_CONST(cephes_log_p3, 1.44989225341610930846E1);
PD_CONST(cephes_log_p4, 1.79368678507819816313E1);
PD_CONST(cephes_log_p5, 7.70838733755885391666E0);
PD_CONST(cephes_log_q1, 1.12873587189167450590E1);
PD_CONST(cephes_log_q2, 4.52279145837532221105E1);
PD_CONST(cephes_log_q3, 8.29875266912776603211E1);
PD_CONST(cephes_log_q4, 7.11544750618563894466E1);
PD_CONST(cephes_log_q5, 2.31251620126765340583E1);
PD_CONST(cephes_log_ec1,2.121944400546905827679e-4);
PD_CONST(cephes_log_ec2,0.693359375);


/** Return the natural logarithm of x **/
__m128d log_sse(__m128d x)
{
/* Exponent */
__m128d e;

/* Check for positiveness */
__m128d invalid_mask = _mm_cmple_pd(x, _mm_setzero_pd());

/* Avoid denormalized numbers */
x = _mm_max_pd(x, *(__m128d*)pi64_min_norm);

x = frexp_sse(x, &e);

/*
  if( x < SQRTHF ) {
     e -= 1;
     x = x + x - 1.0;
   } else { x = x - 1.0; }
*/
__m128d mask = _mm_cmplt_pd(x, *(__m128d*)pd_cephes_SQRTHF);
__m128d q = _mm_and_pd(x, mask);
x = _mm_sub_pd(x, *(__m128d*)pd_one);
e = _mm_sub_pd(e, _mm_and_pd(*(__m128d*)pd_one, mask));
x = _mm_add_pd(x, q);

/* x^2 */
__m128d z = _mm_mul_pd(x,x);

/* Numeratore P(x) */
__m128d p = _mm_mul_pd(x, *(__m128d*)pd_cephes_log_p0);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_log_p1);
p = _mm_mul_pd(p, x);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_log_p2);
p = _mm_mul_pd(p, x);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_log_p3);
p = _mm_mul_pd(p, x);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_log_p4);
p = _mm_mul_pd(p, x);
p = _mm_add_pd(p, *(__m128d*)pd_cephes_log_p5);

/* Denominatore Q(x) */
q = _mm_add_pd(x, *(__m128d*)pd_cephes_log_q1);
q = _mm_mul_pd(q, x);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_log_q2);
q = _mm_mul_pd(q, x);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_log_q3);
q = _mm_mul_pd(q, x);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_log_q4);
q = _mm_mul_pd(q, x);
q = _mm_add_pd(q, *(__m128d*)pd_cephes_log_q5);

/* P(x)/Q(x) */
p = _mm_div_pd(p, q);

/* x^3 * P(x)/Q(x) */
p = _mm_mul_pd(p, z);
p = _mm_mul_pd(p, x);

/* y = y - e*ec - z/2 */
q = _mm_mul_pd(e, *(__m128d*)pd_cephes_log_ec1);
z = _mm_mul_pd(z, *(__m128d*)pd_half);
p = _mm_sub_pd(p, q);
p = _mm_sub_pd(p, z);

q = _mm_mul_pd(e, *(__m128d*)pd_cephes_log_ec2);
z = _mm_add_pd(x,p);
z = _mm_add_pd(z,q);

z = _mm_or_pd(z, invalid_mask); // negative arg will be NAN

return z;
}
