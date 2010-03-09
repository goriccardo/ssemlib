/* Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>  *
 * Copyright (c) 2007  Julien Pommier                       *
 *                                                          *
 * This file is part of ssemlib and it's released under the *
 * BSD license, see LICENSE                                 */

#ifndef SSE_MATH_H
#define SSE_MATH_H

#include <emmintrin.h>


/* The *CONST macros, the frexp, and log functions are a    *
 * conversion of the excellent Julien Pommier's sse math    *
 * library: http://gruntthepeon.free.fr/ssemath/            */

#define PD_CONST(Name, Val)                                                  \
  static const double pd_##Name[2] __attribute__((aligned(16))) = { Val, Val }
#define PI64_CONST(Name, Val)                                                   \
  static const int64_t pi64_##Name[2] __attribute__((aligned(16))) = { Val, Val }
#define PI32_CONST(Name, Val)                                                         \
  static const int pi32_##Name[4] __attribute__((aligned(16))) = { Val, Val, Val, Val }

__m128d asinh_sse(__m128d x);
__m128d atan_sse(__m128d x);
// __m128d atan2_sse(__m128d x, __m128d y);
__m128d frexp_sse(__m128d x, __m128d *e);
__m128d log_sse(__m128d x);

#endif //SSE_MATH_H

