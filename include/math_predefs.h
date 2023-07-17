#ifndef MATH_PREDEFS_H_INCLUDED
#define MATH_PREDEFS_H_INCLUDED

#include <math.h>
#include <string.h>
#include <limits>
#include "assert.h"

namespace math
{
#ifdef SMALLMATH_DOUBLE_PRECISION
    typedef double real;
#else
    typedef float real;
#endif

    static const real EPS = std::sqrt(std::numeric_limits<real>::epsilon());

    inline real sign(real s)
    {
        return s < 0 ? -1.0f : 1.0f;
    }
} // namespace math

#endif // MATH_PREDEFS_H_INCLUDED
