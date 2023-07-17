#ifndef MATH_PREDEFS_H_INCLUDED
#define MATH_PREDEFS_H_INCLUDED

#include <math.h>
#include <string.h>
#include <limits>
#include "assert.h"

namespace math
{
#ifdef SMALLMATH_DOUBLE_PRECISION
    using real = double;
#else
    using real = float;
#endif

    constexpr real EPS = std::sqrt(std::numeric_limits<real>::epsilon());

    constexpr real sign(real s)
    {
        return s < 0 ? -1.0f : 1.0f;
    }
} // namespace math

#endif // MATH_PREDEFS_H_INCLUDED
