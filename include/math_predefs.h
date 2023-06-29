#ifndef MATH_PREDEFS_H_INCLUDED
#define MATH_PREDEFS_H_INCLUDED

#include <math.h>
#include <cstring>

namespace math
{
    using real = float;

    constexpr real EPS = 1e-04f; //std::numeric_limits<real>::epsilon();

    constexpr real sign(real s)
    {
        return s < 0 ? -1.0f : 1.0f;
    }

} // namespace math

#endif // MATH_PREDEFS_H_INCLUDED
