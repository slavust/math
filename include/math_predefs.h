#ifndef MATH_PREDEFS_H_INCLUDED
#define MATH_PREDEFS_H_INCLUDED

#include <math.h>
#include <string.h>
#include <limits>
#include "assert.h"

namespace math
{
    template <typename T>
    constexpr T sign(T s)
    {
        return s < 0 ? -1 : 1;
    }

    template <typename T>
    T eps()
    {
        return sqrt(std::numeric_limits<T>::epsilon());
    }
} // namespace math

#endif // MATH_PREDEFS_H_INCLUDED
