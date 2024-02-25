#ifndef ANGLE_H_INCLUDED
#define ANGLE_H_INCLUDED

#include "math_predefs.h"
#include <cstring>

namespace math
{
    template<typename T> constexpr T pi()
    {
        return static_cast<T>(M_PI); // todo: implement for higher precisions?
    }

    /// \brief Convert degrees to radians
    ///
    /// \param degrees real: value in degrees
    /// \return real: value in radians
    ///
    template <typename T>
    constexpr T to_radians(T degrees)
    {
        constexpr T _PI_180 = pi<T>() / 180;
        return degrees * _PI_180;
    }

    /// \brief Convert radians to degrees
    ///
    /// \param radians real: value in radians
    /// \return real: value in degrees
    ///
    template <typename T>
    constexpr T to_degrees(T radians)
    {
        constexpr T _180_PI = 180 / pi<T>();
        return radians * _180_PI;
    }
}

#endif // ANGLE_H_INCLUDED
