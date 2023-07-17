#ifndef ANGLE_H_INCLUDED
#define ANGLE_H_INCLUDED

#include "math_predefs.h"
#include <cstring>

namespace math
{
    constexpr real PI = static_cast<real>(3.14159265358979323846);
    constexpr real TWO_PI = PI*2.0f;
    constexpr real PI_TWO = PI/2.0f;

    /// \brief Convert degrees to radians
    ///
    /// \param degrees real: value in degrees
    /// \return real: value in radians
    ///
    constexpr real toRadians(real degrees)
    {
        constexpr real _PI_180 = PI / 180.0f;
        return degrees * _PI_180;
    }

    /// \brief Convert radians to degrees
    ///
    /// \param radians real: value in radians
    /// \return real: value in degrees
    ///
    constexpr real toDegrees(real radians)
    {
        constexpr real _180_PI = 180.0f / PI;
        return radians * _180_PI;
    }
}

#endif // ANGLE_H_INCLUDED
