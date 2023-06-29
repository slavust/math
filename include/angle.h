#ifndef ANGLE_H_INCLUDED
#define ANGLE_H_INCLUDED

#include "math_predefs.h"

namespace math
{
    static const real PI = 3.14159265358979323846f;
    static const real TWO_PI = PI*2.0f;
    static const real PI_TWO = PI/2.0f;

    /// \brief Convert degrees to radians
    ///
    /// \param degrees real: value in degrees
    /// \return real: value in radians
    ///
    inline real toRadians(real degrees)
    {
        static const real _PI_180 = PI / 180.0f;
        return degrees * _PI_180;
    }

    /// \brief Convert radians to degrees
    ///
    /// \param radians real: value in radians
    /// \return real: value in degrees
    ///
    inline real toDegrees(real radians)
    {
        static const real _180_PI = 180.0f / PI;
        return radians * _180_PI;
    }
}

#endif // ANGLE_H_INCLUDED
