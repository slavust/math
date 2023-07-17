#include "vector3.h"
#include "cylindrical.h"
#include "spherical.h"

namespace math
{
    const vector3 vector3::ZERO(0.0f, 0.0f, 0.0f);
    const vector3 vector3::UNIT_X(1.0f, 0.0f, 0.0f);
    const vector3 vector3::UNIT_Y(0.0f, 1.0f, 0.0f);
    const vector3 vector3::UNIT_Z(0.0f, 0.0f, 1.0f);
    const vector3 vector3::UNIT_SCALE(1.0f, 1.0f, 1.0f);

    vector3::operator cylindrical () const
    {
        if(x == 0.0f && y == 0.0f) return cylindrical(0.0f, 0.0f, z);
        else return cylindrical(sqrt(x*x + y*y), atan2(y, x), z);
    }

    vector3::operator spherical () const
    {
        real r = magnitude();

        real theta, phi;

        theta = acos(z / r);

        if(x == 0.0f && y == 0.0f) phi = 0.0f;
        else phi = atan2(y, x);

        return spherical(r, theta, phi);
    }

} // namespace math
