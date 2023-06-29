#include "vector2.h"
#include "polar.h"

namespace math
{
    const vector2 vector2::ZERO(0.0f, 0.0f);
    const vector2 vector2::UNIT_X(1.0f, 0.0f);
    const vector2 vector2::UNIT_Y(0.0f, 1.0f);
    const vector2 vector2::UNIT_SCALE(1.0f, 1.0f);

    vector2::operator polar () const
    {
        if(x == 0.0f && y == 0.0f) return polar(0.0f, 0.0f);
        else return polar(sqrt(x*x + y*y), atan2(y, x));
    }
}// namespace math
