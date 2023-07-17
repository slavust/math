#include "polar.h"
#include "vector2.h"

namespace math
{
    polar::operator vector2 () const
    {
        return vector2((real)rho*cos(phi), (real)rho*sin(phi));
    }

    void polar::canonize()
    {
        if(rho == 0.0f) // at the origin
        {
            // at the origin phi = 0.0f
            phi = 0.0f;
        }
        else
        {
            if(rho < 0.0f) // negative distance
            {
                //make positive distance
                rho = -rho;
                phi += PI;
            }
            if(fabs(phi) > PI) // phi out of range
            {
                //make phi within range (-PI, PI]

                // offset by PI
                // maps phi from -PI*n < phi < PI*n to -TWO_PI*(n-1)/2 < phi < TWO_PI*(n+1)/2
                //phi += PI;

                // subtract signed TWO_PI portions
                //phi -= floor((phi) / TWO_PI) * TWO_PI; // 0<=phi<TWO_PI

                // undo offset
                // maps from 0<=phi<TWO_PI to -PI<=phi<PI
                //phi -= PI;

                phi -= floor((phi + PI) / TWO_PI) * TWO_PI;
            }
        }
    }
} // namespace math
