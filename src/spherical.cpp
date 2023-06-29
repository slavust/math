#include "spherical.h"
#include "vector3.h"

namespace math
{
    spherical::operator vector3 () const
    {
        real _sin_theta = sin(theta);
        return vector3(r*cos(phi)*_sin_theta, r*sin(phi)*_sin_theta, r*cos(theta));
    }

    void spherical::canonize()
    {
        if(r == 0.0f) // at the origin
        {
            // at r == 0, theta = phi = 0
            theta = 0.0f;
            phi = 0.0f;
        }
        else
        {
            if(r < 0.0f) // negative radius
            {
                // make positive radius
                r = -r;
                phi += PI;
            }

            // make theta within range [0, PI]
            if(theta < 0.0f)
            {
                // make positive theta
                theta = -theta;
                phi += PI;
            }
            if(theta > PI)
            {
                // subtract from TWO_PI theta without TWO_PI portions. 0<=theta<=PI
                theta = TWO_PI - theta + floor(theta / TWO_PI) * TWO_PI;
                phi += PI;
            }

            // at theta = 0 and theta = PI, phi = 0
            if(theta < EPS || fabs(theta - PI) < EPS) phi = 0.0f;
            else if(fabs(phi) > PI)
            {
                //make phi within range [-PI, PI)

                // offset by PI
                // maps from -PI*n > phi > PI*n to -TWO_PI*(n-1)/2 > phi > TWO_PI*(n+1)/2
                //phi += PI;

                // subtract signed TWO_PI portions
                //phi -= floor(phi / TWO_PI) * TWO_PI; // 0<=phi<TWO_PI
                phi -= floor((phi + PI) / TWO_PI) * TWO_PI;

                // undo offset
                // maps from 0<=phi<TWO_PI to -PI<=phi<PI
                //phi -= PI;
            }
        }
    }

} // namespace math
