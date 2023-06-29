#include <math.h>

#include "matrix3x3.h"
#include "vector3.h"
#include "euler.h"

namespace math
{
    const matrix3x3 matrix3x3::IDENTITY(1.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f);

    matrix3x3 matrix3x3::operator * (real scalar) const
    {
        return matrix3x3(val[0][0] * scalar,
                         val[1][0] * scalar,
                         val[2][0] * scalar,

                         val[0][1] * scalar,
                         val[1][1] * scalar,
                         val[2][1] * scalar,

                         val[0][2] * scalar,
                         val[1][2] * scalar,
                         val[2][2] * scalar);
    }

    matrix3x3& matrix3x3::operator *= (real scalar)
    {
        val[0][0] *= scalar;
        val[1][0] *= scalar;
        val[2][0] *= scalar;

        val[0][1] *= scalar;
        val[1][1] *= scalar;
        val[2][1] *= scalar;

        val[0][2] *= scalar;
        val[1][2] *= scalar;
        val[2][2] *= scalar;

        return *this;
    }

    matrix3x3 matrix3x3::operator / (real scalar) const
    {
        return matrix3x3(val[0][0] / scalar,
                         val[1][0] / scalar,
                         val[2][0] / scalar,

                         val[0][1] / scalar,
                         val[1][1] / scalar,
                         val[2][1] / scalar,

                         val[0][2] / scalar,
                         val[1][2] / scalar,
                         val[2][2] / scalar);
    }

    matrix3x3& matrix3x3::operator /= (real scalar)
    {
        val[0][0] /= scalar;
        val[1][0] /= scalar;
        val[2][0] /= scalar;

        val[0][1] /= scalar;
        val[1][1] /= scalar;
        val[2][1] /= scalar;

        val[0][2] /= scalar;
        val[1][2] /= scalar;
        val[2][2] /= scalar;

        return *this;
    }

    matrix3x3 matrix3x3::operator * (const matrix3x3& b) const
    {

        return matrix3x3(val[0][0]*b[0][0] + val[1][0]*b[0][1] + val[2][0]*b[0][2],
                          val[0][1]*b[0][0] + val[1][1]*b[0][1] + val[2][1]*b[0][2],
                          val[0][2]*b[0][0] + val[1][2]*b[0][1] + val[2][2]*b[0][2],

                          val[0][0]*b[1][0] + val[1][0]*b[1][1] + val[2][0]*b[1][2],
                          val[0][1]*b[1][0] + val[1][1]*b[1][1] + val[2][1]*b[1][2],
                          val[0][2]*b[1][0] + val[1][2]*b[1][1] + val[2][2]*b[1][2],

                          val[0][0]*b[2][0] + val[1][0]*b[2][1] + val[2][0]*b[2][2],
                          val[0][1]*b[2][0] + val[1][1]*b[2][1] + val[2][1]*b[2][2],
                          val[0][2]*b[2][0] + val[1][2]*b[2][1] + val[2][2]*b[2][2]
                          );
    }

    matrix3x3 matrix3x3::adjoint() const
    {
        return matrix3x3(val[1][1]*val[2][2] - val[2][1]*val[1][2],
                         val[2][1]*val[0][2] - val[0][1]*val[2][2],
                         val[0][1]*val[1][2] - val[1][1]*val[0][2],

                         val[2][0]*val[1][2] - val[1][0]*val[2][2],
                         val[0][0]*val[2][2] - val[2][0]*val[0][2],
                         val[1][0]*val[0][2] - val[0][0]*val[1][2],

                         val[1][0]*val[2][1] - val[2][0]*val[1][1],
                         val[2][0]*val[0][1] - val[0][0]*val[2][1],
                         val[0][0]*val[1][1] - val[1][0]*val[0][1]);
    }

    void matrix3x3::orthonormalize()
    {
        vector3 p(val[0][0], val[0][1], val[0][2]);
        vector3 q(val[1][0], val[1][1], val[1][2]);
        vector3 r(val[2][0], val[2][1], val[2][2]);

        p.normalize();

        q = q - p.dot(q)*p;
        q.normalize();

        r = r - p.dot(r)*p - q.dot(r)*q;
        r.normalize();

        val[0][0] = p.x;
        val[0][1] = p.y;
        val[0][2] = p.z;

        val[1][0] = q.x;
        val[1][1] = q.y;
        val[1][2] = q.z;

        val[2][0] = r.x;
        val[2][1] = r.y;
        val[2][2] = r.z;
    }

    euler matrix3x3::toEuler() const
    {
        euler ret;

        real sin_p = -val[2][1];
        ret.pitch = asin(sin_p);
        real cos_p = cos(ret.pitch);

        if(fabs(cos_p) < EPS) // cos(ret.pitch) == 0
        {
            ret.pitch = PI_TWO;
            ret.roll = 0.0f;

            // sin_r = 0.0f;
            // cos_r = 1.0f;
            // tan(ret.yaw) = -val[0][3] / val[0][0];
            ret.yaw = atan2(-val[0][2], val[0][0]);
        }
        else// cos(ret.pitch) != 0
        {

            //real sin_y = val[2][0] / cos_p;
            //real cos_y = val[2][2] / cos_p;

            //assert(cos_y != 0.0f || sin_y != 0.0f);
            //ret.yaw = atan2(sin_y, cos_y);
            ret.yaw = atan2(val[2][0], val[2][2]);

            // cos(ret.pitch) != 0
            //real sin_r = val[0][1] / cos_p;
            //real cos_r = val[1][1] / cos_p;

            //assert(cos_r != 0.0f || sin_r != 0.0f);
            //ret.roll = atan2(sin_r, cos_r);
            ret.roll = atan2(val[0][1], val[1][1]);
        }

        return ret;
    }

    quaternion matrix3x3::toQuaternion() const
    {
        // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
        // article "Quaternion Calculus and Fast Animation".

        quaternion ret;

        real fTrace = val[0][0]+val[1][1]+val[2][2];
        real fRoot;

        if ( fTrace > 0 ) // > 0
        {
            // |w| > 1/2, may as well choose w > 1/2
            fRoot = sqrt(fTrace + 1.0f);  // 2w
            ret.w = 0.5f*fRoot;
            fRoot = 0.5f/fRoot;  // 1/(4w)
            ret.x = (val[1][2]-val[2][1])*fRoot;
            ret.y = (val[2][0]-val[0][2])*fRoot;
            ret.z = (val[0][1]-val[1][0])*fRoot;
        }
        else
        {
            // |w| <= 1/2
            static size_t s_iNext[3] = { 1, 2, 0 };
            size_t i = 0;
            if ( val[1][1] > val[0][0] )
                i = 1;
            if ( val[2][2] > val[i][i] )
                i = 2;
            size_t j = s_iNext[i];
            size_t k = s_iNext[j];

            fRoot = sqrt(val[i][i]-val[j][j]-val[k][k] + 1.0f);
            real* apkQuat[3] = { &ret.x, &ret.y, &ret.z };
            *apkQuat[i] = 0.5f*fRoot;
            fRoot = 0.5f/fRoot;
            ret.w = (val[j][k]-val[k][j])*fRoot;
            *apkQuat[j] = (val[i][j]+val[j][i])*fRoot;
            *apkQuat[k] = (val[i][k]+val[k][i])*fRoot;
        }
        return ret;
    }

    matrix3x3::operator euler() const
    {
        return toEuler();
    }

    matrix3x3::operator quaternion() const
    {
        return toQuaternion();
    }



} // namespace math
