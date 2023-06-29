#include "matrix4x4.h"
#include "euler.h"
#include "quaternion.h"

namespace math
{

    const matrix4x4 matrix4x4::IDENTITY(1.0f, 0.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, 0.0f, 1.0f);

    matrix4x4 matrix4x4::operator * (real scalar) const
    {
        return matrix4x4(val[0][0] * scalar,
                         val[0][1] * scalar,
                         val[0][2] * scalar,
                         val[0][3] * scalar,

                         val[1][0] * scalar,
                         val[1][1] * scalar,
                         val[1][2] * scalar,
                         val[1][3] * scalar,

                         val[2][0] * scalar,
                         val[2][1] * scalar,
                         val[2][2] * scalar,
                         val[2][3] * scalar,

                         val[3][0] * scalar,
                         val[3][1] * scalar,
                         val[3][2] * scalar,
                         val[3][3] * scalar);
    }

    matrix4x4& matrix4x4::operator *= (real scalar)
    {
        val[0][0] *= scalar;
        val[0][1] *= scalar;
        val[0][2] *= scalar;
        val[0][3] *= scalar;

        val[1][0] *= scalar;
        val[1][1] *= scalar;
        val[1][2] *= scalar;
        val[1][3] *= scalar;

        val[2][0] *= scalar;
        val[2][1] *= scalar;
        val[2][2] *= scalar;
        val[2][3] *= scalar;

        val[3][0] *= scalar;
        val[3][1] *= scalar;
        val[3][2] *= scalar;
        val[3][3] *= scalar;

        return *this;
    }

    matrix4x4 matrix4x4::operator / (real scalar) const
    {
        return matrix4x4(val[0][0] / scalar,
                         val[0][1] / scalar,
                         val[0][2] / scalar,
                         val[0][3] / scalar,

                         val[1][0] / scalar,
                         val[1][1] / scalar,
                         val[1][2] / scalar,
                         val[1][3] / scalar,

                         val[2][0] / scalar,
                         val[2][1] / scalar,
                         val[2][2] / scalar,
                         val[2][3] / scalar,

                         val[3][0] / scalar,
                         val[3][1] / scalar,
                         val[3][2] / scalar,
                         val[3][3] / scalar);
    }

    matrix4x4& matrix4x4::operator /= (real scalar)
    {
        val[0][0] /= scalar;
        val[0][1] /= scalar;
        val[0][2] /= scalar;
        val[0][3] /= scalar;

        val[1][0] /= scalar;
        val[1][1] /= scalar;
        val[1][2] /= scalar;
        val[1][3] /= scalar;

        val[2][0] /= scalar;
        val[2][1] /= scalar;
        val[2][2] /= scalar;
        val[2][3] /= scalar;

        val[3][0] /= scalar;
        val[3][1] /= scalar;
        val[3][2] /= scalar;
        val[3][3] /= scalar;

        return *this;
    }

    matrix4x4 matrix4x4::operator * (const matrix4x4& b) const
    {

        return matrix4x4(
                          val[0][0]*b[0][0] + val[1][0]*b[0][1] + val[2][0]*b[0][2] + val[3][0]*b[0][3],
                          val[0][1]*b[0][0] + val[1][1]*b[0][1] + val[2][1]*b[0][2] + val[3][1]*b[0][3],
                          val[0][2]*b[0][0] + val[1][2]*b[0][1] + val[2][2]*b[0][2] + val[3][2]*b[0][3],
                          val[0][3]*b[0][0] + val[1][3]*b[0][1] + val[2][3]*b[0][2] + val[3][3]*b[0][3],

                          val[0][0]*b[1][0] + val[1][0]*b[1][1] + val[2][0]*b[1][2] + val[3][0]*b[1][3],
                          val[0][1]*b[1][0] + val[1][1]*b[1][1] + val[2][1]*b[1][2] + val[3][1]*b[1][3],
                          val[0][2]*b[1][0] + val[1][2]*b[1][1] + val[2][2]*b[1][2] + val[3][2]*b[1][3],
                          val[0][3]*b[1][0] + val[1][3]*b[1][1] + val[2][3]*b[1][2] + val[3][3]*b[1][3],

                          val[0][0]*b[2][0] + val[1][0]*b[2][1] + val[2][0]*b[2][2] + val[3][0]*b[2][3],
                          val[0][1]*b[2][0] + val[1][1]*b[2][1] + val[2][1]*b[2][2] + val[3][1]*b[2][3],
                          val[0][2]*b[2][0] + val[1][2]*b[2][1] + val[2][2]*b[2][2] + val[3][2]*b[2][3],
                          val[0][3]*b[2][0] + val[1][3]*b[2][1] + val[2][3]*b[2][2] + val[3][3]*b[2][3],

                          val[0][0]*b[3][0] + val[1][0]*b[3][1] + val[2][0]*b[3][2] + val[3][0]*b[3][3],
                          val[0][1]*b[3][0] + val[1][1]*b[3][1] + val[2][1]*b[3][2] + val[3][1]*b[3][3],
                          val[0][2]*b[3][0] + val[1][2]*b[3][1] + val[2][2]*b[3][2] + val[3][2]*b[3][3],
                          val[0][3]*b[3][0] + val[1][3]*b[3][1] + val[2][3]*b[3][2] + val[3][3]*b[3][3]);
    }


    real matrix4x4::determinant() const
    {
        return val[0][0] * (val[1][1]*(val[2][2]*val[3][3] - val[3][2]*val[2][3])
                            + val[2][1]*(val[3][2]*val[1][3] - val[1][2]*val[3][3])
                            + val[3][1]*(val[1][2]*val[2][3] - val[2][2]*val[1][3]))

             - val[1][0] * (val[0][1]*(val[2][2]*val[3][3] - val[3][2]*val[2][3])
                            + val[2][1]*(val[3][2]*val[0][3] - val[0][2]*val[3][3])
                            + val[3][1]*(val[0][2]*val[2][3] - val[2][2]*val[0][3]))

             + val[2][0] * (val[0][1]*(val[1][2]*val[3][3] - val[3][2]*val[1][3])
                            + val[1][1]*(val[3][2]*val[0][3] - val[0][2]*val[3][3])
                            + val[3][1]*(val[0][2]*val[1][3] - val[1][2]*val[0][3]))

             - val[3][0] * (val[0][1]*(val[1][2]*val[2][3] - val[2][2]*val[1][3])
                            + val[1][1]*(val[2][2]*val[0][3] - val[0][2]*val[2][3])
                            + val[2][1]*(val[0][2]*val[1][3] - val[1][2]*val[0][3]));
    }

    matrix4x4 matrix4x4::adjoint() const
    {
        return matrix4x4(
             val[1][1]*val[2][2]*val[3][3] + val[2][1]*val[3][2]*val[1][3] + val[3][1]*val[1][2]*val[2][3] - val[3][1]*val[2][2]*val[1][3] - val[2][1]*val[1][2]*val[3][3] - val[1][1]*val[3][2]*val[2][3],
            -val[0][1]*val[2][2]*val[3][3] - val[2][1]*val[3][2]*val[0][3] - val[3][1]*val[0][2]*val[2][3] + val[3][1]*val[2][2]*val[0][3] + val[2][1]*val[0][2]*val[3][3] + val[0][1]*val[3][2]*val[2][3],
             val[0][1]*val[1][2]*val[3][3] + val[1][1]*val[3][2]*val[0][3] + val[3][1]*val[0][2]*val[1][3] - val[3][1]*val[1][2]*val[0][3] - val[1][1]*val[0][2]*val[3][3] - val[0][1]*val[3][2]*val[1][3],
            -val[0][1]*val[1][2]*val[2][3] - val[1][1]*val[2][2]*val[0][3] - val[2][1]*val[0][2]*val[1][3] + val[2][1]*val[1][2]*val[0][3] + val[1][1]*val[0][2]*val[2][3] + val[0][1]*val[2][2]*val[1][3],

            -val[1][0]*val[2][2]*val[3][3] - val[2][0]*val[3][2]*val[1][3] - val[3][0]*val[1][2]*val[2][3] + val[3][0]*val[2][2]*val[1][3] + val[2][0]*val[1][2]*val[3][3] + val[1][0]*val[3][2]*val[2][3],
             val[0][0]*val[2][2]*val[3][3] + val[2][0]*val[3][2]*val[0][3] + val[3][0]*val[0][2]*val[2][3] - val[3][0]*val[2][2]*val[0][3] - val[2][0]*val[0][2]*val[3][3] - val[0][0]*val[3][2]*val[2][3],
            -val[0][0]*val[1][2]*val[3][3] - val[1][0]*val[3][2]*val[0][3] - val[3][0]*val[0][2]*val[1][3] + val[3][0]*val[1][2]*val[0][3] + val[1][0]*val[0][2]*val[3][3] + val[0][0]*val[3][2]*val[1][3],
             val[0][0]*val[1][2]*val[2][3] + val[1][0]*val[2][2]*val[0][3] + val[2][0]*val[0][2]*val[1][3] - val[2][0]*val[1][2]*val[0][3] - val[1][0]*val[0][2]*val[2][3] - val[0][0]*val[2][2]*val[1][3],

             val[1][0]*val[2][1]*val[3][3] + val[2][0]*val[3][1]*val[1][3] + val[3][0]*val[1][1]*val[2][3] - val[3][0]*val[2][1]*val[1][3] - val[2][0]*val[1][1]*val[3][3] - val[1][0]*val[3][1]*val[2][3],
            -val[0][0]*val[2][1]*val[3][3] - val[2][0]*val[3][1]*val[0][3] - val[3][0]*val[0][1]*val[2][3] + val[3][0]*val[2][1]*val[0][3] + val[2][0]*val[0][1]*val[3][3] + val[0][0]*val[3][1]*val[2][3],
             val[0][0]*val[1][1]*val[3][3] + val[1][0]*val[3][1]*val[0][3] + val[3][0]*val[0][1]*val[1][3] - val[3][0]*val[1][1]*val[0][3] - val[1][0]*val[0][1]*val[3][3] - val[0][0]*val[3][1]*val[1][3],
            -val[0][0]*val[1][1]*val[2][3] - val[1][0]*val[2][1]*val[0][3] - val[2][0]*val[0][1]*val[1][3] + val[2][0]*val[1][1]*val[0][3] + val[1][0]*val[0][1]*val[2][3] + val[0][0]*val[2][1]*val[1][3],

            -val[1][0]*val[2][1]*val[3][2] - val[2][0]*val[3][1]*val[1][2] - val[3][0]*val[1][1]*val[2][2] + val[3][0]*val[2][1]*val[1][2] + val[2][0]*val[1][1]*val[3][2] + val[1][0]*val[3][1]*val[2][2],
             val[0][0]*val[2][1]*val[3][2] + val[2][0]*val[3][1]*val[0][2] + val[3][0]*val[0][1]*val[2][2] - val[3][0]*val[2][1]*val[0][2] - val[2][0]*val[0][1]*val[3][2] - val[0][0]*val[3][1]*val[2][2],
            -val[0][0]*val[1][1]*val[3][2] - val[1][0]*val[3][1]*val[0][2] - val[3][0]*val[0][1]*val[1][2] + val[3][0]*val[1][1]*val[0][2] + val[1][0]*val[0][1]*val[3][2] + val[0][0]*val[3][1]*val[1][2],
             val[0][0]*val[1][1]*val[2][2] + val[1][0]*val[2][1]*val[0][2] + val[2][0]*val[0][1]*val[1][2] - val[2][0]*val[1][1]*val[0][2] - val[1][0]*val[0][1]*val[2][2] - val[0][0]*val[2][1]*val[1][2]);

    }

    matrix4x4::operator euler() const
    {
        return toEuler();
    }

    matrix4x4::operator quaternion() const
    {
        return toQuaternion();
    }

    euler matrix4x4::toEuler() const
    {
        euler ret;

        real sin_p = -val[2][1];
        ret.pitch = asin(sin_p);
        real cos_p = cos(ret.pitch);

        if(fabs(cos_p) < EPS) // cos(ret.pitch) == 0
        {
            ret.pitch = PI_TWO;
            ret.roll = 0.0f;
            ret.yaw = atan2(-val[0][2], val[0][0]);
        }
        else// cos(ret.pitch) != 0
        {
            ret.yaw = atan2(val[2][0], val[2][2]);
            ret.roll = atan2(val[0][1], val[1][1]);
        }

        return ret;
    }

    quaternion matrix4x4::toQuaternion() const
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
} // namespace math
