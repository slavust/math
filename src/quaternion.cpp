#include "quaternion.h"
#include "euler.h"

namespace math
{
    const quaternion quaternion::ZERO(0, 0, 0, 0);
    const quaternion quaternion::IDENTITY(1, 0, 0, 0);


    void  quaternion::fromExponentialMap(const vector3& exp_map)
    {
        vector3 axis = exp_map;
        real angle = axis.normalize();

        fromAxisAngle(axis, angle);
    }

    void  quaternion::fromAxisAngle(const vector3& axis, real angle)
    {
        w = cos(angle / 2.0f);
        real s = sin(angle / 2.0f);

        x = s * axis.x;
        y = s * axis.y;
        z = s * axis.z;
    }

    void  quaternion::toAxisAngle(vector3& axis, real& angle) const
    {
        real length2 = x*x + y*y + z*z;

        if(length2 > 0.0f)
        {
            angle = 2.0f * acos(w);
            real inv_length = 1.0f / sqrt(length2);

            axis.x = x*inv_length;
            axis.y = y*inv_length;
            axis.z = z*inv_length;
        }
        else
        {
            angle = 0.0f;
            axis = vector3::UNIT_X;
        }
    }

    void  quaternion::toRotationMatrix(matrix3x3& rot) const
    {
        rot[0][0] = 1.0f - 2.0f*y*y - 2.0f*z*z;
        rot[0][1] = 2.0f*x*y + 2.0f*w*z;
        rot[0][2] = 2.0f*x*z - 2.0f*w*y;

        rot[1][0] = 2.0f*x*y - 2.0f*w*z;
        rot[1][1] = 1.0f - 2.0f*x*x - 2.0f*z*z;
        rot[1][2] = 2.0f*y*z + 2.0f*w*x;

        rot[2][0] = 2.0f*x*z + 2.0f*w*y;
        rot[2][1] = 2.0f*y*z - 2.0f*w*x;
        rot[2][2] = 1.0f - 2.0f*x*x - 2.0f*y*y;
    }

    void  quaternion::toRotationMatrix(matrix4x4& rot) const
    {
        rot[0][0] = 1.0f - 2.0f*y*y - 2.0f*z*z;
        rot[0][1] = 2.0f*x*y + 2.0f*w*z;
        rot[0][2] = 2.0f*x*z - 2.0f*w*y;
        rot[0][3] = 0.0f;

        rot[1][0] = 2.0f*x*y - 2.0f*w*z;
        rot[1][1] = 1.0f - 2.0f*x*x - 2.0f*z*z;
        rot[1][2] = 2.0f*y*z + 2.0f*w*x;
        rot[1][3] = 0.0f;

        rot[2][0] = 2.0f*x*z + 2.0f*w*y;
        rot[2][1] = 2.0f*y*z - 2.0f*w*x;
        rot[2][2] = 1.0f - 2.0f*x*x - 2.0f*y*y;
        rot[2][3] = 0.0f;

        rot[3][0] = 0.0f;
        rot[3][1] = 0.0f;
        rot[3][2] = 0.0f;
        rot[3][3] = 1.0f;
    }

    euler  quaternion::toEuler () const
    {
        euler ret;

        ret.pitch = asin(-2*(y*z - w*x));

        if(fabs(fabs(ret.pitch) - PI_TWO) > EPS) // cos(pitch) != 0
        {
            ret.yaw = atan2(x*z + w*y, 0.5f - x*x - y*y);
            ret.roll = atan2(x*y + w*z, 0.5f - x*x - z*z);
        }
        else // cos(pitch) == 0
        {
            ret.yaw = atan2(-x*z + w*y, 0.5f - y*y - z*z);
            ret.roll = 0.0f;
        }

        return ret;
    }

    quaternion quaternion::log() const
    {
        if(fabs(w) < 1.0f)
        {
            real a = acos(w);
            real s = sin(a);
            if(fabs(s) >= EPS)
            {
                real coeff = a / s;
                return quaternion(0, x*coeff, y*coeff, z*coeff);
            }
        }

        return quaternion(0, x, y, z);
    }

    quaternion quaternion::exp() const
    {
        real a = sqrt(x*x + y*y + z*z);
        real s = sin(a);

        if(fabs(s) >= EPS)
        {
            real coeff = s / a;
            return quaternion(cos(a), x*coeff, y*coeff, z*coeff);
        }

        return quaternion(cos(a), x, y, z);
    }

    quaternion quaternion::slerp(quaternion q, real t) const
    {
        quaternion ret;

        real d = dot(q);

        //choose signs of q1 and q2 such as q1.dot(q2) >= 0;
        if(d < 0)
        {
            q = -q;
            d = -d;
        }

        real a = acos(dot(q));

        real s = sin(d);

        if(s >= EPS)
        {
            ret = *this * sin(a*(1.0f - t)) / s + q * sin(a*t) / s;
        }
        else
        {
            quaternion q0 = *this;
            q0.inverse();
            quaternion e = ((q*q0).log()*t).exp();
            ret = e * (*this);
        }

        return ret;
    }

    quaternion::operator euler () const
    {
        return toEuler();
    }
} // namespace math
