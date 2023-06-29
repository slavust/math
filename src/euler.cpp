#include "euler.h"
#include "matrix3x3.h"
#include "matrix4x4.h"
#include "matrix_transforms.h"

namespace math
{
    void euler::canonize()
    {
        if(fabs(pitch) > PI)
        {
            //subtract TWO_PI portions
            pitch -= floor((pitch + PI) / TWO_PI) * TWO_PI;
        }
        if(fabs(pitch) > PI_TWO)
        {
            pitch = PI - pitch;
            roll += PI;
        }
        else if(fabs(fabs(pitch) - PI/2) < EPS)
        {
            yaw += -sign(pitch)*roll;
            roll = 0.0f;
        }
        if(fabs(roll) > PI)
        {
            roll -= floor((roll + PI) / TWO_PI) * TWO_PI;
        }

        if(fabs(yaw) > PI)
        {
            yaw -= floor((yaw + PI) / TWO_PI) * TWO_PI;
        }
    }

    void euler::toRotationMatrix(matrix3x3& rot) const
    {
        real c1 = cos(yaw);
        real c2 = cos(pitch);
        real c3 = cos(roll);
        real s1 = sin(yaw);
        real s2 = sin(pitch);
        real s3 = sin(roll);

        rot[0][0] = c1*c3 + s1*s2*s3;
        rot[1][0] = c3*s1*s2 - c1*s1;
        rot[2][0] = c2*s1;

        rot[0][1] = c2*s3;
        rot[1][1] = c2*c3;
        rot[2][1] = -s2;

        rot[0][2] = c1*s2*s3 - c3*s1;
        rot[1][2] = s1*s3 + c1*c3*s2;
        rot[2][2] = c1*c2;
    }

    void euler::toRotationMatrix(matrix4x4& rot) const
    {
        rot = rotate3Dhy(yaw) * rotate3Dhx(pitch) * rotate3Dhz(roll);
/*
        real c1 = cos(yaw);
        real c2 = cos(pitch);
        real c3 = cos(roll);
        real s1 = sin(yaw);
        real s2 = sin(pitch);
        real s3 = sin(roll);

        rot[0][0] = c1*c3 + s1*s2*s3;
        rot[1][0] = c3*s1*s2 - c1*s1;
        rot[2][0] = c2*s1;
        rot[3][0] = 0.0f;

        rot[0][1] = c2*s3;
        rot[1][1] = c2*c3;
        rot[2][1] = -s2;
        rot[3][1] = 0.0f;

        rot[0][2] = c1*s2*s3 - c3*s1;
        rot[1][2] = s1*s3 + c1*c3*s2;
        rot[2][2] = c1*c2;
        rot[3][2] = 0.0f;

        rot[0][3] = rot[1][3] = rot[2][3] = 0.0f;
        rot[3][3] = 1.0f;

        //rot = matrix4x4(c1*c3 + s1*s2*s3, c3*s1*s2 - c1*s1, c2*s1, 0.0f,
        //                            c2*s3,            c2*c3,   -s2, 0.0f,
        //                 c1*s2*s3 - c3*s1, s1*s3 + c1*c3*s2, c1*c2, 0.0f,
        //                             0.0f,             0.0f,  0.0f, 1.0f);
*/
    }

    quaternion euler::toQuaternion() const
    {
        real chy = cos(yaw * 0.5f);
        real chp = cos(pitch * 0.5f);
        real chr = cos(roll * 0.5f);

        real shy = sin(yaw * 0.5f);
        real shp = sin(pitch * 0.5f);
        real shr = sin(roll * 0.5f);

        return quaternion(chy*chp*chr + shy*shp*shr,
                          chy*shp*chr + shy*chp*shr,
                          shy*chp*chr - chy*shp*shr,
                          chy*chp*shr - shy*shp*chr);

        // todo: compare with:
        //return quaternion(vector3::UNIT_Y, yaw) * quaternion(vector3::UNIT_X, pitch)
        //    * quaternion(vector3::UNIT_Z, roll);
    }

} // namespace math
