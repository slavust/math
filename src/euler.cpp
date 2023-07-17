#include "euler.h"
#include "matrix3x3.h"
#include "matrix4x4.h"
#include "matrix_transforms.h"

namespace math {
    void euler::canonize()
    {
        if (fabs(pitch) > PI)
        {
            //subtract TWO_PI portions
            pitch -= floor((pitch + PI) / TWO_PI) * TWO_PI;
        }
        if (fabs(pitch) > PI_TWO)
        {
            pitch = PI - pitch;
            roll += PI;
        }
        else if (fabs(fabs(pitch) - PI / 2) < EPS)
        {
            yaw += -sign(pitch) * roll;
            roll = 0.0f;
        }
        if (fabs(roll) > PI)
        {
            roll -= floor((roll + PI) / TWO_PI) * TWO_PI;
        }

        if (fabs(yaw) > PI)
        {
            yaw -= floor((yaw + PI) / TWO_PI) * TWO_PI;
        }
    }

    void euler::toRotationMatrix(matrix3x3 &rot) const
    {
        rot = rotate3Dy(yaw) * rotate3Dx(pitch) * rotate3Dz(roll);
    }

    void euler::toRotationMatrix(matrix4x4 &rot) const
    {
        rot = rotate3Dhy(yaw) * rotate3Dhx(pitch) * rotate3Dhz(roll);
    }

    quaternion euler::toQuaternion() const
    {
        return quaternion(vector3::UNIT_Y, yaw) * quaternion(vector3::UNIT_X, pitch)
               * quaternion(vector3::UNIT_Z, roll);
    }

} // namespace math
