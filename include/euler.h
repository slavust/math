#ifndef EULER_H_INCLUDED
#define EULER_H_INCLUDED

#include "angle.h"
#include "math_predefs.h"
#include "matrix3x3.h"
#include "matrix4x4.h"
#include "quaternion.h"

namespace math
{
    /// \brief Euler angles (Tait-Bryan angles).
    /// First, rotate system around y axis (yaw),
    /// then rotate around derived x axis (pitch)
    /// and then rotate around derived z axis (roll)
    ///
    class euler
    {
    public:
        real yaw; ///< amount of rotation around y axis (radians)
        real pitch;///< amount of rotation around x axis (radians)
        real roll;///< amount of rotation around z axis (radians)

        /// \brief Default constructor.
        /// initializes yaw, pitch, roll with zero.
        ///
        euler() : yaw(0.0f), pitch(0.0f), roll(0.0f)
        {
        }

        /// \brief Initializes yaw, pitch, roll with given values
        ///
        euler(real yaw, real pitch, real roll) : yaw(yaw), pitch(pitch), roll(roll)
        {
        }

        /// \brief Convert to 3x3 rotation matrix
        ///
        operator matrix3x3 () const
        {
            matrix3x3 ret;
            toRotationMatrix(ret);
            return ret;
        }

        /// \brief Convert to 4x4 rotation matrix
        ///
        operator matrix4x4 () const
        {
            matrix4x4 ret;
            toRotationMatrix(ret);
            return ret;
        }

        /// \brief Convert to rotation quaternion
        ///
        operator quaternion () const
        {
            return toQuaternion();
        }

        /// Canonize Euler angles.
        /// Canonical form:
        /// -PI < yaw <= PI,
        /// -PI/2 <= pitch <= PI/2,
        /// -PI < roll < PI,
        /// fabs(pitch) == PI/2 --> yaw += -sign(pitch)*roll, roll = 0.
        ///
        void canonize();


        /// \brief Convert to 3x3 rotation matrix (by link)
        ///
        /// \param rot matrix3x3&: output matrix
        ///
        void toRotationMatrix(matrix3x3& rot) const;


        /// \brief Convert to 4x4 rotation matrix (by link)
        ///
        /// \param rot matrix4x4&: output matrix
        ///
        void toRotationMatrix(matrix4x4& rot) const;


        /// \brief Convert to rotation quaternion
        ///
        /// \param rot: output quaternion
        ///
        quaternion toQuaternion() const;
    };
} // namespace math

#endif // EULER_H_INCLUDED
