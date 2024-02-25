#ifndef EULER_H_INCLUDED
#define EULER_H_INCLUDED

#include "angle.h"
#include "math_predefs.h"

namespace math
{
    /// \brief Euler angles (Tait-Bryan angles).
    /// First, rotate system around y axis (yaw),
    /// then rotate around derived x axis (pitch)
    /// and then rotate around derived z axis (roll)
    ///
    template<typename T>
    class euler
    {
    public:
        T yaw; ///< amount of rotation around y axis (radians)
        T pitch;///< amount of rotation around x axis (radians)
        T roll;///< amount of rotation around z axis (radians)

        /// \brief Default constructor.
        /// initializes yaw, pitch, roll with zero.
        ///
        constexpr euler() : yaw(0), pitch(0), roll(0)
        {
        }

        /// \brief Initializes yaw, pitch, roll with given values
        ///
        constexpr euler(T yaw, T pitch, T roll) : yaw(yaw), pitch(pitch), roll(roll)
        {
        }


        /// Canonize Euler angles.
        /// Canonical form:
        /// -PI < yaw <= PI,
        /// -PI/2 <= pitch <= PI/2,
        /// -PI < roll < PI
        ///
        /// TODO: tests show failure, to fix
        /*void canonize()
        {
            yaw -= floor(yaw / (pi<T>() * 2)) * pi<T>() * 2;
            if(yaw > pi<T>())
                yaw -= pi<T>() * 2;
            else if(yaw < -pi<T>())
                yaw += pi<T>() * 2;

            pitch -= floor(pitch / (pi<T>() * 2)) * pi<T>() * 2;
            if (pitch) > pi<T>() / 2)
            {
                pitch = pi<T>() - pitch;
                roll += pi<T>();
            }

            roll -= floor(roll / (pi<T>() * 2)) * pi<T>() * 2;
            if(roll > pi<T>())
                roll -= pi<T>() * 2;
            else if(yaw < -pi<T>())
                roll += pi<T>() * 2;
        }*/
    };
} // namespace math

#endif // EULER_H_INCLUDED
