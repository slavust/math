#ifndef VECTOR4_H_INCLUDED
#define VECTOR4_H_INCLUDED

#include "math_predefs.h"
#include "vector3.h"

namespace math
{
    /// \brief Vector in 3D homogeneous coordinate system
    ///
    class vector4
    {
    public:
        real x, y, z, w; ///< components

        static const vector4 ZERO; ///< zero vector
        static const vector4 UNIT_X; ///< unit vector along X axis
        static const vector4 UNIT_Y; ///< unit vector along Y axis
        static const vector4 UNIT_Z; ///< unit vector along Z axis
        static const vector4 UNIT_W; ///< unit vector along W axis
        static const vector4 UNIT_SCALE;  ///< vector with all components unit


        /// \brief Default constructor
        ///
        /// Initializes components with zero
        ///
        vector4() : x(0), y(0), z(0), w(0)
        {
        }

        /// \brief Initializes all components by scalar
        ///
        /// \param scalar - scalar
        ///
        vector4(real scalar) : x(scalar), y(scalar), z(scalar), w(scalar)
        {
        }


        /// \brief Initializes components with x, y, z
        ///
        /// \param x - x component
        /// \param y - y component
        /// \param z - z component
        /// \param w - w component
        ///
        vector4(real x, real y, real z, real w) : x(x), y(y), z(z), w(w)
        {
        }


        /// \brief Initializes components x, y, z, w with corresponding array elements
        ///
        /// \param src - array of three elements
        ///
        vector4(real src[4]) : x(src[0]), y(src[1]), z(src[2]), w(src[3])
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src - source 3D homogeneous vector
        ///
        vector4(const vector4& src) : x(src.x), y(src.y), z(src.z), w(src.w)
        {
        }


        /// \brief Initializes components x, y, z, w with corresponding src elements and w = 1.0f
        ///
        /// \param src - source 3D vector
        ///
        vector4(const vector3& src) : x(src.x), y(src.y), z(src.z), w(1)
        {
        }


        /// \brief Pointer to array of elements
        ///
        /// \return real*
        ///
        /// Useful for copying with memcpy etc.
        ///
        real* ptr()
        {
            return &x;
        }


        /// \brief Pointer to constant array of elements
        ///
        /// \return const real*
        ///
        /// Useful for copying with memcpy etc.
        ///
        const real* ptr() const
        {
            return &x;
        }


        real& operator [] (size_t indx)
        {
            return ptr()[indx];
        }

        real operator [] (size_t indx) const
        {
            return ptr()[indx];
        }

        bool operator == (const vector4& v) const
        {
            return x == v.x && y == v.y && z == v.z && w == v.w;
        }

        bool operator != (const vector4& v) const
        {
            return x != v.x || y != v.y || z != v.z || w != v.w;
        }

        vector4& operator = (const vector4& v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
            w = v.w;
            return *this;
        }

        vector4 operator * (real scalar) const
        {
            return vector4(x * scalar, y * scalar, z * scalar, w * scalar);
        }

        vector4 operator / (real scalar) const
        {
            return vector4(x / scalar, y / scalar, z / scalar, w / scalar);
        }

        vector4& operator *= (real scalar)
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            w *= scalar;
            return *this;
        }

        vector4& operator /= (real scalar)
        {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            w /= scalar;
            return *this;
        }


        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 4D homogeneous vector
        ///
        vector4 operator * (const vector4& v) const
        {
            return vector4(x*v.x, y*v.y, z*v.z, w*v.w);
        }


        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 4D homogeneous vector
        ///
        vector4& operator *= (const vector4& v)
        {
            x *= v.x;
            y *= v.y;
            z *= v.z;
            w *= v.w;

            return *this;
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 4D homogeneous vector
        ///
        vector4 operator / (const vector4& v) const
        {
            return vector4(x/v.x, y/v.y, z/v.z, w/v.w);
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 4D homogeneous vector
        ///
        vector4& operator /= (const vector4& v)
        {
            x /= v.x;
            y /= v.y;
            z /= v.z;
            w /= v.w;

            return *this;
        }

        const vector4& operator + () const
        {
            return *this;
        }

        vector4 operator - () const
        {
            return vector4(-x, -y, -z, -w);
        }

        vector4 operator + (const vector4& v) const
        {
            return vector4(x + v.x, y + v.y, z + v.z, w + v.w);
        }

        vector4 operator - (const vector4& v) const
        {
            return vector4(x - v.x, y - v.y, z - v.z, w - v.w);
        }

        vector4& operator += (const vector4& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            w += v.w;
            return *this;
        }

        vector4& operator -= (const vector4& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            w -= v.w;
            return *this;
        }
    };

    inline vector4 operator * (real scalar, const vector4& v)
    {
        return v*scalar;
    }

} // namespace math

#endif // VECTOR4_H_INCLUDED
