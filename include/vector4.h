#ifndef VECTOR4_H_INCLUDED
#define VECTOR4_H_INCLUDED

#include "math_predefs.h"
#include "vector3.h"
#include <array>

namespace math
{
    /// \brief Vector in 3D homogeneous coordinate system
    ///
    template <typename T>
    class vector4
    {
    public:
        T x, y, z, w; ///< components

        static constexpr vector4<T> ZERO {0, 0, 0, 0}; ///< zero vector
        static constexpr vector4<T> UNIT_X {1, 0, 0, 1}; ///< unit vector along X axis
        static constexpr vector4<T> UNIT_Y {0, 1, 0, 1}; ///< unit vector along Y axis
        static constexpr vector4<T> UNIT_Z {0, 0, 1, 1}; ///< unit vector along Z axis
        static constexpr vector4<T> UNIT_W {0, 0, 0, 1}; ///< unit vector along W axis
        static constexpr vector4<T> UNIT_SCALE {1, 1, 1, 1};  ///< vector with all components unit


        /// \brief Default constructor
        ///
        /// Initializes components with zero (except unit w)
        ///
        constexpr vector4() : x(0), y(0), z(0), w(1)
        {
        }


        /// \brief Initializes components with x, y, z
        ///
        /// \param x - x component
        /// \param y - y component
        /// \param z - z component
        /// \param w - w component
        ///
        vector4(T x, T y, T z, T w) : x(x), y(y), z(z), w(w)
        {
        }


        /// \brief Initializes components x, y, z, w with corresponding array elements
        ///
        /// \param src - array of three elements
        ///
        constexpr vector4(const std::array<T, 4>& src) : x(src[0]), y(src[1]), z(src[2]), w(src[3])
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src - source 3D homogeneous vector
        ///
        constexpr vector4(const vector4<T>& src) : x(src.x), y(src.y), z(src.z), w(src.w)
        {
        }


        /// \brief Initializes components x, y, z, w with corresponding src elements and w = 1
        ///
        /// \param src - source 3D vector
        ///
        constexpr vector4(const vector3<T>& src) : x(src.x), y(src.y), z(src.z), w(1)
        {
        }

        /// convertion operator to vector3
        constexpr operator vector3<T>() const
        {
            return vector3<T>{x / w, y / w, z / w};
        }

        constexpr bool operator == (const vector4<T>& v) const
        {
            return x == v.x && y == v.y && z == v.z && w == v.w;
        }

        constexpr bool operator != (const vector4<T>& v) const
        {
            return x != v.x || y != v.y || z != v.z || w != v.w;
        }

        vector4& operator = (const vector4<T>& v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
            w = v.w;
            return *this;
        }

        constexpr vector4<T> operator * (T scalar) const
        {
            return vector4<T>(x * scalar, y * scalar, z * scalar, w * scalar);
        }

        constexpr vector4<T> operator / (T scalar) const
        {
            return vector4<T>(x / scalar, y / scalar, z / scalar, w / scalar);
        }

        vector4<T>& operator *= (T scalar)
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            w *= scalar;
            return *this;
        }

        vector4<T>& operator /= (T scalar)
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
        constexpr vector4<T> operator * (const vector4<T>& v) const
        {
            return vector4<T>(x*v.x, y*v.y, z*v.z, w*v.w);
        }


        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 4D homogeneous vector
        ///
        vector4<T>& operator *= (const vector4<T>& v)
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
        constexpr vector4 operator / (const vector4<T>& v) const
        {
            return vector4<T>(x/v.x, y/v.y, z/v.z, w/v.w);
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 4D homogeneous vector
        ///
        vector4<T>& operator /= (const vector4<T>& v)
        {
            x /= v.x;
            y /= v.y;
            z /= v.z;
            w /= v.w;

            return *this;
        }

        constexpr vector4<T> operator + () const
        {
            return this;
        }

        constexpr vector4<T> operator - () const
        {
            return vector4<T>(-x, -y, -z, -w);
        }

        constexpr vector4 operator + (const vector4& v) const
        {
            return vector4<T>(x + v.x, y + v.y, z + v.z, w + v.w);
        }

        constexpr vector4<T> operator - (const vector4& v) const
        {
            return vector4<T>(x - v.x, y - v.y, z - v.z, w - v.w);
        }

        vector4<T>& operator += (const vector4<T>& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            w += v.w;
            return *this;
        }

        vector4<T>& operator -= (const vector4<T>& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            w -= v.w;
            return *this;
        }
    };

    template <typename T>
    inline vector4<T> operator * (T scalar, const vector4<T>& v)
    {
        return v * scalar;
    }
} // namespace math

#endif // VECTOR4_H_INCLUDED
