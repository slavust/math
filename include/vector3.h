#ifndef VECTOR3_H_INCLUDED
#define VECTOR3_H_INCLUDED

#include "math_predefs.h"
#include <array>

namespace math
{
    /// \brief Vector in 3D Cartesian
    /// or 2D homogeneous coordinate system
    ///
    template<typename T>
    class vector3
    {
    public:
        T x, y, z; ///< components

        static constexpr vector3<T> ZERO {0, 0, 0}; ///< zero vector
        static constexpr vector3<T> UNIT_X {1, 0, 0}; ///< unit vector along X axis
        static constexpr vector3<T> UNIT_Y {0, 1, 0}; ///< unit vector along Y axis
        static constexpr vector3<T> UNIT_Z {0, 0, 1}; ///< unit vector along Z axis
        static constexpr vector3<T> UNIT_SCALE {1, 1, 1};  ///< vector with all components unit


        /// \brief Default constructor
        ///
        /// Initializes components with zero
        ///
        constexpr vector3() : x(0), y(0), z(0) {}


        /// \brief Initializes components x, y, z with corresponding array elements
        ///
        /// \param src - array of three elements
        ///
        vector3(const std::array<T, 3>& src) : x(src[0]), y(src[1]), z(src[2])
        {
        }


        /// \brief Initializes components with x, y, z
        ///
        /// \param x - x component
        /// \param y - y component
        /// \param z - z component
        ///
        constexpr vector3(T x, T y, T z) : x(x), y(y), z(z) {}


        /// \brief Copy constructor
        ///
        /// \param src - source 3D vector
        ///
        constexpr vector3(const vector3<T>& src) : x(src.x), y(src.y), z(src.z) {}


        /// \brief Magnitude of vector
        ///
        /// \return magnitude
        ///
        /// Geometrycally magnitude of vector is it's length
        ///
        constexpr T magnitude() const
        {
            return sqrt(x*x + y*y + z*z);
        }


        /// \brief Distance to vector b
        ///
        /// \param b: 3D vector
        /// \return distance
        ///
        constexpr T distance(const vector3<T>& b) const
        {
            return (b - *this).magnitude();
        }


        /// \brief Normalize vector
        ///
        /// \return vector's magnitude before normalization
        ///
        /// Result is vector with the same direction and unit magnitude
        ///
        T normalize()
        {
            const T m = magnitude();

            if(*this != ZERO)
            {
                *this /= m;
            }
            return m;
        }

        constexpr vector3<T> normalized() const
        {
            const T m = magnitude();

            if(*this != ZERO)
            {
                return *this / m;
            }
            else
            {
                return ZERO;
            }
        }


        /// \brief Dot product of two vectors
        ///
        /// \param v: second vector
        /// \return dot product
        ///
        /// Dot product of two 3D vectors a, b is (a.x*b.x + a.y*b.y + a.z*b.z).
        /// Geometrically, it is the product of the magnitudes of
        /// a, b and the cosine of the angle between them.
        ///
        constexpr T dot(const vector3<T>& v) const
        {
            return x*v.x + y*v.y + z*v.z;
        }


        /// \brief Cross product
        ///
        /// \param v: 3D vector
        /// \return 3D vector
        ///
        /// The result of cross product of two 3D vectors a, b
        /// is perpependicular to both a and b.
        ///
        constexpr vector3<T> cross(const vector3<T>& v) const
        {
            return vector3<T>(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
        }

        /// conversion operator to vector2 in case if vector3 is used for translation in 2d
        constexpr operator vector2<T>() const
        {
            return {x / z, y / z};
        }

        constexpr bool operator == (const vector3<T>& v) const
        {
            return x == v.x && y == v.y && z == v.z;
        }

        constexpr bool operator != (const vector3<T>& v) const
        {
            return x != v.x || y != v.y || z != v.z;
        }


        vector3<T>& operator = (const vector3<T>& v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
            return *this;
        }

        constexpr vector3<T> operator * (T scalar) const
        {
            return vector3<T>(x * scalar, y * scalar, z * scalar);
        }

        constexpr vector3<T> operator / (T scalar) const
        {
            return vector3<T>(x / scalar, y / scalar, z / scalar);
        }

        vector3<T>& operator *= (T scalar)
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        vector3<T>& operator /= (T scalar)
        {
            x /= scalar;
            y /= scalar;
            z /= scalar;
            return *this;
        }

        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        constexpr vector3<T> operator * (const vector3<T>& v) const
        {
            return vector3<T>(x*v.x, y*v.y, z*v.z);
        }

        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        vector3<T>& operator *= (const vector3<T>& v)
        {
            x *= v.x;
            y *= v.y;
            z *= v.z;

            return *this;
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        constexpr vector3<T> operator / (const vector3<T>& v) const
        {
            return vector3<T>(x/v.x, y/v.y, z/v.z);
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        vector3<T>& operator /= (const vector3<T>& v)
        {
            x /= v.x;
            y /= v.y;
            z /= v.z;

            return *this;
        }

        constexpr vector3<T> operator + () const
        {
            return *this;
        }

        constexpr vector3<T> operator - () const
        {
            return vector3<T>(-x, -y, -z);
        }

        constexpr vector3<T> operator + (const vector3<T>& v) const
        {
            return vector3<T>(x + v.x, y + v.y, z + v.z);
        }

        constexpr vector3 operator - (const vector3<T>& v) const
        {
            return vector3<T>(x - v.x, y - v.y, z - v.z);
        }

        vector3<T>& operator += (const vector3<T>& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        vector3<T>& operator -= (const vector3<T>& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
    };

    template<typename T>
    constexpr vector3<T> operator * (T scalar, const vector3<T>& v)
    {
        return v*scalar;
    }
} // namespace math

#endif // VECTOR3_H_INCLUDED
