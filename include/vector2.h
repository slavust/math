#pragma once

#include "math_predefs.h"
#include <array>

namespace math
{
    template<typename T>
    class vector2
    {
    public:
        T x, y; ///< components

        static constexpr vector2<T> ZERO{0, 0}; ///< zero vector
        static constexpr vector2<T> UNIT_X{1, 0}; ///< unit vector along X axis
        static constexpr vector2<T> UNIT_Y{0, 1}; ///< unit vector along Y axis
        static constexpr vector2<T> UNIT_SCALE{1, 1}; ///< unit vector along Z axis


        /// \brief Default constructor
        ///
        constexpr vector2(): x(0), y(0)
        {
        }

        /// \brief Initializes components x, y with corresponding array elements
        ///
        /// \param src - array of two elements
        ///
        constexpr vector2(const std::array<T, 2>& src) : x(src[0]), y(src[1])
        {
        }


        /// \brief Initializes components with x, y
        ///
        /// \param x - x component
        /// \param y - y component
        ///
        constexpr vector2(T x, T y) : x(x), y(y)
        {
        }


        /// \brief Magnitude of vector
        ///
        /// \return magnitude
        ///
        /// Geometrically magnitude of vector is it's length
        ///
        constexpr T magnitude() const
        {
            return sqrt(x*x + y*y);
        }


        /// \brief Distance to vector b
        ///
        /// \param b: 2D vector
        /// \return distance
        ///
        constexpr T distance(const vector2<T>& b) const
        {
            return (b - *this).magnitude();
        }


        /// \brief Normalize vector
        ///
        /// \return vector's magnitude before normalization
        ///
        /// Result is vector with the same direction and unit magnitude
        ///
        constexpr T normalize()
        {
            const T m = magnitude();

            if(*this != ZERO)
            {
                *this /= m;
            }
            return magnitude();
        }

        /// \brief Normalized copy of vector
        ///
        /// \return vector normalized vector
        ///
        /// Result is vector with the same direction and unit magnitude
        ///
        constexpr vector2<T> normalized() const
        {
            const T m = magnitude();

            if(*this != ZERO)
                return *this / m;

            return *this;
        }


        /// \brief Dot product of two vectors
        ///
        /// \param v: second vector
        /// \return dot product
        ///
        /// Dot product of two 2D vectors a, b is (a.x*b.x + a.y*b.y).
        /// Geometrically, it is the product of the magnitudes of
        /// a, b and the cosine of the angle between them.
        ///
        constexpr T dot(const vector2<T>& v) const
        {
            return x*v.x + y*v.y;
        }

        constexpr vector2<T>& operator = (const vector2<T>& v)
        {
            x = v.x;
            y = v.y;
            return *this;
        }

        constexpr bool operator == (const vector2<T>& v) const
        {
            return x == v.x && y == v.y;
        }

        constexpr bool operator != (const vector2<T>& v) const
        {
            return x != v.x || y != v.y;
        }

        constexpr vector2<T> operator * (T scalar) const
        {
            return vector2<T>(x * scalar, y * scalar);
        }

        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 2D vector
        ///
        constexpr vector2<T> operator * (const vector2<T>& v) const
        {
            return vector2<T>(x*v.x, y*v.y);
        }
        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 2D vector
        ///
        vector2<T>& operator *= (const vector2<T>& v)
        {
            x *= v.x;
            y *= v.y;
            return *this;
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        constexpr vector2<T> operator / (const vector2<T>& v) const
        {
            return vector2<T>(x/v.x, y/v.y);
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        vector2<T>& operator /= (const vector2<T>& v)
        {
            x /= v.x;
            y /= v.y;

            return *this;
        }

        constexpr vector2<T> operator / (T scalar) const
        {
            return vector2<T>(x / scalar, y / scalar);
        }

        vector2<T>& operator *= (T scalar)
        {
            x *= scalar;
            y *= scalar;
            return *this;
        }

        vector2<T>& operator /= (T scalar)
        {
            x /= scalar;
            y /= scalar;
            return *this;
        }

        constexpr vector2<T> operator + () const
        {
            return *this;
        }

        constexpr vector2<T> operator - () const
        {
            return vector2<T>(-x, -y);
        }

        constexpr vector2<T> operator + (const vector2<T>& v) const
        {
            return vector2<T>(x + v.x, y + v.y);
        }

        constexpr vector2<T> operator - (const vector2<T>& v) const
        {
            return vector2<T>(x - v.x, y - v.y);
        }

        vector2<T>& operator += (const vector2<T>& v)
        {
            x += v.x;
            y += v.y;
            return *this;
        }

        vector2<T>& operator -= (const vector2<T>& v)
        {
            x -= v.x;
            y -= v.y;
            return *this;
        }
    };

    template<typename T>
    constexpr vector2<T> operator * (T scalar, const vector2<T>& v)
    {
        return v * scalar;
    }

} // namespace math