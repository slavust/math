#ifndef VECTOR2_H_INCLUDED
#define VECTOR2_H_INCLUDED

#include "math_predefs.h"
#include <array>

namespace math
{
    class polar;

    /// \brief Vector in 2D Cartesian coordinate system
    ///
    class vector2
    {
    public:
        real x, y; ///< components

        static const vector2 ZERO; ///< zero vector
        static const vector2 UNIT_X; ///< unit vector along X axis
        static const vector2 UNIT_Y; ///< unit vector along Y axis
        static const vector2 UNIT_SCALE; ///< unit vector along Z axis


        /// \brief Default constructor
        ///
        vector2()
        {
        }


        /// \brief Initializes all components by scalar
        ///
        /// \param scalar - scalar
        ///
        vector2(real scalar) : x(scalar), y(scalar)
        {
        }


        /// \brief Initializes components x, y with corresponding array elements
        ///
        /// \param src - array of two elements
        ///
        vector2(const std::array<real, 2>& src) : x(src[0]), y(src[1])
        {
        }


        /// \brief Initializes components with x, y
        ///
        /// \param x - x component
        /// \param y - y component
        ///
        vector2(real x, real y) : x(x), y(y)
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src - source 2D vector
        ///
        vector2(const vector2& src) : x(src.x), y(src.y)
        {
        }


        /// \brief Conversion to 2D polar coordinate system
        ///
        operator polar () const;


        /// \brief Magnitude of vector
        ///
        /// \return magnitude
        ///
        /// Geometrycally magnitude of vector is it's length
        ///
        real magnitude() const
        {
            return sqrt(x*x + y*y);
        }


        /// \brief Distance to vector b
        ///
        /// \param b: 2D vector
        /// \return distance
        ///
        real distance(const vector2& b) const
        {
            return (b - *this).magnitude();
        }


        /// \brief Normalize vector
        ///
        /// \return vector's magnitude before normalization
        ///
        /// Result is vector with the same direction and unit magnitude
        ///
        real normalize()
        {
            real m = magnitude();

            if(*this != ZERO)
            {
                *this /= m;
            }
            return magnitude();
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
        real dot(const vector2& v) const
        {
            return x*v.x + y*v.y;
        }

        vector2& operator = (const vector2& v)
        {
            x = v.x;
            y = v.y;
            return *this;
        }

        bool operator == (const vector2& v) const
        {
            return x == v.x && y == v.y;
        }

        bool operator != (const vector2& v) const
        {
            return x != v.x || y != v.y;
        }

        vector2 operator * (real scalar) const
        {
            return vector2(x * scalar, y * scalar);
        }

        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 2D vector
        ///
        vector2 operator * (const vector2& v) const
        {
            return vector2(x*v.x, y*v.y);
        }

        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 2D vector
        ///
        vector2& operator *= (const vector2& v)
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
        vector2 operator / (const vector2& v) const
        {
            return vector2(x/v.x, y/v.y);
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        vector2& operator /= (const vector2& v)
        {
            x /= v.x;
            y /= v.y;

            return *this;
        }

        vector2 operator / (real scalar) const
        {
            return vector2(x / scalar, y / scalar);
        }

        vector2& operator *= (real scalar)
        {
            x *= scalar;
            y *= scalar;
            return *this;
        }

        vector2& operator /= (real scalar)
        {
            x /= scalar;
            y /= scalar;
            return *this;
        }

        const vector2& operator + () const
        {
            return *this;
        }

        vector2 operator - () const
        {
            return vector2(-x, -y);
        }

        vector2 operator + (const vector2& v) const
        {
            return vector2(x + v.x, y + v.y);
        }

        vector2 operator - (const vector2& v) const
        {
            return vector2(x - v.x, y - v.y);
        }

        vector2& operator += (const vector2& v)
        {
            x += v.x;
            y += v.y;
            return *this;
        }

        vector2& operator -= (const vector2& v)
        {
            x -= v.x;
            y -= v.y;
            return *this;
        }
    };

    inline vector2 operator * (real scalar, const vector2& v)
    {
        return v*scalar;
    }

} // namespace math

#endif // VECTOR2_H_INCLUDED
