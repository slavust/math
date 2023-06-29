#ifndef VECTOR3_H_INCLUDED
#define VECTOR3_H_INCLUDED

#include "math_predefs.h"

namespace math
{

    class cylindrical;
    class spherical;

    /// \brief Vector in 3D Cartesian
    /// or 2D homogeneous coordinate system
    ///
    class vector3
    {
    public:
        real x, y, z; ///< components

        static const vector3 ZERO; ///< zero vector
        static const vector3 UNIT_X; ///< unit vector along X axis
        static const vector3 UNIT_Y; ///< unit vector along Y axis
        static const vector3 UNIT_Z; ///< unit vector along Z axis
        static const vector3 UNIT_SCALE;  ///< vector with all components unit


        /// \brief Default constructor
        ///
        /// Initializes components with zero
        ///
        vector3() : x(0.0f), y(0.0f), z(0.0f)
        {
        }


        /// \brief Initializes all components by scalar
        ///
        /// \param scalar - scalar
        ///
        vector3(real scalar) : x(scalar), y(scalar), z(scalar)
        {
        }


        /// \brief Initializes components x, y, z with corresponding array elements
        ///
        /// \param src - array of three elements
        ///
        vector3(const real src[3]) : x(src[0]), y(src[1]), z(src[2])
        {
        }


        /// \brief Initializes components with x, y, z
        ///
        /// \param x - x component
        /// \param y - y component
        /// \param z - z component
        ///
        vector3(real x, real y, real z) : x(x), y(y), z(z)
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src - source 3D vector
        ///
        vector3(const vector3& src) : x(src.x), y(src.y), z(src.z)
        {
        }


        /// \brief Convertion to cylindrical coordinate space
        ///
        operator cylindrical () const;

        /// \brief Conversial to spherical coordinate space
        ///
        operator spherical () const;


        /// \brief Magnitude of vector
        ///
        /// \return magnitude
        ///
        /// Geometrycally magnitude of vector is it's length
        ///
        real magnitude() const
        {
            return sqrt(x*x + y*y + z*z);
        }


        /// \brief Distance to vector b
        ///
        /// \param b: 3D vector
        /// \return distance
        ///
        real distance(const vector3& b) const
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
            return m;
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
        real dot(const vector3& v) const
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
        vector3 cross(const vector3& v) const
        {
            return vector3(y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);
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

        bool operator == (const vector3& v) const
        {
            return x == v.x && y == v.y && z == v.z;
        }

        bool operator != (const vector3& v) const
        {
            return x != v.x || y != v.y || z != v.z;
        }


        vector3& operator = (const vector3& v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
            return *this;
        }

        vector3 operator * (real scalar) const
        {
            return vector3(x * scalar, y * scalar, z * scalar);
        }

        vector3 operator / (real scalar) const
        {
            return vector3(x / scalar, y / scalar, z / scalar);
        }

        vector3& operator *= (real scalar)
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        vector3& operator /= (real scalar)
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
        vector3 operator * (const vector3& v) const
        {
            return vector3(x*v.x, y*v.y, z*v.z);
        }

        /// \brief Component-wise multiplication
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        vector3& operator *= (const vector3& v)
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
        vector3 operator / (const vector3& v) const
        {
            return vector3(x/v.x, y/v.y, z/v.z);
        }

        /// \brief Component-wise division
        ///
        /// \param v: second operand
        /// \return 3D vector
        ///
        vector3& operator /= (const vector3& v)
        {
            x /= v.x;
            y /= v.y;
            z /= v.z;

            return *this;
        }

        const vector3& operator + () const
        {
            return *this;
        }

        vector3 operator - () const
        {
            return vector3(-x, -y, -z);
        }

        vector3 operator + (const vector3& v) const
        {
            return vector3(x + v.x, y + v.y, z + v.z);
        }

        vector3 operator - (const vector3& v) const
        {
            return vector3(x - v.x, y - v.y, z - v.z);
        }

        vector3& operator += (const vector3& v)
        {
            x += v.x;
            y += v.y;
            z += v.z;
            return *this;
        }

        vector3& operator -= (const vector3& v)
        {
            x -= v.x;
            y -= v.y;
            z -= v.z;
            return *this;
        }
    };

    inline vector3 operator * (real scalar, const vector3& v)
    {
        return v*scalar;
    }

} // namespace math

#endif // VECTOR3_H_INCLUDED
