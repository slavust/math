#ifndef QUATERNION_H_INCLUDED
#define QUATERNION_H_INCLUDED

#include "math_predefs.h"
#include "vector3.h"
#include "matrix3x3.h"
#include "matrix4x4.h"

namespace math
{
    class euler;

    /// \brief Quaternion.
    /// Useful for rotations in 3D space.
    ///
    class quaternion
    {
        public:
        real w, x, y, z;

        static const quaternion ZERO;
        static const quaternion IDENTITY;

        /// \brief Default constructor.
        /// Initializes w, x, y, z with zero
        ///
        quaternion() : w(1.0f), x(0.0f), y(0.0f), z(0.0f)
        {
        }

        /// \brief initializes w, x, y, z with given values
        ///
        /// \param w real
        /// \param x real
        /// \param y real
        /// \param z real
        ///
        quaternion(real w, real x, real y, real z) : w(w), x(x), y(y), z(z)
        {
        }

        /// \brief Initializes rotation quaternion from axis-angle form
        ///
        /// \param axis: rotation axis
        /// \param angle: amount of rotation (radians)
        ///
        quaternion(const vector3& axis, real angle)
        {
            fromAxisAngle(axis, angle);
        }

        /// \brief Initializes rotation quaternion from exponential map form
        ///
        /// \param exp_map: exponential map
        ///
        /// See fromExponentialMap
        ///
        quaternion(const vector3& exp_map)
        {
            fromExponentialMap(exp_map);
        }

        /// \brief Convertion to Euler angles
        ///
        operator euler () const;


        /// \brief Convertion to 3x3 rotation matrix
        ///
        operator matrix3x3 () const
        {
            matrix3x3 ret;
            toRotationMatrix(ret);
            return ret;
        }


        /// \brief Convertion to 4x4 rotation matrix
        ///
        operator matrix4x4 () const
        {
            matrix4x4 ret;
            toRotationMatrix(ret);
            return ret;
        }


        /// \brief Convert from exponential map rotation representation to quaternion
        ///
        /// \param exp_map: 3D vector representing exponential map
        /// \return void
        ///
        /// Amount of rotation in radians is exp_map.length(),
        /// rotation axis is normalized exp_map.
        /// See fromAxisAngle().
        ///
        void fromExponentialMap(const vector3& exp_map);


        /// \brief Convert from axis angle rotation representation to quaternion
        ///
        /// \param angle: umount of rotation (radians)
        /// \param axis: unit rotation axis
        /// \return void
        ///
        /// Rotation by angle Q around unit axis N is
        /// quaternion(cos(Q/2), sin(Q/2)*N.x, sin(Q/2)*N.y, sin(Q/2)*N.z)
        ///
        void fromAxisAngle(const vector3& axis, real angle);


        /// \brief Convert from quaternion rotation representation to axis angle
        ///
        /// \param axis: rotation axis
        /// \param angle: amount of rotation
        /// \return void
        ///
        void toAxisAngle(vector3& axis, real& angle) const;


        /// \brief Convert from quaternion rotation representation to matrix
        ///
        /// \param rot: rotation matrix
        /// \return void
        ///
        ///
        void toRotationMatrix(matrix3x3& rot) const;


        /// \brief Convert from quaternion rotation representation to matrix (homogeneous space)
        ///
        /// \param rot: rotation matrix
        /// \return void
        ///
        ///
        void toRotationMatrix(matrix4x4& rot) const;


        /// \brief Convert from quaternion rotation representation to Euler angles
        ///
        /// \return euler
        ///
        euler toEuler () const;


        /// \brief Quaternion negation.
        /// Represents the same angular displacement: -q = q.
        ///
        quaternion operator - () const
        {
            return quaternion(-w, -x, -y, -z);
        }

        quaternion operator + () const
        {
            return quaternion(w, x, y, z);
        }

        /// \brief Magnitude of quaternion.
        /// Magnitude of rotation quaternion is 1.
        ///
        /// \return real
        ///
        real magnitude() const
        {
            return sqrt(w*w + x*x + y*y + z*z);
        }

        /// \brief Inverse of quaternion.
        /// Inverse of rotation quaterniojn is rotation in opposite direction
        ///
        /// \return quaternion
        ///
        quaternion inverse() const
        {
            return quaternion(w, -x, -y, -z) / magnitude();
        }

        quaternion operator * (real scalar) const
        {
            return quaternion(w*scalar, x*scalar, y*scalar, z*scalar);
        }

        quaternion& operator *= (real scalar)
        {
            w *= scalar;
            x *= scalar;
            y *= scalar;
            z *= scalar;

            return *this;
        }

        quaternion operator * (const quaternion& q) const
        {
            return quaternion(w*q.w - x*q.x - y*q.y - z*q.z,
                              q.y*z - q.z*y + q.w*x + w*q.x,
                              q.z*x - q.x*z + q.w*y + w*q.y,
                              q.x*y - q.y*x + q.w*z + w*q.z);
        }

        /// \brief Quaternion by vector multiplication.
        /// Rotates vector by quaternion: q * v * q^(-1)
        /// \param
        /// \param
        /// \return
        ///
        ///
        vector3 operator * (const vector3& v) const
        {
            return vector3(v.y*z - v.z*y + x + w*v.x,
                           v.z*x - v.x*z + y + w*v.y,
                           v.x*y - v.y*x + z + w*v.z);
        }

        friend vector3 operator * (const vector3& v, const quaternion& q)
        {
            return vector3(q.y*v.z - q.z*v.y + q.w*v.x + q.x,
                           q.z*v.x - q.x*v.z + q.w*v.y + q.y,
                           q.x*v.y - q.y*v.x + q.w*v.z + q.z);
        }

        quaternion& operator *= (const quaternion q)
        {
            quaternion tmp(w, x, y, z);

            w = tmp.w*q.w - tmp.x*q.x - tmp.y*q.y - tmp.z*q.z;
            x = tmp.w*q.x + tmp.x*q.w + tmp.y*q.z - tmp.z*q.y;
            y = tmp.w*q.y + tmp.y*q.w + tmp.z*q.x - tmp.x*q.z;
            z = tmp.w*q.z + tmp.z*q.w + tmp.x*q.y - tmp.y*q.x;

            return *this;
        }

        quaternion operator / (real scalar) const
        {
            return quaternion(w/scalar, x/scalar, y/scalar, z/scalar);
        }

        quaternion& operator /= (real scalar)
        {
            w /= scalar;
            x /= scalar;
            y /= scalar;
            z /= scalar;

            return *this;
        }

        quaternion operator + (const quaternion& q) const
        {
            return quaternion(w + q.w, x + q.x, y + q.y, z + q.z);
        }

        quaternion& operator += (const quaternion& q)
        {
            w += q.w;
            x += q.x;
            y += q.y;
            z += q.z;

            return *this;
        }

        quaternion operator - (const quaternion& q) const
        {
            return quaternion(w - q.w, x - q.x, y - q.y, z - q.z);
        }

        quaternion& operator -= (const quaternion& q)
        {
            w -= q.w;
            x -= q.x;
            y -= q.y;
            z -= q.z;

            return *this;
        }

        /// \brief Quaternion dot product
        /// Geometrically: the larger absolute value of rotation quaternions
        /// dot product the more similary their rotations are.
        ///
        /// \param q: second operand
        /// \return real
        ///
        real dot(const quaternion& q) const
        {
            return w*q.w + x*q.x + y*q.y + z*q.z;
        }

        /// \brief Logarythm of quaternion
        /// If a = Q/2, N - unit vector then
        /// logarithm of quaternion (cos(a), N*sin(a)) is quaternion (0, a*N)
        ///
        /// \return quaternion
        ///
        quaternion log() const;


        /// \brief Quaternion exponential.
        /// Exponential of quaternion (0, a*N) is quaternion (cos(a), N*sin(a))
        ///
        /// \return quaternion
        ///
        quaternion exp() const;


        /// \brief Quaternion exponetiation.
        /// As t varies from [0, 1], quaternion exponentiation varies
        /// from (1, 0) to source quaternion. Useful for extracting
        /// a fraction of rotation.
        /// \param t real
        /// \return quaternion
        ///
        ///
        quaternion exponentiation(real t)
        {
            return (log()*t).exp();
        }

        /// \brief Spherical linear interpolation between two quaternions.
        /// Allows to smoothly interpolate between two orientations.
        /// As t varies from 0 to 1, the slerp varies from source quaternion to q.
        ///
        /// \param q: second quaternion
        /// \param t: interpolation factor
        /// \return quaternion
        ///
        quaternion slerp(quaternion q, real t) const;
    };

} // namespace math

#endif // QUATERNION_H_INCLUDED
