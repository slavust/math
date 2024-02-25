#ifndef MATRIX_TRANSFORMS_H_INCLUDED
#define MATRIX_TRANSFORMS_H_INCLUDED

#include "math_predefs.h"
#include "vector2.h"
#include "vector3.h"
#include "vector4.h"
#include "matrix2x2.h"
#include "matrix3x3.h"
#include "matrix4x4.h"

namespace math
{
    // -------------- CONSTANTS --------------

    /// Reflection about X axis matrix (2D)
    template<typename T>
    constexpr matrix2x2<T> reflect_2Dx() {
        return matrix2x2<T>(-1, 0,
                            0, 1);
    }


    /// Reflection about Y axis matrix (2D)
    template<typename T>
    constexpr matrix2x2<T> reflect2Dy() {
        return matrix2x2<T>(1, 0,
                            0, -1);
    }


    /// Reflection about X axis matrix (2D homogeneous)
    template<typename T>
    constexpr matrix3x3<T> reflect2Dhx() {
        return matrix3x3<T>(-1, 0, 0,
                            0, 1, 0,
                            0, 0, 1);
    }


    /// Reflection about Y axis matrix (2D homogeneous)
    template<typename T>
    constexpr matrix3x3<T> reflect2Dhy() {
        return matrix3x3<T>(1, 0, 0,
                            0, -1, 0,
                            0, 0, 1);
    }


    /// Reflection about XY plane matrix (3D)
    template<typename T>
    constexpr matrix3x3<T> reflect3Dxy() {
        return matrix3x3<T>(1, 0, 0,
                            0, 1, 0,
                            0, 0, -1);
    }


    /// Reflection about XZ plane matrix (3D)
    template<typename T>
    constexpr matrix3x3<T> reflect3Dxz() {
        return matrix3x3<T>(1, 0, 0,
                            0, -1, 0,
                            0, 0, 1);
    }

    /// Reflection about YZ plane matrix (3D)
    template<typename T>
    constexpr matrix3x3<T> reflect3Dyz() {
        return matrix3x3<T>(-1, 0, 0,
                            0, 1, 0,
                            0, 0, 1);
    }


    /// Projection to XY plane matrix (3D)
    template<typename T>
    constexpr matrix3x3<T> project3Dxy() {
        return matrix3x3<T>(1, 0, 0,
                            0, 1, 0,
                            0, 0, 0);
    }


    /// Projection to XZ plane matrix (3D)
    template<typename T>
    constexpr matrix3x3<T> project3Dxz() {
        return matrix3x3<T>(1, 0, 0,
                            0, 0, 0,
                            0, 0, 1);
    }

    /// Projection to YZ plane matrix (3D)
    template<typename T>
    constexpr matrix3x3<T> project3Dyz() {
        return matrix3x3<T>(0, 0, 0,
                            0, 1, 0,
                            0, 0, 1);
    }


    // -------------- FUNCTIONS --------------

    /// \brief Construct 2x2 matrix for 2D rotation clockwise
    ///
    /// \param angle: amount of rotation in radians
    /// \return 2x2 rotation matrix
    ///
    template <typename T>
    matrix2x2<T> rotate2D(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix2x2(_cos, -_sin,
                          _sin, _cos);
    }

    /// \brief Construct 3x3 matrix for rotation clockwise in homogeneous 2D space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 2D rotation matrix
    ///
    template <typename T>
    matrix3x3<T> rotate2Dh(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix3x3<T>(_cos, -_sin, 0,
                         _sin, _cos, 0,
                         0, 0, 1);
    }


    /// \brief Construct 2x2 matrix for 2D scale
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \return 2x2 scale matrix
    ///
    template<typename T>
    constexpr matrix2x2<T> scale2D(T scale_x, T scale_y)
    {
        return matrix2x2(scale_x, 0,
                         0, scale_y);
    }


    /// \brief Construct 2x2 matrix for 2D scale
    ///
    /// \param scale: amount of scale along X and Y axes
    /// \return 2x2 scale matrix
    ///
    template <typename T>
    constexpr matrix2x2<T> scale2D(const vector2<T>& scale)
    {
        return scale2D<T>(scale.x, scale.y);
    }

    /// \brief Construct 3x3 matrix for scale in homogeneous 2D space
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \return 3x3 scale matrix
    ///
    template <typename T>
    constexpr matrix3x3<T> scale2Dh(T scale_x, T scale_y)
    {
        return matrix3x3<T>(scale_x, 0, 0,
                         0, scale_y, 0,
                         0,    0, 1);
    }


    /// \brief Construct 3x3 matrix for scale in homogeneous 2D space
    ///
    /// \param scale: amount of scale along X and Y axes
    /// \return 3x3 scale matrix
    ///
    template <typename T>
    constexpr matrix3x3<T> scale2Dh(const vector2<T>& scale)
    {
        return scale2Dh<T>(scale.x, scale.y);
    }


    /// \brief Construct 3x3 matrix for translation in homogeneous 2D space
    ///
    /// \param translate_x: translation along X axis
    /// \param translate_y: translation along Y axis
    /// \return 3x3 translation matrix
    ///
    template <typename T>
    constexpr matrix3x3<T> translate2D(T translate_x, T translate_y)
    {
        return matrix3x3<T>(1, 0, 0,
                         0, 1, 0,
                         translate_x, translate_y, 1);
    }


    /// \brief Construct 3x3 matrix for translation in homogeneous 2D space
    ///
    /// \param translate: translation along X and Y axes
    /// \return 3x3 translation matrix
    ///
    template <typename T>
    constexpr matrix3x3<T> translate2D(const vector2<T>& translate)
    {
        return translate2D<T>(translate.x, translate.y);
    }


    /// \brief Construct 3x3 matrix for 3D rotation around X axis
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 rotation matrix
    ///
    template <typename T>
    matrix3x3<T> rotate3Dx(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix3x3<T>(1, 0, 0,
                         0, _cos, _sin,
                         0, -_sin, _cos);
    }

    /// \brief Construct 4x4 matrix for rotation around X axis in 3d homogeneous space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 4x4 rotation matrix
    ///
    template <typename T>
    inline matrix4x4<T> rotate3Dhx(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix4x4<T>(1,  0, 0, 0,
                         0, _cos, _sin, 0,
                         0, -_sin, _cos, 0,
                         0,  0, 0, 1);
    }


    /// \brief Construct 3x3 matrix for 3D rotation around Y axis
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 rotation matrix
    ///
    template <typename T>
    matrix3x3<T> rotate3Dy(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix3x3<T>(_cos, 0, -_sin,
                         0, 1,  0,
                         _sin, 0,  _cos);
    }


    /// \brief Construct 4x4 matrix for rotation around Y axis in 3D homogeneous space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 4x4 rotation matrix
    ///
    template <typename T>
    matrix4x4<T> rotate3Dhy(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix4x4<T>(_cos, 0, -_sin, 0,
                         0, 1,  0, 0,
                         _sin, 0,  _cos, 0,
                         0, 0,  0, 1);
    }


    /// \brief Construct 3x3 matrix for 3D rotation around Z axis
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 rotation matrix
    ///
    template <typename T>
    matrix3x3<T> rotate3Dz(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix3x3<T>( _cos, _sin, 0,
                           -_sin, _cos, 0,
                           0, 0, 1);
    }


    /// \brief Construct 4x4 matrix for rotation around Z axis in homogeneous 3D space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 4x4 rotation matrix
    ///
    template <typename T>
    matrix4x4<T> rotate3Dhz(T angle)
    {
        const T _cos = cos(angle);
        const T _sin = sin(angle);

        return matrix4x4<T>( _cos, _sin, 0, 0,
                          -_sin, _cos, 0, 0,
                          0, 0, 1, 0,
                          0, 0, 0, 1);
    }


    /// \brief Construct 3x3 matrix for 3D scale
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \param scale_z: amount of scale along Z axis
    /// \return 3x3 scale matrix
    ///
    template <typename T>
    constexpr matrix3x3<T> scale3D(T scale_x, T scale_y, T scale_z)
    {
        return matrix3x3<T>(scale_x, 0, 0,
                         0, scale_y, 0,
                         0, 0, scale_z);
    }


    /// \brief Construct 3x3 matrix for 3D scale
    ///
    /// \param scale: amount of scale along X, Y, Z axes
    /// \return 3x3 scale matrix
    ///
    template <typename T>
    constexpr matrix3x3<T> scale3D(const vector3<T>& scale)
    {
        return scale3D<T>(scale.x, scale.y, scale.z);
    }


    /// \brief Construct 4x4 matrix for scale in homogeneous 3D space
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \param scale_z: amount of scale along Z axis
    /// \return 4x4 scale matrix
    ///
    template <typename T>
    constexpr matrix4x4<T> scale3Dh(T scale_x, T scale_y, T scale_z)
    {
        return matrix4x4<T>(scale_x, 0, 0, 0,
                         0, scale_y, 0, 0,
                         0, 0, scale_z, 0,
                         0, 0,    0, 1);
    }


    /// \brief Construct 4x4 matrix for scale in homogeneous 3D space
    ///
    /// \param scale: amount of scale along X, Y and Z axes
    /// \return 4x4 scale matrix
    ///
    template <typename T>
    constexpr matrix4x4<T> scale3Dh(const vector3<T>& scale)
    {
        return scale3Dh<T>(scale.x, scale.y, scale.z);
    }


    /// \brief Construct 4x4 matrix for translation in homogeneous 3D space
    ///
    /// \param x: translation along x axis
    /// \param y: translation along y axis
    /// \param z: translation along z axis
    /// \return 4x4 translation matrix
    ///
    template <typename T>
    constexpr matrix4x4<T> translate3D(T x, T y, T z)
    {
        return matrix4x4<T>(1, 0, 0, 0,
                         0, 1, 0, 0,
                         0, 0, 1, 0,
                            x,    y,    z, 1);
    }


    /// \brief Construct 4x4 matrix for translation in homogeneous 3D space
    ///
    /// \param translate: translation along X, Y, and Z axes
    /// \return 4x4 translation matrix
    ///
    template <typename T>
    constexpr matrix4x4<T> translate3D(const vector3<T>& translate)
    {
        return translate3D<T>(translate.x, translate.y, translate.z);
    }


    /// \brief Construct direction cosines matrix for
    /// rotation from one orthonormal basis to another (2D)
    ///
    /// \param p1: 1st basis vector of 1st basis
    /// \param q1: 2nd basis vector of 1st basis
    /// \param p2: 1st basis vector of 2nd basis
    /// \param q2: 2nd basis vector of 2nd basis
    /// \return matrix2x2
    ///
    /// With given basis vectors (that form two orthonormal basises
    /// B1={p1,q1}, B2={p2,q2}) expressed using the same coordinate
    /// space, constructs rotation matrix from B1 to B2.
    ///
    template <typename T>
    constexpr matrix2x2<T> direction_cosines(const vector2<T>& p1, const vector2<T>& q1,
                                const vector2<T>& p2, const vector2<T>& q2)
    {
        return matrix2x2<T>(p1.dot(p2), q1.dot(p2),
                          p1.dot(q2), q1.dot(q2));
    }


    /// \brief Construct direction cosines matrix for
    /// rotation from one orthonormal basis to another (3D)
    ///
    /// \param p1: 1st basis vector of 1st basis
    /// \param q1: 2nd basis vector of 1st basis
    /// \param r1: 3rd basis vector of 1st basis
    /// \param p2: 1st basis vector of 2nd basis
    /// \param q2: 2nd basis vector of 2nd basis
    /// \param r2: 3rd basis vector of 2nd basis
    /// \return matrix2x2
    ///
    /// With given basis vectors (that form two orthonormal basises
    /// B1={p1,q1,r1}, B2={p2,q2,r2}) expressed using the same coordinate
    /// space, constructs rotation matrix from B1 to B2.
    ///
    template <typename T>
    constexpr matrix3x3<T> direction_cosines(const vector3<T>& p1, const vector3<T>& q1, const vector3<T>& r1,
                                const vector3<T>& p2, const vector3<T>& q2, const vector3<T>& r2)
    {
        return matrix3x3<T>(p1.dot(p2), q1.dot(p2), r1.dot(p2),
                         p1.dot(q2), q1.dot(q2), r1.dot(q2),
                         p1.dot(r2), q1.dot(r2), r1.dot(r2));
    }
} // namespace math

#endif // MATRIX_TRANSFORMS_H_INCLUDED
