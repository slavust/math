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
    static const matrix2x2 reflect2Dx(-1.0f, 0.0f,
                                        0.0f, 1.0f);


    /// Reflection about Y axis matrix (2D)
    static const matrix2x2 reflect2Dy(1.0f, 0.0f,
                                       0.0f, -1.0f);


    /// Reflection about X axis matrix (2D homogeneous)
    static const matrix3x3 reflect2Dhx(-1.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f);


    /// Reflection about Y axis matrix (2D homogeneous)
    static const matrix3x3 reflect2Dhy(1.0f, 0.0f, 0.0f,
                                        0.0f, -1.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f);


    /// Reflection about XY plane matrix (3D)
    static const matrix3x3 reflect3Dxy(1.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, -1.0f);


    /// Reflection about XZ plane matrix (3D)
    static const matrix3x3 reflect3Dxz(1.0f, 0.0f, 0.0f,
                                       0.0f, -1.0f, 0.0f,
                                       0.0f, 0.0f, 1.0f);


    /// Reflection about YZ plane matrix (3D)
    static const matrix3x3 reflect3Dyz(-1.0f, 0.0f, 0.0f,
                                        0.0f, 1.0f, 0.0f,
                                        0.0f, 0.0f, 1.0f);


    /// Projection to XY plane matrix (3D)
    static const matrix3x3 project3Dxy(1.0f, 0.0f, 0.0f,
                                       0.0f, 1.0f, 0.0f,
                                       0.0f, 0.0f, 0.0f);


    /// Projection to XZ plane matrix (3D)
    static const matrix3x3 project3Dxz(1.0f, 0.0f, 0.0f,
                                       0.0f, 0.0f, 0.0f,
                                       0.0f, 0.0f, 1.0f);

    /// Projection to YZ plane matrix (3D)
    static const matrix3x3 project3Dyz(0.0f, 0.0f, 0.0f,
                                       0.0f, 1.0f, 0.0f,
                                       0.0f, 0.0f, 1.0f);


    // -------------- FUNCTIONS --------------

    /// \brief Construct 2x2 matrix for 2D rotation clockwise
    ///
    /// \param angle: amount of rotation in radians
    /// \return 2x2 rotation matrix
    ///
    inline matrix2x2 rotate2D(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix2x2(_cos, -_sin,
                          _sin, _cos);
    }

    /// \brief Construct 3x3 matrix for rotation clockwise in homogeneous 2D space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 2D rotation matrix
    ///
    inline matrix3x3 rotate2Dh(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix3x3(_cos, -_sin, 0.0f,
                         _sin, _cos, 0.0f,
                         0.0f, 0.0f, 1.0f);
    }


    /// \brief Construct 2x2 matrix for 2D scale
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \return 2x2 scale matrix
    ///
    inline matrix2x2 scale2D(real scale_x, real scale_y)
    {
        return matrix2x2(scale_x, 0.0f,
                         0.0f, scale_y);
    }


    /// \brief Construct 2x2 matrix for 2D scale
    ///
    /// \param scale: amount of scale along X and Y axes
    /// \return 2x2 scale matrix
    ///
    inline matrix2x2 scale2D(const vector2& scale)
    {
        return scale2D(scale.x, scale.y);
    }

    /// \brief Construct 3x3 matrix for scale in homogeneous 2D space
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \return 3x3 scale matrix
    ///
    inline matrix3x3 scale2Dh(real scale_x, real scale_y)
    {
        return matrix3x3(scale_x, 0.0f, 0.0f,
                         0.0f, scale_y, 0.0f,
                         0.0f,    0.0f, 1.0f);
    }


    /// \brief Construct 3x3 matrix for scale in homogeneous 2D space
    ///
    /// \param scale: amount of scale along X and Y axes
    /// \return 3x3 scale matrix
    ///
    inline matrix3x3 scale2Dh(const vector2& scale)
    {
        return scale2Dh(scale.x, scale.y);
    }


    /// \brief Construct 3x3 matrix for translation in homogeneous 2D space
    ///
    /// \param translate_x: translation along X axis
    /// \param translate_y: translation along Y axis
    /// \return 3x3 translation matrix
    ///
    inline matrix3x3 translate2D(real translate_x, real translate_y)
    {
        return matrix3x3(1.0f, 0.0f, 0.0f,
                         0.0f, 1.0f, 0.0f,
                         translate_x, translate_y, 1.0f);
    }


    /// \brief Construct 3x3 matrix for translation in homogeneous 2D space
    ///
    /// \param translate: translation along X and Y axes
    /// \return 3x3 translation matrix
    ///
    inline matrix3x3 translate2D(const vector2& translate)
    {
        return translate2D(translate.x, translate.y);
    }


    /// \brief Construct 3x3 matrix for 3D rotation around X axis
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 rotation matrix
    ///
    inline matrix3x3 rotate3Dx(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix3x3(1.0f, 0.0f, 0.0f,
                         0.0f, _cos, _sin,
                         0.0f, -_sin, _cos);
    }

    /// \brief Construct 4x4 matrix for rotation around X axis in 3d homogeneous space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 4x4 rotation matrix
    ///
    inline matrix4x4 rotate3Dhx(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix4x4(1.0f,  0.0f, 0.0f, 0.0f,
                         0.0f, _cos, _sin, 0.0f,
                         0.0f, -_sin, _cos, 0.0f,
                         0.0f,  0.0f, 0.0f, 1.0f);
    }


    /// \brief Construct 3x3 matrix for 3D rotation around Y axis
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 rotation matrix
    ///
    inline matrix3x3 rotate3Dy(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix3x3(_cos, 0.0f, -_sin,
                         0.0f, 1.0f,  0.0f,
                         _sin, 0.0f,  _cos);
    }


    /// \brief Construct 4x4 matrix for rotation around Y axis in 3D homogeneous space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 4x4 rotation matrix
    ///
    inline matrix4x4 rotate3Dhy(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix4x4(_cos, 0.0f, -_sin, 0.0f,
                         0.0f, 1.0f,  0.0f, 0.0f,
                         _sin, 0.0f,  _cos, 0.0f,
                         0.0f, 0.0f,  0.0f, 1.0f);
    }


    /// \brief Construct 3x3 matrix for 3D rotation around Z axis
    ///
    /// \param angle: amount of rotation in radians
    /// \return 3x3 rotation matrix
    ///
    inline matrix3x3 rotate3Dz(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix3x3( _cos, _sin, 0.0f,
                           -_sin, _cos, 0.0f,
                           0.0f, 0.0f, 1.0f);
    }


    /// \brief Construct 4x4 matrix for rotation around Z axis in homogeneous 3D space
    ///
    /// \param angle: amount of rotation in radians
    /// \return 4x4 rotation matrix
    ///
    inline matrix4x4 rotate3Dhz(real angle)
    {
        real _cos = cos(angle);
        real _sin = sin(angle);

        return matrix4x4( _cos, _sin, 0.0f, 0.0f,
                          -_sin, _cos, 0.0f, 0.0f,
                          0.0f, 0.0f, 1.0f, 0.0f,
                          0.0f, 0.0f, 0.0f, 1.0f);
    }


    /// \brief Construct 3x3 matrix for 3D scale
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \param scale_z: amount of scale along Z axis
    /// \return 3x3 scale matrix
    ///
    inline matrix3x3 scale3D(real scale_x, real scale_y, real scale_z)
    {
        return matrix3x3(scale_x, 0.0f, 0.0f,
                         0.0f, scale_y, 0.0f,
                         0.0f, 0.0f, scale_z);
    }


    /// \brief Construct 3x3 matrix for 3D scale
    ///
    /// \param scale: amount of scale along X, Y, Z axes
    /// \return 3x3 scale matrix
    ///
    inline matrix3x3 scale3D(const vector3& scale)
    {
        return scale3D(scale.x, scale.y, scale.z);
    }


    /// \brief Construct 4x4 matrix for scale in homogeneous 3D space
    ///
    /// \param scale_x: amount of scale along X axis
    /// \param scale_y: amount of scale along Y axis
    /// \param scale_z: amount of scale along Z axis
    /// \return 4x4 scale matrix
    ///
    inline matrix4x4 scale3Dh(real scale_x, real scale_y, real scale_z)
    {
        return matrix4x4(scale_x, 0.0f, 0.0f, 0.0f,
                         0.0f, scale_y, 0.0f, 0.0f,
                         0.0f, 0.0f, scale_z, 0.0f,
                         0.0f, 0.0f,    0.0f, 1.0f);
    }


    /// \brief Construct 4x4 matrix for scale in homogeneous 3D space
    ///
    /// \param scale: amount of scale along X, Y and Z axes
    /// \return 4x4 scale matrix
    ///
    inline matrix4x4 scale3Dh(const vector3& scale)
    {
        return scale3Dh(scale.x, scale.y, scale.z);
    }


    /// \brief Construct 4x4 matrix for translation in homogeneous 3D space
    ///
    /// \param x: translation along x axis
    /// \param y: translation along y axis
    /// \param z: translation along z axis
    /// \return 4x4 translation matrix
    ///
    inline matrix4x4 translate3D(real x, real y, real z)
    {
        return matrix4x4(1.0f, 0.0f, 0.0f, 0.0f,
                         0.0f, 1.0f, 0.0f, 0.0f,
                         0.0f, 0.0f, 1.0f, 0.0f,
                            x,    y,    z, 1.0f);
    }


    /// \brief Construct 4x4 matrix for translation in homogeneous 3D space
    ///
    /// \param translate: translation along X, Y, and Z axes
    /// \return 4x4 translation matrix
    ///
    inline matrix4x4 translate3D(const vector3& translate)
    {
        return translate3D(translate.x, translate.y, translate.z);
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
    inline matrix2x2 direction_cosines(const vector2& p1, const vector2& q1,
                                const vector2& p2, const vector2& q2)
    {
        return matrix2x2(p1.dot(p2), q1.dot(p2),
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
    inline matrix3x3 direction_cosines(const vector3& p1, const vector3& q1, const vector3& r1,
                                const vector3& p2, const vector3& q2, const vector3& r2)
    {
        return matrix3x3(p1.dot(p2), q1.dot(p2), r1.dot(p2),
                         p1.dot(q2), q1.dot(q2), r1.dot(q2),
                         p1.dot(r2), q1.dot(r2), r1.dot(r2));
    }
} // namespace math

#endif // MATRIX_TRANSFORMS_H_INCLUDED
