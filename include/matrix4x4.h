#ifndef MATRIX4X4_H_INCLUDED
#define MATRIX4X4_H_INCLUDED

#include "math_predefs.h"
#include "vector4.h"

namespace math
{
    class euler;
    class quaternion;

    /// \brief 4x4 Matrix
    ///
    /// Used for affine transformations in 3D homogeneous space
    ///
    class matrix4x4
    {
    protected:
        real val[4][4]; ///< matrix elements

    public:
        static const matrix4x4 IDENTITY; ///< Identity matrix

        /// \brief Default constructor
        ///
        /// Lefts matrix elements unitialized
        ///
        matrix4x4();


        /// \brief Copy constructor
        ///
        /// \param src: source 4x4 matrix
        ///
        matrix4x4(const matrix4x4& src);


        /// \brief Initializes matrix elements with corresponding array elements
        ///
        /// \param src: 4x4 array
        ///
        matrix4x4(const real src[4][4]);


        /// \brief Initializes matrix with a00-a22 elements
        ///
        /// \param a00-a33: elements of matrix starting with zero
        ///
        matrix4x4(real a00, real a10, real a20, real a30,
                  real a01, real a11, real a21, real a31,
                  real a02, real a12, real a22, real a32,
                  real a03, real a13, real a23, real a33);


        /// \brief Convertion from rotation matrix to Euler angles
        ///
        operator euler () const;


        /// \brief Convertion from rotation matrix to quaternion
        ///
        operator quaternion () const;


        /// \brief Constructs matrix transpose
        ///
        /// \return 4x4 matrix
        ///
        matrix4x4 transpose() const;


        /// \brief Computes determinant of matrix
        ///
        /// \return determinant
        ///
        real determinant() const;


        /// \brief Constructs classical adjoint
        ///
        /// \return 4x4 matrix
        ///
        /// Classical adjoint of matrix M is transpose of
        /// the matrix of cofactors of M.
        ///
        matrix4x4 adjoint() const;


        /// \brief Constructs inverse of matrix
        ///
        /// \return 4x4 matrix
        ///
        /// Inverse of matrix M is computed as M.adjoint()/M.determinant().
        /// If determinant of M is zero, then M is non-invertible matrix.
        /// In this case exception ET_NON_INVERTIBLE_MATRIX is occurred.
        ///
        matrix4x4 inverse() const
        {
            real d = determinant();
            return adjoint() / d;
        }


        matrix4x4& operator = (const real src[4][4])
        {
            memcpy(val, src, sizeof(real)*16);
            return *this;
        }

        matrix4x4& operator = (const matrix4x4& src)
        {
            memcpy(val, src.val, sizeof(real)*16);
            return *this;
        }

        bool operator == (const real op2[4][4]) const
        {
            return !memcmp(val, op2, sizeof(real)*16);
        }

        bool operator == (const matrix4x4& op2) const
        {
            return !memcmp(val, op2.val, sizeof(real)*16);
        }

        matrix4x4 operator * (real scalar) const;

        matrix4x4& operator *= (real scalar);

        matrix4x4 operator / (real scalar) const;

        matrix4x4& operator /= (real scalar);


        /// \brief Multiple by column vector
        ///
        /// \param v: 3D column vector
        /// \return 3D column vector
        ///
        vector4 operator * (const vector4& v) const;

        matrix4x4 operator * (const matrix4x4& b) const;

        matrix4x4& operator *= (const matrix4x4& b)
        {
            *this = *this * b;

            return *this;
        }

        real* operator [] (size_t indx)
        {
            return val[indx];
        }

        const real* operator [] (size_t indx) const
        {
            return val[indx];
        }


        /// \brief Convert rotation matrix to Euler angles
        ///
        /// \return euler
        ///
        euler toEuler() const;


        /// \brief Convert rotation matrix to quaternion
        ///
        /// \return quaternion
        ///
        quaternion toQuaternion() const;
    };

    inline matrix4x4 operator * (real scalar, const matrix4x4& m)
    {
        return m*scalar;
    }


    inline matrix4x4::matrix4x4()
    {
        //memcpy(val, IDENTITY.val, sizeof(real)*16);
    }

    inline matrix4x4::matrix4x4(const real src[4][4])
    {
        memcpy(val, src, sizeof(real)*16);
    }

    inline matrix4x4::matrix4x4(const matrix4x4& src)
    {
        memcpy(val, src.val, sizeof(real)*16);
    }

    inline matrix4x4::matrix4x4(real a00, real a01, real a02, real a03,
                                real a10, real a11, real a12, real a13,
                                real a20, real a21, real a22, real a23,
                                real a30, real a31, real a32, real a33)
    {
        val[0][0] = a00;
        val[0][1] = a01;
        val[0][2] = a02;
        val[0][3] = a03;

        val[1][0] = a10;
        val[1][1] = a11;
        val[1][2] = a12;
        val[1][3] = a13;

        val[2][0] = a20;
        val[2][1] = a21;
        val[2][2] = a22;
        val[2][3] = a23;

        val[3][0] = a30;
        val[3][1] = a31;
        val[3][2] = a32;
        val[3][3] = a33;
    }

    inline vector4 matrix4x4::operator * (const vector4& v) const
    {
        return vector4(
            val[0][0]*v.x + val[1][0]*v.y + val[2][0]*v.z + val[3][0]*v.w,
            val[0][1]*v.x + val[1][1]*v.y + val[2][1]*v.z + val[3][1]*v.w,
            val[0][2]*v.x + val[1][2]*v.y + val[2][2]*v.z + val[3][2]*v.w,
            val[0][3]*v.x + val[1][3]*v.y + val[2][3]*v.z + val[3][3]*v.w);
    }

    inline matrix4x4 matrix4x4::transpose() const
    {
        return matrix4x4(val[0][0], val[1][0], val[2][0], val[3][0],
                         val[0][1], val[1][1], val[2][1], val[3][1],
                         val[0][2], val[1][2], val[2][2], val[3][2],
                         val[0][3], val[1][3], val[2][3], val[3][3]);
    }

} // namespace math

#endif // MATRIX4X4_H_INCLUDED
