#ifndef MATRIX3X3_H_INCLUDED
#define MATRIX3X3_H_INCLUDED

#include "math_predefs.h"
#include "vector3.h"

namespace math
{
    class euler;
    class quaternion;

    /// \brief 3x3 Matrix
    ///
    /// Used for linear transformations in 3D space and
    /// affine transformations in 2D homogeneous space
    ///
    class matrix3x3
    {
    protected:
        real val[3][3]; ///< matrix elements

    public:
        static const matrix3x3 IDENTITY; ///< Identity matrix


        /// \brief Default constructor
        ///
        /// Lefts matrix elements unitialized
        ///
        matrix3x3();


        /// \brief Copy constructor
        ///
        /// \param src: source 3x3 matrix
        ///
        matrix3x3(const matrix3x3& src);

        /// \brief construct matrix from three column vectors as follows:
        /// [p.x q.x r.x
        ///  p.y q.y r.y
        ///  p.z q.z r.z]
        ///
        /// \param p: 1st vector
        /// \param q: 2nd vector
        /// \param r: 3rd vector
        ///
        matrix3x3(const vector3& p, const vector3& q, const vector3& r);

        /// \brief Initializes matrix elements with corresponding array elements
        ///
        /// \param src: 3x3 array
        ///
        matrix3x3(const real src[3][3]);


        /// \brief Initializes matrix with a00-a22 elements
        ///
        /// \param a00-a22: elements of matrix starting with zero
        ///
        matrix3x3(real a00, real a01, real a02,
                  real a10, real a11, real a12,
                  real a20, real a21, real a22);


        /// \brief Convert rotation matrix to Euler angles
        ///
        operator euler() const;


        /// \brief Convert rotation matrix to quaternion
        ///
        operator quaternion() const;


        /// \brief Constructs matrix transpose
        ///
        /// \return 3x3 matrix
        ///
        matrix3x3 transpose() const;


        /// \brief Computes determinant of matrix
        ///
        /// \return determinant
        ///
        real determinant() const;


        /// \brief Constructs classical adjoint
        ///
        /// \return 3x3 matrix
        ///
        /// Classical adjoint of matrix M is transpose of
        /// the matrix of cofactors of M.
        ///
        matrix3x3 adjoint() const;


        /// \brief Constructs inverse of matrix
        ///
        /// \return 3x3 matrix
        ///
        /// Inverse of matrix M is computed as M.adjoint()/M.determinant().
        /// If determinant of M is zero, then M is non-invertible matrix.
        /// In this case exception ET_NON_INVERTIBLE_MATRIX is occurred.
        ///
        matrix3x3 inverse() const
        {
            real d = determinant();

            if(d == 0.0f) MATH_EXCEPTION(ET_NON_INVERTIBLE_MATRIX);

            return adjoint() / d;
        }


        /// \brief Gram-Shmidt orthonormalization
        ///
        /// Biased towards OX orthonormalization
        ///
        void orthonormalize();


        matrix3x3& operator = (const real src[3][3])
        {
            memcpy(val, src, sizeof(real)*9);
            return *this;
        }

        matrix3x3& operator = (const matrix3x3& src)
        {
            memcpy(val, src.val, sizeof(real)*9);
            return *this;
        }

        bool operator == (const real op2[3][3]) const
        {
            return !memcmp(val, op2, sizeof(real)*9);
        }

        bool operator == (const matrix3x3& op2) const
        {
            return !memcmp(val, op2.val, sizeof(real)*9);
        }

        matrix3x3 operator * (real scalar) const;

        matrix3x3& operator *= (real scalar);

        matrix3x3 operator / (real scalar) const;

        matrix3x3& operator /= (real scalar);


        /// \brief Multiple by column vector
        ///
        /// \param v: 3D column vector
        /// \return 3D column vector
        ///
        vector3 operator * (const vector3& v) const;

        matrix3x3 operator * (const matrix3x3& b) const;

        matrix3x3& operator *= (const matrix3x3& b)
        {
            *this = *this * b;
            return *this;
        }

#ifdef MATH_CHECK_INDEX_BOUNDS
        ArrayHolder<real, 3> operator [] (size_t indx)
        {
            if(indx >= 3) MATH_EXCEPTION(ET_INDEX_OUT_OF_BOUNDS);
            return ArrayHolder<real, 3>(val[indx]);
        }

        ArrayHolder<const real, 3> operator [] (size_t indx) const
        {
            if(indx >= 3) MATH_EXCEPTION(ET_INDEX_OUT_OF_BOUNDS);
            return ArrayHolder<const real, 3>(val[indx]);
        }
#else
        real* operator [] (size_t indx)
        {
            return val[indx];
        }

        const real* operator [] (size_t indx) const
        {
            return val[indx];
        }
#endif


        /// \brief Pointer to array of elements
        ///
        /// \return real*
        ///
        /// Useful for copying with memcpy etc.
        ///
        real* ptr()
        {
            return reinterpret_cast<real*>(val);
        }


        /// \brief Pointer to constant array of elements
        ///
        /// \return real*
        ///
        /// Useful for copying with memcpy etc.
        ///
        const real* ptr() const
        {
            return reinterpret_cast<const real*>(val);
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

    inline matrix3x3 operator * (real scalar, const matrix3x3& m)
    {
        return m*scalar;
    }


    inline matrix3x3::matrix3x3()
    {
        //memcpy(val, IDENTITY.val, sizeof(real)*9);
    }

    inline matrix3x3::matrix3x3(const real src[3][3])
    {
        memcpy(val, src, sizeof(real)*9);
    }

    inline matrix3x3::matrix3x3(const matrix3x3& src)
    {
        memcpy(val, src.val, sizeof(real)*9);
    }

    inline matrix3x3::matrix3x3(const vector3& p, const vector3& q, const vector3& r)
    {
        val[0][0] = p.x;
        val[0][1] = p.y;
        val[0][2] = p.z;

        val[1][0] = q.x;
        val[1][1] = q.y;
        val[1][2] = q.z;

        val[2][0] = r.x;
        val[2][1] = r.y;
        val[2][2] = r.z;
    }

    inline matrix3x3::matrix3x3(real a00, real a01, real a02,
              real a10, real a11, real a12,
              real a20, real a21, real a22)
    {
        val[0][0] = a00;
        val[0][1] = a01;
        val[0][2] = a02;

        val[1][0] = a10;
        val[1][1] = a11;
        val[1][2] = a12;

        val[2][0] = a20;
        val[2][1] = a21;
        val[2][2] = a22;
    }

    inline vector3 matrix3x3::operator * (const vector3& v) const
    {
        return vector3(val[0][0]*v.x + val[1][0]*v.y + val[2][0]*v.z,
                       val[0][1]*v.x + val[1][1]*v.y + val[2][1]*v.z,
                       val[0][2]*v.x + val[1][2]*v.y + val[2][2]*v.z);
    }

    inline matrix3x3 matrix3x3::transpose() const
    {
        return matrix3x3(val[0][0], val[1][0], val[2][0],
                         val[0][1], val[1][1], val[2][1],
                         val[0][2], val[1][2], val[2][2]);
    }

    inline real matrix3x3::determinant() const
    {
        return val[0][0]*(val[1][1]*val[2][2] - val[2][1]*val[1][2])
                + val[1][0]*(val[2][1]*val[0][2] - val[0][1]*val[2][2])
                + val[2][0]*(val[0][1]*val[1][2] - val[1][1]*val[0][2]);
    }

} // namespace math

#endif // MATRIX3X3_H_INCLUDED
