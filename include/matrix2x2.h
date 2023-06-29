#ifndef MATRIX2X2_H_INCLUDED
#define MATRIX2X2_H_INCLUDED

#include "math_predefs.h"
#include "vector2.h"


namespace math
{
    /// \brief 2x2 matrix
    ///
    /// Used for 2D linear transformations
    ///
    class matrix2x2
    {
    protected:
        real val[2][2]; ///< matrix elements

    public:
        static const matrix2x2 IDENTITY; ///< Identity matrix


        /// \brief Default constructor
        ///
        /// Initializes matrix elements with zero
        ///
        matrix2x2();


        /// \brief Copy constructor
        ///
        /// \param src: source 2x2 matrix
        ///
        matrix2x2(const matrix2x2& src);

        /// \brief construct matrix from two column vectors as follows:
        /// [p.x q.x
        ///  p.y q.y]
        ///
        /// \param p: 1st vector
        /// \param q: 2nd vector
        ///
        matrix2x2(const vector2& p, const vector2& q);

        /// \brief Initializes matrix elements with corresponding array elements
        ///
        /// \param src: 2x2 array
        ///
        matrix2x2(const real src[2][2]);


        /// \brief Initializes matrix with a00-a11 elements
        ///
        /// \param a00-a11: elements of matrix starting with zero
        ///
        matrix2x2(real a00, real a01,
                  real a10, real a11);


        /// \brief Constructs matrix transpose
        ///
        /// \return 2x2 matrix
        ///
        matrix2x2 transpose() const;


        /// \brief Computes determinant of matrix
        ///
        /// \return determinant
        ///
        real determinant() const
        {
            return val[0][0]*val[1][1] - val[0][1]*val[1][0];
        }


        /// \brief Constructs classical adjoint
        ///
        /// \return 2x2 matrix
        ///
        /// Classical adjoint of matrix M is transpose of
        /// the matrix of cofactors of M.
        ///
        matrix2x2 adjoint() const
        {
            return matrix2x2(val[1][1], -val[0][1],
                            -val[1][0],  val[0][0]);
        }


        /// \brief Constructs inverse of matrix
        ///
        /// \return 2x2 matrix
        ///
        /// Inverse of matrix M is computed as M.adjoint()/M.determinant().
        /// If determinant of M is zero, then M is non-invertible matrix.
        /// In this case exception ET_NON_INVERTIBLE_MATRIX is occurred.
        ///
        matrix2x2 inverse() const
        {
            real d = determinant();
            return adjoint() / d;
        }


        /// \brief Gram-Shmidt orthonormalization
        ///
        /// Biased towards OX orthonormalization
        ///
        void orthonormalize();

        matrix2x2 operator * (real scalar) const;

        matrix2x2& operator *= (real scalar);

        matrix2x2 operator / (real scalar) const;

        matrix2x2& operator /= (real scalar);


        /// \brief Multiple by column vector
        ///
        /// \param v: 2D column vector
        /// \return 2D column vector
        ///
        inline vector2 operator * (const vector2& v) const
        {
            return vector2(val[0][0] * v.x + val[1][0] * v.y,
                           val[0][1] * v.x + val[1][1] * v.y);
        }


        matrix2x2 operator * (const matrix2x2& b) const;

        matrix2x2& operator *= (const matrix2x2& b)
        {
            *this = *this * b;
            return *this;
        }

        matrix2x2& operator = (real src[2][2])
        {
            memcpy(val, src, sizeof(real)*4);
            return *this;
        }

        matrix2x2& operator = (const matrix2x2& src)
        {
            memcpy(val, src.val, sizeof(real)*4);
            return *this;
        }

        bool operator == (const real op2[2][2]) const
        {
            return !memcmp(val, op2, sizeof(real)*4);
        }

        bool operator == (const matrix2x2& op2) const
        {
            return !memcmp(val, op2.val, sizeof(real)*4);
        }
        real* operator [] (size_t indx)
        {
            return val[indx];
        }

        const real* operator [] (size_t indx) const
        {
            return val[indx];
        }

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
    };

    inline matrix2x2 operator * (real scalar, const matrix2x2& m)
    {
        return m*scalar;
    }


    inline matrix2x2::matrix2x2()
    {
        memcpy(val, IDENTITY.val, sizeof(real)*4);
    }

    inline matrix2x2::matrix2x2(const vector2& p, const vector2& q)
    {
        val[0][0] = p.x;
        val[1][0] = p.y;

        val[0][1] = q.x;
        val[1][1] = q.y;
    }

    inline matrix2x2::matrix2x2(const real src[2][2])
    {
        memcpy(val, src, sizeof(real)*4);
    }

    inline matrix2x2::matrix2x2(real a00, real a01,
                                real a10, real a11)
    {
        val[0][0] = a00;
        val[0][1] = a01;
        val[1][0] = a10;
        val[1][1] = a11;
    }

    inline matrix2x2 matrix2x2::transpose() const
    {
        return matrix2x2(val[0][0], val[1][0], val[0][1], val[1][1]);
    }

    inline matrix2x2 matrix2x2::operator * (real scalar) const
    {
        return matrix2x2(val[0][0]*scalar, val[0][1]*scalar,
                         val[1][0]*scalar, val[1][1]*scalar);
    }

    inline matrix2x2& matrix2x2::operator *= (real scalar)
    {
        val[0][0] *= scalar;
        val[0][1] *= scalar;
        val[1][0] *= scalar;
        val[1][1] *= scalar;

        return *this;
    }

    inline matrix2x2 matrix2x2::operator / (real scalar) const
    {
#ifdef MATH_CHECK_DIVISION
            if(scalar == 0.0f) MATH_EXCEPTION(ET_DIVIDE_BY_ZERO);
#endif // MATH_CHECK_DIVISION
        return matrix2x2(val[0][0] / scalar, val[0][1] / scalar,
                         val[1][0] / scalar, val[1][1] / scalar);
    }

    inline matrix2x2& matrix2x2::operator /= (real scalar)
    {
#ifdef MATH_CHECK_DIVISION
            if(scalar == 0.0f) MATH_EXCEPTION(ET_DIVIDE_BY_ZERO);
#endif // MATH_CHECK_DIVISION
        val[0][0] /= scalar;
        val[0][1] /= scalar;
        val[1][0] /= scalar;
        val[1][1] /= scalar;

        return *this;
    }

    inline matrix2x2 matrix2x2::operator * (const matrix2x2& b) const
    {
        return matrix2x2(val[0][0] * b[0][0] + val[1][0] * b[0][1],
                         val[0][1] * b[0][0] + val[1][1] * b[0][1],
                         val[0][0] * b[1][0] + val[1][0] * b[1][1],
                         val[0][1] * b[1][0] + val[1][1] * b[1][1]);
    }

} // namespace math

#endif // MATRIX2X2_H_INCLUDED
