#ifndef MATRIX4X4_H_INCLUDED
#define MATRIX4X4_H_INCLUDED

#include "math_predefs.h"
#include "vector4.h"

namespace math
{
    /// \brief 4x4 Matrix
    ///
    /// Used for affine transformations in 3D homogeneous space
    ///
    template <typename T>
    class matrix4x4
    {
    protected:
        std::array<std::array<T, 4>, 4> val; ///< matrix elements

    public:
          static constexpr matrix4x4<T> IDENTITY = matrix4x4<T>(1, 0, 0, 0,
                                                        0, 1, 0, 0,
                                                        0, 0, 1, 0,
                                                        0, 0, 0, 1); ///< Identity matrix

        /// \brief Default constructor
        ///
        /// Lefts matrix elements unitialized
        ///
        constexpr matrix4x4() : matrix4x4(IDENTITY)
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src: source 4x4 matrix
        ///
        constexpr matrix4x4(const matrix4x4<T>& src) : val(src.val)
        {
        }


        /// \brief Initializes matrix with a00-a22 elements
        ///
        /// \param a00-a33: elements of matrix starting with zero
        ///
        constexpr matrix4x4(T a00, T a01, T a02, T a03,
                            T a10, T a11, T a12, T a13,
                            T a20, T a21, T a22, T a23,
                            T a30, T a31, T a32, T a33)
            : val({a00, a01, a02, a03,
                   a10, a11, a12, a13,
                   a20, a21, a22, a23,
                   a30, a31, a32, a33})
        {
        }


        /// \brief Constructs matrix transpose
        ///
        /// \return 4x4 matrix
        ///
        constexpr matrix4x4<T> transposed() const
        {
            return matrix4x4<T>(
                val[0][0], val[1][0], val[2][0], val[3][0],
                val[0][1], val[1][1], val[2][1], val[3][1],
                val[0][2], val[1][2], val[2][2], val[3][2],
                val[0][3], val[1][3], val[2][3], val[3][3]
            );
        }


        /// \brief Computes determinant of matrix
        ///
        /// \return determinant
        ///
        constexpr T determinant() const
        {
            return val[0][0] * (val[1][1]*(val[2][2]*val[3][3] - val[3][2]*val[2][3])
                            + val[2][1]*(val[3][2]*val[1][3] - val[1][2]*val[3][3])
                            + val[3][1]*(val[1][2]*val[2][3] - val[2][2]*val[1][3]))

             - val[1][0] * (val[0][1]*(val[2][2]*val[3][3] - val[3][2]*val[2][3])
                            + val[2][1]*(val[3][2]*val[0][3] - val[0][2]*val[3][3])
                            + val[3][1]*(val[0][2]*val[2][3] - val[2][2]*val[0][3]))

             + val[2][0] * (val[0][1]*(val[1][2]*val[3][3] - val[3][2]*val[1][3])
                            + val[1][1]*(val[3][2]*val[0][3] - val[0][2]*val[3][3])
                            + val[3][1]*(val[0][2]*val[1][3] - val[1][2]*val[0][3]))

             - val[3][0] * (val[0][1]*(val[1][2]*val[2][3] - val[2][2]*val[1][3])
                            + val[1][1]*(val[2][2]*val[0][3] - val[0][2]*val[2][3])
                            + val[2][1]*(val[0][2]*val[1][3] - val[1][2]*val[0][3]));

        }


        /// \brief Constructs classical adjoint
        ///
        /// \return 4x4 matrix
        ///
        /// Classical adjoint of matrix M is transpose of
        /// the matrix of cofactors of M.
        ///
        constexpr matrix4x4<T> adjoint() const
        {
        return matrix4x4<T>(
             val[1][1]*val[2][2]*val[3][3] + val[2][1]*val[3][2]*val[1][3] + val[3][1]*val[1][2]*val[2][3] - val[3][1]*val[2][2]*val[1][3] - val[2][1]*val[1][2]*val[3][3] - val[1][1]*val[3][2]*val[2][3],
            -val[0][1]*val[2][2]*val[3][3] - val[2][1]*val[3][2]*val[0][3] - val[3][1]*val[0][2]*val[2][3] + val[3][1]*val[2][2]*val[0][3] + val[2][1]*val[0][2]*val[3][3] + val[0][1]*val[3][2]*val[2][3],
             val[0][1]*val[1][2]*val[3][3] + val[1][1]*val[3][2]*val[0][3] + val[3][1]*val[0][2]*val[1][3] - val[3][1]*val[1][2]*val[0][3] - val[1][1]*val[0][2]*val[3][3] - val[0][1]*val[3][2]*val[1][3],
            -val[0][1]*val[1][2]*val[2][3] - val[1][1]*val[2][2]*val[0][3] - val[2][1]*val[0][2]*val[1][3] + val[2][1]*val[1][2]*val[0][3] + val[1][1]*val[0][2]*val[2][3] + val[0][1]*val[2][2]*val[1][3],

            -val[1][0]*val[2][2]*val[3][3] - val[2][0]*val[3][2]*val[1][3] - val[3][0]*val[1][2]*val[2][3] + val[3][0]*val[2][2]*val[1][3] + val[2][0]*val[1][2]*val[3][3] + val[1][0]*val[3][2]*val[2][3],
             val[0][0]*val[2][2]*val[3][3] + val[2][0]*val[3][2]*val[0][3] + val[3][0]*val[0][2]*val[2][3] - val[3][0]*val[2][2]*val[0][3] - val[2][0]*val[0][2]*val[3][3] - val[0][0]*val[3][2]*val[2][3],
            -val[0][0]*val[1][2]*val[3][3] - val[1][0]*val[3][2]*val[0][3] - val[3][0]*val[0][2]*val[1][3] + val[3][0]*val[1][2]*val[0][3] + val[1][0]*val[0][2]*val[3][3] + val[0][0]*val[3][2]*val[1][3],
             val[0][0]*val[1][2]*val[2][3] + val[1][0]*val[2][2]*val[0][3] + val[2][0]*val[0][2]*val[1][3] - val[2][0]*val[1][2]*val[0][3] - val[1][0]*val[0][2]*val[2][3] - val[0][0]*val[2][2]*val[1][3],

             val[1][0]*val[2][1]*val[3][3] + val[2][0]*val[3][1]*val[1][3] + val[3][0]*val[1][1]*val[2][3] - val[3][0]*val[2][1]*val[1][3] - val[2][0]*val[1][1]*val[3][3] - val[1][0]*val[3][1]*val[2][3],
            -val[0][0]*val[2][1]*val[3][3] - val[2][0]*val[3][1]*val[0][3] - val[3][0]*val[0][1]*val[2][3] + val[3][0]*val[2][1]*val[0][3] + val[2][0]*val[0][1]*val[3][3] + val[0][0]*val[3][1]*val[2][3],
             val[0][0]*val[1][1]*val[3][3] + val[1][0]*val[3][1]*val[0][3] + val[3][0]*val[0][1]*val[1][3] - val[3][0]*val[1][1]*val[0][3] - val[1][0]*val[0][1]*val[3][3] - val[0][0]*val[3][1]*val[1][3],
            -val[0][0]*val[1][1]*val[2][3] - val[1][0]*val[2][1]*val[0][3] - val[2][0]*val[0][1]*val[1][3] + val[2][0]*val[1][1]*val[0][3] + val[1][0]*val[0][1]*val[2][3] + val[0][0]*val[2][1]*val[1][3],

            -val[1][0]*val[2][1]*val[3][2] - val[2][0]*val[3][1]*val[1][2] - val[3][0]*val[1][1]*val[2][2] + val[3][0]*val[2][1]*val[1][2] + val[2][0]*val[1][1]*val[3][2] + val[1][0]*val[3][1]*val[2][2],
             val[0][0]*val[2][1]*val[3][2] + val[2][0]*val[3][1]*val[0][2] + val[3][0]*val[0][1]*val[2][2] - val[3][0]*val[2][1]*val[0][2] - val[2][0]*val[0][1]*val[3][2] - val[0][0]*val[3][1]*val[2][2],
            -val[0][0]*val[1][1]*val[3][2] - val[1][0]*val[3][1]*val[0][2] - val[3][0]*val[0][1]*val[1][2] + val[3][0]*val[1][1]*val[0][2] + val[1][0]*val[0][1]*val[3][2] + val[0][0]*val[3][1]*val[1][2],
             val[0][0]*val[1][1]*val[2][2] + val[1][0]*val[2][1]*val[0][2] + val[2][0]*val[0][1]*val[1][2] - val[2][0]*val[1][1]*val[0][2] - val[1][0]*val[0][1]*val[2][2] - val[0][0]*val[2][1]*val[1][2]);
        }


        /// \brief Constructs inverse of matrix
        ///
        /// \return 4x4 matrix
        ///
        /// Inverse of matrix M is computed as M.adjoint()/M.determinant().
        /// If determinant of M is zero, then M is non-invertible matrix.
        /// In this case exception ET_NON_INVERTIBLE_MATRIX is occurred.
        ///
        constexpr matrix4x4<T> inverse() const
        {
            const T d = determinant();
            return adjoint() / d;
        }


        matrix4x4<T>& operator = (const std::array<std::array<T, 4>, 4>& src)
        {
            val = src.val;
            return *this;
        }

        matrix4x4& operator = (const matrix4x4& src)
        {
            val = src.val;
            return *this;
        }

        constexpr bool operator == (const std::array<std::array<T, 4>, 4>& op2) const
        {
            return val == op2;
        }

        constexpr bool operator == (const matrix4x4& op2) const
        {
            return val == op2.val;
        }

        constexpr matrix4x4<T> operator * (T scalar) const
        {
            return matrix4x4<T>(
                scalar * val[0][0], scalar * val[0][1], scalar * val[0][2], scalar * val[0][3],
                scalar * val[1][0], scalar * val[1][1], scalar * val[1][2], scalar * val[1][3],
                scalar * val[2][0], scalar * val[2][1], scalar * val[2][2], scalar * val[2][3],
                scalar * val[3][0], scalar * val[3][1], scalar * val[3][2], scalar * val[3][3]
            );
        }

        matrix4x4& operator *= (T scalar)
        {
            for(auto& arr : val)
                for(T& e : arr)
                    e *= scalar;
        }

        constexpr matrix4x4<T> operator / (T scalar) const
        {
            return matrix4x4<T>(
                val[0][0] / scalar, val[0][1] / scalar, val[0][2] / scalar, val[0][3] / scalar,
                val[1][0] / scalar, val[1][1] / scalar, val[1][2] / scalar, val[1][3] / scalar,
                val[2][0] / scalar, val[2][1] / scalar, val[2][2] / scalar, val[2][3] / scalar,
                val[3][0] / scalar, val[3][1] / scalar, val[3][2] / scalar, val[3][3] / scalar
            );
        }

        matrix4x4& operator /= (T scalar)
        {
            for(auto& arr : val)
                for (T& e : arr)
                    e /= scalar;
        }


        /// \brief Multiple by column vector
        ///
        /// \param v: 3D column vector
        /// \return 3D column vector
        ///
        constexpr vector4<T> operator * (const vector4<T>& v) const
        {
            return vector4<T>(
                val[0][0]*v.x + val[1][0]*v.y + val[2][0]*v.z + val[3][0]*v.w,
                val[0][1]*v.x + val[1][1]*v.y + val[2][1]*v.z + val[3][1]*v.w,
                val[0][2]*v.x + val[1][2]*v.y + val[2][2]*v.z + val[3][2]*v.w,
                val[0][3]*v.x + val[1][3]*v.y + val[2][3]*v.z + val[3][3]*v.w);
        }

        constexpr matrix4x4<T> operator * (const matrix4x4<T>& m) const
        {
        return matrix4x4(
                          m[0][0]*val[0][0] + m[1][0]*val[0][1] + m[2][0]*val[0][2] + m[3][0]*val[0][3],
                          m[0][1]*val[0][0] + m[1][1]*val[0][1] + m[2][1]*val[0][2] + m[3][1]*val[0][3],
                          m[0][2]*val[0][0] + m[1][2]*val[0][1] + m[2][2]*val[0][2] + m[3][2]*val[0][3],
                          m[0][3]*val[0][0] + m[1][3]*val[0][1] + m[2][3]*val[0][2] + m[3][3]*val[0][3],

                          m[0][0]*val[1][0] + m[1][0]*val[1][1] + m[2][0]*val[1][2] + m[3][0]*val[1][3],
                          m[0][1]*val[1][0] + m[1][1]*val[1][1] + m[2][1]*val[1][2] + m[3][1]*val[1][3],
                          m[0][2]*val[1][0] + m[1][2]*val[1][1] + m[2][2]*val[1][2] + m[3][2]*val[1][3],
                          m[0][3]*val[1][0] + m[1][3]*val[1][1] + m[2][3]*val[1][2] + m[3][3]*val[1][3],

                          m[0][0]*val[2][0] + m[1][0]*val[2][1] + m[2][0]*val[2][2] + m[3][0]*val[2][3],
                          m[0][1]*val[2][0] + m[1][1]*val[2][1] + m[2][1]*val[2][2] + m[3][1]*val[2][3],
                          m[0][2]*val[2][0] + m[1][2]*val[2][1] + m[2][2]*val[2][2] + m[3][2]*val[2][3],
                          m[0][3]*val[2][0] + m[1][3]*val[2][1] + m[2][3]*val[2][2] + m[3][3]*val[2][3],

                          m[0][0]*val[3][0] + m[1][0]*val[3][1] + m[2][0]*val[3][2] + m[3][0]*val[3][3],
                          m[0][1]*val[3][0] + m[1][1]*val[3][1] + m[2][1]*val[3][2] + m[3][1]*val[3][3],
                          m[0][2]*val[3][0] + m[1][2]*val[3][1] + m[2][2]*val[3][2] + m[3][2]*val[3][3],
                          m[0][3]*val[3][0] + m[1][3]*val[3][1] + m[2][3]*val[3][2] + m[3][3]*val[3][3]);

        }

        matrix4x4& operator *= (const matrix4x4& b)
        {
            *this = *this * b;
            return *this;
        }

        std::array<T, 4>& operator [] (size_t indx)
        {
            return val[indx];
        }

        constexpr const std::array<T, 4>& operator [] (size_t indx) const
        {
            return val[indx];
        }
    };
} // namespace math

#endif // MATRIX4X4_H_INCLUDED
