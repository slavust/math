#ifndef MATRIX3X3_H_INCLUDED
#define MATRIX3X3_H_INCLUDED

#include "math_predefs.h"
#include "vector3.h"

#include <array>

namespace math
{
    /// \brief 3x3 Matrix
    ///
    /// Used for linear transformations in 3D space and
    /// affine transformations in 2D homogeneous space
    ///
    template <typename T>
    class matrix3x3
    {
    protected:
        std::array<std::array<T, 3>, 3>  val; ///< matrix elements

    public:
        static constexpr matrix3x3<T> IDENTITY {1, 0, 0, 0, 1, 0, 0, 0, 1}; ///< Identity matrix


        /// \brief Default constructor
        ///
        /// Lefts matrix elements unitialized
        ///
        constexpr matrix3x3(): matrix3x3(IDENTITY)
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src: source 3x3 matrix
        ///
        constexpr matrix3x3(const matrix3x3<T>& src) : val(src.val)
        {
        }

        /// \brief construct matrix from three column vectors as follows:
        /// [p.x q.x r.x
        ///  p.y q.y r.y
        ///  p.z q.z r.z]
        ///
        /// \param p: 1st vector
        /// \param q: 2nd vector
        /// \param r: 3rd vector
        ///
        constexpr matrix3x3(const vector3<T>& p, const vector3<T>& q, const vector3<T>& r)
        : val({{p.x, p.y, p.z}, {q.x, q.y, q.z}, {r.x, r.y, r.z}})
        {
        }


        /// \brief Initializes matrix with a00-a22 elements
        ///
        /// \param a00-a22: elements of matrix starting with zero
        ///
        constexpr matrix3x3(T a00, T a01, T a02,
                  T a10, T a11, T a12,
                  T a20, T a21, T a22)
                  : val({a00, a01, a02, a10, a11,a12, a20, a21, a22})
        {
        }


        /// \brief Constructs matrix transpose
        ///
        /// \return 3x3 matrix
        ///
        constexpr matrix3x3<T> transposed() const
        {
        return matrix3x3<T>(val[0][0], val[1][0], val[2][0],
                         val[0][1], val[1][1], val[2][1],
                         val[0][2], val[1][2], val[2][2]);
        }


        /// \brief Computes determinant of matrix
        ///
        /// \return determinant
        ///
        constexpr T determinant() const
        {
            return val[0][0]*(val[1][1]*val[2][2] - val[2][1]*val[1][2])
                + val[1][0]*(val[2][1]*val[0][2] - val[0][1]*val[2][2])
                + val[2][0]*(val[0][1]*val[1][2] - val[1][1]*val[0][2]);
        }


        /// \brief Constructs classical adjoint
        ///
        /// \return 3x3 matrix
        ///
        /// Classical adjoint of matrix M is transpose of
        /// the matrix of cofactors of M.
        ///
        constexpr matrix3x3 adjoint() const
        {
            return matrix3x3<T>(val[1][1]*val[2][2] - val[2][1]*val[1][2],
                            val[2][1]*val[0][2] - val[0][1]*val[2][2],
                            val[0][1]*val[1][2] - val[1][1]*val[0][2],

                            val[2][0]*val[1][2] - val[1][0]*val[2][2],
                            val[0][0]*val[2][2] - val[2][0]*val[0][2],
                            val[1][0]*val[0][2] - val[0][0]*val[1][2],

                            val[1][0]*val[2][1] - val[2][0]*val[1][1],
                            val[2][0]*val[0][1] - val[0][0]*val[2][1],
                            val[0][0]*val[1][1] - val[1][0]*val[0][1]);
        }


        /// \brief Constructs inverse of matrix
        ///
        /// \return 3x3 matrix
        ///
        /// Inverse of matrix M is computed as M.adjoint()/M.determinant().
        /// If determinant of M is zero, then M is non-invertible matrix.
        /// In this case exception ET_NON_INVERTIBLE_MATRIX is occurred.
        ///
        constexpr matrix3x3<T> inverse() const
        {
            const T d = determinant();
            return adjoint() / d;
        }


        /// \brief Gram-Shmidt orthonormalization
        ///
        /// Biased towards OX orthonormalization
        ///
        constexpr matrix3x3<T> orthonormalized() const
        {
            const vector3<T> p(val[0][0], val[0][1], val[0][2]);
            const vector3<T> q(val[1][0], val[1][1], val[1][2]);
            const vector3<T> r(val[2][0], val[2][1], val[2][2]);

            const auto pn = p.normalized();
            const auto qn = (q - pn.dot(q)*pn).normalized();
            const auto rn = (r - pn.dot(r)*pn - qn.dot(r)*qn).normalized();

            return matrix3x3<T>(pn, qn, rn);
        }

        matrix3x3<T>& operator = (const matrix3x3<T>& src)
        {
            val = src.val;
            return *this;
        }

        constexpr bool operator == (const matrix3x3<T>& op2) const
        {
            return val == op2.val;
        }

        constexpr matrix3x3<T> operator * (T scalar) const
        {
            return matrix3x3<T>(
                scalar * val[0][0], scalar * val[0][1], scalar * val[0][2],
                scalar * val[1][0], scalar * val[1][1], scalar * val[1][2],
                scalar * val[2][0], scalar * val[2][1], scalar * val[2][2]
            );
        }

        matrix3x3<T>& operator *= (T scalar)
        {
            for(auto& arr : val)
                for(auto& e : arr)
                    e *= scalar;
                return *this;
        }

        constexpr matrix3x3<T> operator / (T scalar) const
        {
            return matrix3x3<T>(
                val[0][0] / scalar, val[0][1] / scalar, val[0][2] / scalar,
                val[1][0] / scalar, val[1][1] / scalar, val[1][2] / scalar,
                val[2][0] / scalar, val[2][1] / scalar, val[2][2] / scalar
            );
        }

        matrix3x3<T>& operator /= (T scalar)
        {
            for(auto& arr : val)
                for(auto& e : arr)
                    e /= scalar;
        }


        /// \brief Multiple by column vector
        ///
        /// \param v: 3D column vector
        /// \return 3D column vector
        ///
        constexpr vector3<T> operator * (const vector3<T>& v) const
        {
            return vector3<T>(val[0][0]*v.x + val[1][0]*v.y + val[2][0]*v.z,
                        val[0][1]*v.x + val[1][1]*v.y + val[2][1]*v.z,
                        val[0][2]*v.x + val[1][2]*v.y + val[2][2]*v.z);
        }

        constexpr matrix3x3<T> operator * (const matrix3x3<T>& m) const
        {
            return matrix3x3<T>(m[0][0]*val[0][0] + m[1][0]*val[0][1] + m[2][0]*val[0][2],
                            m[0][1]*val[0][0] + m[1][1]*val[0][1] + m[2][1]*val[0][2],
                            m[0][2]*val[0][0] + m[1][2]*val[0][1] + m[2][2]*val[0][2],

                            m[0][0]*val[1][0] + m[1][0]*val[1][1] + m[2][0]*val[1][2],
                            m[0][1]*val[1][0] + m[1][1]*val[1][1] + m[2][1]*val[1][2],
                            m[0][2]*val[1][0] + m[1][2]*val[1][1] + m[2][2]*val[1][2],

                            m[0][0]*val[2][0] + m[1][0]*val[2][1] + m[2][0]*val[2][2],
                            m[0][1]*val[2][0] + m[1][1]*val[2][1] + m[2][1]*val[2][2],
                            m[0][2]*val[2][0] + m[1][2]*val[2][1] + m[2][2]*val[2][2]
                            );
        }

        matrix3x3<T>& operator *= (const matrix3x3<T>& b)
        {
            *this = *this * b;
            return *this;
        }

        std::array<T, 3>& operator [] (size_t indx)
        {
            return val[indx];
        }

        constexpr const std::array<T, 3>& operator [] (size_t indx) const
        {
            return val[indx];
        }
    };
} // namespace math

#endif // MATRIX3X3_H_INCLUDED
