#ifndef MATRIX2X2_H_INCLUDED
#define MATRIX2X2_H_INCLUDED

#include "math_predefs.h"
#include "vector2.h"
#include <array>


namespace math
{
    /// \brief 2x2 matrix
    ///
    /// Used for 2D linear transformations
    ///
    template <typename T>
    class matrix2x2
    {
    protected:
        std::array<std::array<T, 2>, 2> val; ///< matrix elements

    public:
        static constexpr matrix2x2<T> IDENTITY {1, 0, 0, 1}; ///< Identity matrix

        /// \brief Default constructor
        ///
        /// Initializes matrix elements with zero
        ///
        constexpr matrix2x2(): matrix2x2(IDENTITY)
        {
        }


        /// \brief Initializes matrix with a00-a11 elements
        ///
        /// \param a00-a11: elements of matrix starting with zero
        ///
        constexpr matrix2x2(T a00, T a01, T a10, T a11) : val({a00, a01, a10, a11}) 
        {
        }


        /// \brief Copy constructor
        ///
        /// \param src: source 2x2 matrix
        ///
        constexpr matrix2x2(const matrix2x2<T>& src) : val(src.val)
        {
        };


        /// \brief construct matrix from two column vectors as follows:
        /// [p.x q.x
        ///  p.y q.y]
        ///
        /// \param p: 1st vector
        /// \param q: 2nd vector
        ///
        constexpr matrix2x2(const vector2<T>& p, const vector2<T>& q) : val({p.x, p.y}, {q.x, q.y})
        {
        }


        /// \brief Constructs matrix transpose
        ///
        /// \return 2x2 matrix
        ///
        constexpr matrix2x2<T> transposed() const 
        {
            return matrix2x2<T>(val[0][0], val[1][0], val[0][1], val[1][1]);
        }


        /// \brief Computes determinant of matrix
        ///
        /// \return determinant
        ///
        constexpr T determinant() const 
        {
            return val[0][0] * val[1][1] - val[0][1] * val[1][0];
        }  


        /// \brief Constructs classical adjoint
        ///
        /// \return 2x2 matrix
        ///
        /// Classical adjoint of matrix M is transpose of
        /// the matrix of cofactors of M.
        ///
        constexpr matrix2x2<T> adjoint() const
        {
            return matrix2x2<T>(val[1][1], -val[0][1],
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
        constexpr matrix2x2<T> inverse() const
        {
            const T d = determinant();
            return adjoint() / d;
        }


        /// \brief Gram-Shmidt orthonormalization
        ///
        /// Biased towards OX orthonormalization
        ///
        void orthonormalize()
        {
            vector2 ox(val[0][0], val[0][1]);
            vector2 oy(val[1][0], val[1][1]);

            ox.normalize();
            oy = oy - ox.dot(oy)*ox;
            oy.normalize();

            val[0][0] = ox.x;
            val[0][1] = ox.y;
            val[1][0] = oy.x;
            val[1][1] = oy.y;
        }


        constexpr matrix2x2<T> operator*(T scalar) const 
        {
            return matrix2x2(val[0][0] * scalar, val[0][1] * scalar,
                            val[1][0] * scalar, val[1][1] * scalar);
        }

        matrix2x2& operator*=(T scalar) 
        {
            val[0][0] *= scalar;
            val[0][1] *= scalar;
            val[1][0] *= scalar;
            val[1][1] *= scalar;
            return *this;
        }

        constexpr matrix2x2<T> operator/(T scalar) const 
        {
            return matrix2x2(val[0][0] / scalar, val[0][1] / scalar,
                            val[1][0] / scalar, val[1][1] / scalar);
        }

        matrix2x2<T>& operator/=(T scalar) 
        {
            val[0][0] /= scalar;
            val[0][1] /= scalar;
            val[1][0] /= scalar;
            val[1][1] /= scalar;
            return *this;
        }


        /// \brief Multiple by column vector
        ///
        /// \param v: 2D column vector
        /// \return 2D column vector
        ///
        constexpr vector2<T> operator*(const vector2<T>& v) const 
        {
            return vector2<T>(val[0][0] * v.x + val[1][0] * v.y,
                        val[0][1] * v.x + val[1][1] * v.y);
        }

        constexpr matrix2x2<T> operator*(const matrix2x2<T>& b) const 
        {
            return matrix2x2<T>(val[0][0] * b[0][0] + val[1][0] * b[0][1],
                            val[0][1] * b[0][0] + val[1][1] * b[0][1],
                            val[0][0] * b[1][0] + val[1][0] * b[1][1],
                            val[0][1] * b[1][0] + val[1][1] * b[1][1]);
        }

        constexpr matrix2x2<T>& operator*=(const matrix2x2<T>& b) {
            *this = *this * b;
            return *this;
        }

        matrix2x2<T>& operator = (const matrix2x2<T>& src)
        {
            val = src.val;
            return *this;
        }

        constexpr bool operator == (const matrix2x2& op2) const
        {
            return val == op2.val;
        }
        std::array<T, 2>& operator [] (size_t indx)
        {
            return val[indx];
        }

        constexpr std::array<T, 2>& operator [] (size_t indx) const
        {
            return val[indx];
        }
    };

    template<typename C>
    constexpr matrix2x2<C> operator * (C scalar, const matrix2x2<C>& m)
    {
        return m * scalar;
    }
} // namespace math

#endif // MATRIX2X2_H_INCLUDED
