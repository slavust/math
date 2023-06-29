#ifndef MATH_PREDEFS_H_INCLUDED
#define MATH_PREDEFS_H_INCLUDED

#include <math.h>
#include <string.h>
#include <limits>
#include "assert.h"

#define MATH_EXCEPTION(ExceptionType) throw MathException(ExceptionType, __FILE__, __LINE__);
//#define MATH_CHECK_INDEX_BOUNDS
#define MATH_CHECK_DIVISION

namespace math
{
    typedef float real;

    static const real EPS = 1e-04f; //std::numeric_limits<real>::epsilon();

    inline real sign(real s)
    {
        return s < 0 ? -1.0f : 1.0f;
    }


    enum ExceptionType
    {
        ET_INDEX_OUT_OF_BOUNDS,
        ET_DIVIDE_BY_ZERO,
        ET_NON_INVERTIBLE_MATRIX
    };

    class MathException
    {
    protected:
        ExceptionType mType;
        const char* mFile;
        int mLine;

    public:
        MathException(ExceptionType type, const char* file, int line) :
            mType(type), mFile(file), mLine(line)
        {
        }

        ExceptionType getType() const
        {
            return mType;
        }

        const char* getFile() const
        {
            return mFile;
        }

        int getLine() const
        {
            return mLine;
        }
    };

#ifdef MATH_CHECK_INDEX_BOUNDS
    template<typename T, size_t _size> class ArrayHolder
    {
    protected:
        T* _ptr;

    public:
        ArrayHolder(T* arr) : _ptr(arr)
        {
        }

        T& operator [] (size_t indx)
        {
            if(indx >= _size) MATH_EXCEPTION(ET_INDEX_OUT_OF_BOUNDS);
            return _ptr[indx];
        }

        T* ptr()
        {
            return _ptr;
        }
    };
#endif // MATH_CHECK_INDEX_BOUNDS

} // namespace math

#endif // MATH_PREDEFS_H_INCLUDED
