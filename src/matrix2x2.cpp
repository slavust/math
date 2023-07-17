#include "matrix2x2.h"

namespace math
{
    const matrix2x2 matrix2x2::IDENTITY(1.0f, 0.0f,
                                        0.0f, 1.0f);

    void matrix2x2::orthonormalize()
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

} // namespace math

