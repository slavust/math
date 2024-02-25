#pragma once

#include "euler.h"
#include "matrix3x3.h"
#include "matrix4x4.h"
#include "matrix_transforms.h"

#include "vector2.h"
#include "polar.h"

#include "vector3.h"
#include "cylindrical.h"
#include "spherical.h"

namespace math
{
    /// Convert Euler angles to 3x3 rotation matrix
    template<typename T>
    matrix3x3<T> to_rotation_matrix3x3(const euler<T>& e)
    {
        return rotate3Dz<T>(e.roll) * rotate3Dx<T>(e.pitch) * rotate3Dy<T>(e.yaw);
    }

    /// Convert euler angles to 4x4 rotation matrix
    template<typename T>
    matrix4x4<T> to_rotation_matrix4x4(const euler<T>& e)
    {
        return rotate3Dhz<T>(e.roll) * rotate3Dhx<T>(e.pitch) * rotate3Dhy<T>(e.yaw);
    }

    /// Convert polar coordinates to Cartesian 2d vector
    template<typename T>
    vector2<T> to_vector2(const polar<T>& p)
    {
        return vector2<T>(p.rho * cos(p.phi), p.rho * sin(p.phi));
    }

    /// Convert Cartesian 2d vector to polar coordinates
    template<typename T>
    polar<T> to_polar(const vector2<T>& v)
    {
        if(v == vector2<T>::ZERO) 
            return polar(0.0f, 0.0f);
        return polar(v.magnitude(), atan2(v.y, v.x));
    }

    /// Convert cylindrical coordinates to Cartesian 3d vector
    template<typename T>
    vector3<T> to_vector3(const cylindrical<T>& c)
    {
        return vector3<T>(c.rho*cos(c.phi), c.rho*sin(c.phi), c.z);
    }

    /// Convert Cartesian 3d vector to cylindrical coordinates
    template<typename T>
    cylindrical<T> from_vector3(const vector3<T>& v)
    {
        if(v.x == 0 && v.y == 0) 
            return cylindrical(0.0f, 0.0f, v.z);
        return cylindrical(sqrt(v.x*v.x + v.y*v.y), atan2(v.y, v.x), v.z);
    }

    /// Convert spherical coordinates to Cartesian 3d vector
    template<typename T>
    vector3<T> to_vector3(const spherical<T>& s)
    {
        const T sin_theta = sin(s.theta);
        return vector3<T>(
            s.r * cos(s.phi) * sin_theta, 
            s.r * sin(s.phi) * sin_theta, 
            s.r * cos(s.theta));
    }

    /// convert Cartesian 3d vector to spherical coordinates
    template<typename T>
    spherical<T> to_spherical(const vector3<T>& v)
    {
        const T r = v.magnitude();
        const T theta = acos(v.z / r);

        T phi;
        if(v.x == 0 && v.y == 0)
            phi = 0;
        else
            phi = atan2(v.y, v.x);

        return spherical<T>(r, theta, phi);
    }
}