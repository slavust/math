#ifndef SPHERICAL_H_INCLUDED
#define SPHERICAL_H_INCLUDED

#include "math_predefs.h"
#include "angle.h"

namespace math
{
    /// \brief Vector in spherical coordinate system (3D).
    /// Horizontal axis corresponds to +X and
    /// vertical axis corresponds to +Z in Cartezian 3D system.
template<typename T>
class spherical
{
public:
    T r; ///< radius
    T theta; ///< inclination (radians)
    T phi; ///< azimuth (radians)

    /// \brief Default constructor
    ///
    /// Initializes r, theta and phi with zero
    ///
    constexpr spherical() : r(0), theta(0), phi(0)
    {
    }

    /// \brief Initializes r, theta and phi with given values
    ///
    /// \param r T: radius
    /// \param theta T: inclination (radians)
    /// \param phi T: azimuth (radians)
    ///
    constexpr spherical(T r, T theta, T phi) : r(r), theta(theta), phi(phi)
    {
    }

    /// \brief Coordinate negation
    ///
    constexpr spherical operator - () const { return spherical(-r, theta, phi); }

    /// \brief Canonize coordinates
    /// Canonical spherical coordinates are defined as follows:
    /// r >= 0,
    /// 0 <= theta <= PI,
    /// -PI <= phi <= PI,
    /// at r=0, theta=0,
    /// at theta=0 and theta=PI, phi=0.
    ///
    void canonize()
    {
        if(r == 0) // at the origin
        {
            // at r == 0, theta = phi = 0
            theta = 0;
            phi = 0;
        }
        else
        {
            if(r < 0) // negative radius
            {
                // make positive radius
                r = -r;
                phi += pi<T>();
            }

            // make theta within range [0, PI]
            if(theta < 0)
            {
                // make positive theta
                theta = -theta;
                phi += pi<T>();
            }
            // put theta in first period
            theta -= floor(theta / (pi<T>() * 2)) * pi<T>() * 2;
            if(theta > pi<T>())
            {
                // subtract from TWO_PI theta without TWO_PI portions. 0<=theta<=PI
                theta = pi<T>() * 2 - theta;
                phi += pi<T>();
            }

            // at theta = 0 and theta = PI, phi = 0
            if(theta < 0 || abs(theta - pi<T>()) < 0) phi = 0;
            else if(abs(phi) > pi<T>())
            {
                //make phi within range [-PI, PI)
                phi -= floor((phi + pi<T>()) / (pi<T>() * 2)) * pi<T>() * 2;
            }
        }
    }
};
} // namespace math

#endif // SPHERICAL_H_INCLUDED
