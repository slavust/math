#ifndef SPHERICAL_H_INCLUDED
#define SPHERICAL_H_INCLUDED

#include "math_predefs.h"
#include "angle.h"

namespace math
{
    class vector3;

    /// \brief Vector in spherical coordinate system (3D).
    /// Horizontal axis corresponds to +X and
    /// vertical axis corresponds to +Z in Cartezian 3D system.
    class spherical
    {
    public:
        real r; ///< radius
        real theta; ///< inclination (radians)
        real phi; ///< azimuth (radians)

        /// \brief Default constructor
        ///
        /// Initializes r, theta and phi with zero
        ///
        spherical() : r(0.0f), theta(0.0f), phi(0.0f)
        {
        }

        /// \brief Initializes r, theta and phi with given values
        ///
        /// \param r real: radius
        /// \param theta real: inclination (radians)
        /// \param phi real: azimuth (radians)
        ///
        spherical(real r, real theta, real phi) : r(r), theta(theta), phi(phi)
        {
        }

        /// \brief Conversion to 3D Cartesian coordinate system
        ///
        operator vector3 () const;

        /// \brief Coordinate negation
        ///
        spherical operator - () const { return spherical(-r, theta, phi); }

        /// \brief Canonize coordinates
        /// Canonical spherical coordinates are defined as follows:
        /// r >= 0,
        /// 0 <= theta <= PI,
        /// -PI <= phi <= PI,
        /// at r=0, theta=0,
        /// at theta=0 and theta=PI, phi=0.
        ///
        void canonize();
    };
} // namespace math

#endif // SPHERICAL_H_INCLUDED
