#ifndef POLAR_H_INCLUDED
#define POLAR_H_INCLUDED

#include "math_predefs.h"
#include "angle.h"

namespace math
{
    class vector2;

    /// \brief Vector in polar coordinate system (2D).
    /// Axis corresponds to +X in Cartesian 2D system,
    class polar
    {
    public:
        real rho;///< radius
        real phi;///< angle (radians)

        /// \brief Default constructor
        ///
        /// Initializes rho and phi with zero (origin)
        ///
        polar() : rho(0.0f), phi(0.0f)
        {
        }

        /// \brief Initializes vector with given values
        ///
        /// \param rho real: radius
        /// \param phi real: angle (radians)
        ///
        polar(real rho, real phi) : rho(rho), phi(phi)
        {
        }

        /// \brief Convertion to 2D Cartesian coordinate system
        ///
        operator vector2 () const;

        /// \brief Canonize coordinates
        ///
        /// \return void
        ///
        /// Polar coordinates (phi, rho) are in canonical form if:
        /// rho >= 0,
        /// -PI < phi <= PI,
        /// at rho=0 phi=0.
        ///
        void canonize();
    };

} // namespace math

#endif // POLAR_H_INCLUDED
