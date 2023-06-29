#ifndef CYLINDRICAL_H_INCLUDED
#define CYLINDRICAL_H_INCLUDED

#include "math_predefs.h"
#include "angle.h"

namespace math
{
    class vector3;

    /// \brief Vector in cylindrical coordinate system (3D).
    /// Axis corresponds to +X in Cartesian 3D system.
    class cylindrical
    {
    public:
        real rho;///< radius
        real phi;///< azimuth (radians)
        real z;///< height

        /// \brief Default constructor.
        ///
        /// Initializes rho, phi and z with zero
        ///
        cylindrical() : rho(0.0f), phi(0.0f), z(0.0f)
        {
        }

        /// \brief initializes rho, phi, z with given values
        ///
        /// \param rho real: radius
        /// \param phi real: angle (radians)
        /// \param z real: height
        ///
        cylindrical(real rho, real phi, real z) : rho(rho), phi(phi), z(z)
        {
        }

        /// \brief Convertion to 3D Cartesian coordinate system
        ///
        operator vector3 () const;

        /// \brief Canonize coordinates
        ///
        /// \return void
        ///
        /// Canonizes rho, phi as for polar coordinate system:
        /// rho >= 0,
        /// -PI < phi <= PI,
        /// at rho=0 phi=0.
        ///
        void canonize();
    };
} // namespace math

#endif // CYLINDRICAL_H_INCLUDED
