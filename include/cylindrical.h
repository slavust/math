#ifndef CYLINDRICAL_H_INCLUDED
#define CYLINDRICAL_H_INCLUDED

#include "math_predefs.h"
#include "angle.h"

namespace math
{
    /// \brief Vector in cylindrical coordinate system (3D).
    /// Polar axis corresponds to +X in Cartesian 3D system.
    template <typename T>
    class cylindrical
    {
    public:
        T rho;///< radius
        T phi;///< azimuth (radians)
        T z;///< height

        /// \brief Default constructor.
        ///
        /// Initializes rho, phi and z with zero
        ///
        constexpr cylindrical() : rho(0), phi(0), z(0)
        {
        }

        /// \brief initializes rho, phi, z with given values
        ///
        /// \param rho real: radius
        /// \param phi real: angle (radians)
        /// \param z real: height
        ///
        constexpr cylindrical(T rho, T phi, T z) : rho(rho), phi(phi), z(z)
        {
        }

        /// \brief Canonize coordinates
        ///
        /// \return void
        ///
        /// Canonizes rho, phi as for polar coordinate system:
        /// rho >= 0,
        /// -PI < phi <= PI,
        /// at rho=0 phi=0.
        ///
        void canonize()
        {
            if(rho == 0.0f) // at the origin
            {
                // at rho==0.0f phi = 0.0f
                phi = 0.0f;
            }
            else
            {
                if(rho < 0.0f) // negative distance
                {
                    //make positive distance
                    rho = -rho;
                    phi += pi<T>();
                }
                if(abs(phi) > pi<T>()) // phi out of range
                {
                    //make phi within range [-PI, PI)
                    phi -= floor((phi + pi<T>()) / (2 * pi<T>())) * 2 * pi<T>();
                }
            }
        }
    };
} // namespace math

#endif // CYLINDRICAL_H_INCLUDED
