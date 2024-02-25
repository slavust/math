#ifndef POLAR_H_INCLUDED
#define POLAR_H_INCLUDED

#include "math_predefs.h"
#include "angle.h"

namespace math
{
    /// \brief Vector in polar coordinate system (2D).
    /// Axis corresponds to +X in Cartesian 2D system,
    template<typename T>
    class polar
    {
    public:
        T rho;///< radius
        T phi;///< angle (radians)

        /// \brief Default constructor
        ///
        /// Initializes rho and phi with zero (origin)
        ///
        constexpr polar() : rho(0), phi(0)
        {
        }

        /// \brief Initializes vector with given values
        ///
        /// \param rho real: radius
        /// \param phi real: angle (radians)
        ///
        constexpr polar(T rho, T phi) : rho(rho), phi(phi)
        {
        }

        /// \brief Canonize coordinates
        ///
        /// \return void
        ///
        /// Polar coordinates (phi, rho) are in canonical form if:
        /// rho >= 0,
        /// -PI < phi <= PI,
        /// at rho=0 phi=0.
        ///
        
        void canonize()
        {
            if(rho == 0) // at the origin
            {
                // at the origin phi = 0
                phi = 0;
            }
            else
            {
                if(rho < 0) // negative distance
                {
                    //make positive distance
                    rho = -rho;
                    phi += pi<T>();
                }
                if(abs(phi) > pi<T>()) // phi out of range
                {
                    //make phi within range (-PI, PI]
                    phi -= floor((phi + pi<T>()) / (pi<T>() * 2)) * pi<T>() * 2;
                }
            }
        }
    };
} // namespace math

#endif // POLAR_H_INCLUDED
