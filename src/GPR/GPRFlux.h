#ifndef GPRFLUX_H
#define GPRFLUX_H

#include <array>

//! \brief Physical flux vector for the GPR model
struct GPRFlux {
    double mass;                       //!< Mass flux
    std::array<double,3> momentum;     //!< Momentum flux
    double entropy;                    //!< Entropy flux
    std::array<double,3> impulse;      //!< Thermal impulse flux
    double energy;                     //!< Energy flux
};

#endif
