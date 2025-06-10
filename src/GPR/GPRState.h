#ifndef GPRSTATE_H
#define GPRSTATE_H

#include <array>

//! \brief State vector for the GPR model
struct GPRState {
    double rho;                   //!< Density
    std::array<double,3> u;       //!< Velocity vector
    double chi;                   //!< Specific entropy
    std::array<double,3> j;       //!< Thermal impulse per unit mass
    double e;                     //!< Total energy per unit mass
    double p;                     //!< Pressure
    double T;                     //!< Temperature
};

#endif
