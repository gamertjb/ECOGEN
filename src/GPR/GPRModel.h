#ifndef GPRMODEL_H
#define GPRMODEL_H

#include "GPRState.h"
#include "GPRFlux.h"
#include <array>

//! \class GPRModel
//! \brief Utility functions for the GPR governing equations
class GPRModel {
  public:
    //! \brief Compute physical flux in x-direction
    static GPRFlux computeFlux(const GPRState& s, double alpha2);

    //! \brief Compute algebraic source terms
    static GPRState computeSource(const GPRState& s, double alpha2,
                                  double tau, double rho0, double T0);

    //! \brief Compute eigenvalues in direction n
    static std::array<double,8> eigenvalues(const GPRState& s,
                                            const std::array<double,3>& n,
                                            double cs, double ch,
                                            double cp, double cv);
};

#endif
