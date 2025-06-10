#ifndef GPRSOLVER1D_H
#define GPRSOLVER1D_H

#include "GPRModel.h"
#include <vector>
#include <functional>

//! \class GPRSolver1D
//! \brief Simple 1D finite volume solver for the GPR system
class GPRSolver1D {
  public:
    //! \brief Constructor
    GPRSolver1D(int cells, double dx, double dt, double alpha2,
                double tau, double rho0, double T0,
                double cv, double cp);

    //! \brief Set initial condition using a function of position
    void setInitialCondition(const std::function<GPRState(double)> &ic);

    //! \brief Perform one time step
    void step();

    //! \brief Run simulation until timeEnd
    void run(double timeEnd);

    //! \brief Access cell states
    const std::vector<GPRState> &cells() const { return m_cells; }

  private:
    struct Conserv {
        double rho;
        std::array<double,3> mom;
        double chi;
        std::array<double,3> J;
        double rhoE;
    };

    Conserv toConservative(const GPRState &s) const;
    GPRState fromConservative(const Conserv &c) const;

    GPRFlux hllFlux(const GPRState &L, const GPRState &R) const;

    int m_cellsCount; //!< number of cells
    double m_dx; //!< cell size
    double m_dt; //!< time step
    double m_alpha2; //!< alpha squared
    double m_tau; //!< relaxation time
    double m_rho0; //!< reference density
    double m_T0; //!< reference temperature
    double m_cv; //!< specific heat capacity cv
    double m_cp; //!< specific heat capacity cp

    std::vector<GPRState> m_cells; //!< solution values
    std::vector<GPRFlux> m_fluxes; //!< interface fluxes
};

#endif
