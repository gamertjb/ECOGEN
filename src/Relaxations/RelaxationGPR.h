#ifndef RELAXATIONGPR_H
#define RELAXATIONGPR_H

#include "Relaxation.h"
#include "../GPR/GPRModel.h"

//! \class     RelaxationGPR
//! \brief     Robust 3D relaxation using the GPR model
class RelaxationGPR : public Relaxation
{
public:
    //! \brief Constructor with default GPR parameters
    RelaxationGPR(double alpha2 = 1.0, double tau = 1e-3,
                  double rho0 = 1.0, double T0 = 1.0);
    virtual ~RelaxationGPR();

    //! \brief Apply a single GPR relaxation step
    virtual void relaxation(Cell* cell, const double& dt,
                            Prim type = vecPhases);

    //! \brief Return GPR relaxation type
    virtual int getType() const { return GPR; }

private:
    double m_alpha2; //!< GPR alpha^2 parameter
    double m_tau;    //!< relaxation time
    double m_rho0;   //!< reference density
    double m_T0;     //!< reference temperature
};

#endif // RELAXATIONGPR_H
