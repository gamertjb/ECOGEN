#include "RelaxationGPR.h"

RelaxationGPR::RelaxationGPR(double alpha2, double tau,
                             double rho0, double T0)
    : m_alpha2(alpha2), m_tau(tau), m_rho0(rho0), m_T0(T0)
{}

RelaxationGPR::~RelaxationGPR(){}

void RelaxationGPR::relaxation(Cell* cell, const double& dt, Prim type)
{
    // Build a GPR state from mixture values
    Mixture* mix = cell->getMixture(type);
    GPRState s{};
    s.rho = mix->getDensity();
    s.u[0] = mix->getVelocity().getX();
    s.u[1] = mix->getVelocity().getY();
    s.u[2] = mix->getVelocity().getZ();
    s.e = mix->getTotalEnergy();
    s.p = mix->getPressure();
    s.T = mix->getTemperature();
    s.chi = 0.0;
    s.j = {0.0, 0.0, 0.0};

    // Two stage Runge-Kutta integration of the relaxation source terms
    GPRState k1 = GPRModel::computeSource(s, m_alpha2, m_tau, m_rho0, m_T0);

    GPRState mid = s;
    mid.chi += 0.5 * dt * k1.chi / s.rho;
    for (int k=0;k<3;++k)
        mid.j[k] += 0.5 * dt * k1.j[k] / s.rho;

    GPRState k2 = GPRModel::computeSource(mid, m_alpha2, m_tau, m_rho0, m_T0);

    s.chi += dt * k2.chi / s.rho;
    for (int k=0;k<3;++k)
        s.j[k] += dt * k2.j[k] / s.rho;

    // Store back updated pressure and temperature for the mixture
    mix->setPressure(s.p);
    mix->setTemperature(s.T);
    mix->setTotalEnergy(s.e);
    cell->fulfillState(type);
}
