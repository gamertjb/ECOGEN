#include "GPRSolver1D.h"
#include <cmath>

GPRSolver1D::GPRSolver1D(int cells, double dx, double dt, double alpha2,
                         double tau, double rho0, double T0,
                         double cv, double cp)
    : m_cellsCount(cells), m_dx(dx), m_dt(dt), m_alpha2(alpha2), m_tau(tau),
      m_rho0(rho0), m_T0(T0), m_cv(cv), m_cp(cp)
{
    m_cells.resize(cells);
    m_fluxes.resize(cells + 1);
}

void GPRSolver1D::setInitialCondition(const std::function<GPRState(double)> &ic)
{
    for (int i = 0; i < m_cellsCount; ++i) {
        double x = (i + 0.5) * m_dx;
        m_cells[i] = ic(x);
    }
}

GPRSolver1D::Conserv GPRSolver1D::toConservative(const GPRState &s) const
{
    Conserv c{};
    c.rho = s.rho;
    for (int i = 0; i < 3; ++i) {
        c.mom[i] = s.rho * s.u[i];
        c.J[i] = s.rho * s.j[i];
    }
    c.chi = s.rho * s.chi;
    c.rhoE = s.rho * s.e;
    return c;
}

GPRState GPRSolver1D::fromConservative(const Conserv &c) const
{
    GPRState s{};
    s.rho = c.rho;
    for (int i = 0; i < 3; ++i) {
        s.u[i] = c.mom[i] / c.rho;
        s.j[i] = c.J[i] / c.rho;
    }
    s.chi = c.chi / c.rho;
    s.e = c.rhoE / c.rho;
    double j2 = s.j[0]*s.j[0] + s.j[1]*s.j[1] + s.j[2]*s.j[2];
    double kinetic = 0.5*(s.u[0]*s.u[0] + s.u[1]*s.u[1] + s.u[2]*s.u[2]);
    double eps = s.e - kinetic - 0.5 * m_alpha2 * j2;
    double gamma = m_cp / m_cv;
    s.T = eps / m_cv;
    s.p = (gamma - 1.0) * s.rho * eps;
    return s;
}

GPRFlux GPRSolver1D::hllFlux(const GPRState &L, const GPRState &R) const
{
    GPRFlux fL = GPRModel::computeFlux(L, m_alpha2);
    GPRFlux fR = GPRModel::computeFlux(R, m_alpha2);
    double gamma = m_cp / m_cv;
    double aL = std::sqrt(gamma * L.p / L.rho);
    double aR = std::sqrt(gamma * R.p / R.rho);
    double sL = std::min(L.u[0] - aL, R.u[0] - aR);
    double sR = std::max(L.u[0] + aL, R.u[0] + aR);

    if (sL >= 0.0)
        return fL;
    if (sR <= 0.0)
        return fR;

    Conserv UL = toConservative(L);
    Conserv UR = toConservative(R);
    GPRFlux flux{};
    flux.mass = (sR*fL.mass - sL*fR.mass + sL*sR*(UR.rho-UL.rho))/(sR - sL);
    for (int i=0;i<3;++i)
        flux.momentum[i] = (sR*fL.momentum[i] - sL*fR.momentum[i] +
                            sL*sR*(UR.mom[i]-UL.mom[i]))/(sR - sL);
    flux.entropy = (sR*fL.entropy - sL*fR.entropy + sL*sR*(UR.chi-UL.chi))/(sR - sL);
    for (int i=0;i<3;++i)
        flux.impulse[i] = (sR*fL.impulse[i] - sL*fR.impulse[i] +
                           sL*sR*(UR.J[i]-UL.J[i]))/(sR - sL);
    flux.energy = (sR*fL.energy - sL*fR.energy + sL*sR*(UR.rhoE-UL.rhoE))/(sR - sL);
    return flux;
}

void GPRSolver1D::step()
{
    // compute interface fluxes
    for (int i = 0; i <= m_cellsCount; ++i) {
        const GPRState &L = (i == 0) ? m_cells[0] : m_cells[i-1];
        const GPRState &R = (i == m_cellsCount) ? m_cells[m_cellsCount-1] : m_cells[i];
        m_fluxes[i] = hllFlux(L, R);
    }

    std::vector<GPRState> next = m_cells;
    for (int i = 0; i < m_cellsCount; ++i) {
        Conserv U = toConservative(m_cells[i]);
        Conserv Unew = U;
        Unew.rho -= m_dt/m_dx * (m_fluxes[i+1].mass - m_fluxes[i].mass);
        for (int k=0;k<3;++k)
            Unew.mom[k] -= m_dt/m_dx * (m_fluxes[i+1].momentum[k] - m_fluxes[i].momentum[k]);
        Unew.chi -= m_dt/m_dx * (m_fluxes[i+1].entropy - m_fluxes[i].entropy);
        for (int k=0;k<3;++k)
            Unew.J[k] -= m_dt/m_dx * (m_fluxes[i+1].impulse[k] - m_fluxes[i].impulse[k]);
        Unew.rhoE -= m_dt/m_dx * (m_fluxes[i+1].energy - m_fluxes[i].energy);
        next[i] = fromConservative(Unew);
        // add source terms
        GPRState src = GPRModel::computeSource(next[i], m_alpha2, m_tau, m_rho0, m_T0);
        next[i].chi += m_dt * src.chi / next[i].rho;
        for (int k=0;k<3;++k)
            next[i].j[k] += m_dt * src.j[k] / next[i].rho;
        double j2 = next[i].j[0]*next[i].j[0] + next[i].j[1]*next[i].j[1] + next[i].j[2]*next[i].j[2];
        double kinetic = 0.5*(next[i].u[0]*next[i].u[0]+next[i].u[1]*next[i].u[1]+next[i].u[2]*next[i].u[2]);
        double eps = next[i].e - kinetic - 0.5*m_alpha2*j2;
        next[i].T = eps / m_cv;
        double gamma = m_cp/m_cv;
        next[i].p = (gamma - 1.0)*next[i].rho*eps;
    }
    m_cells.swap(next);
}

void GPRSolver1D::run(double timeEnd)
{
    int steps = static_cast<int>(std::ceil(timeEnd / m_dt));
    for (int s = 0; s < steps; ++s) {
        step();
    }
}

