#include "GPRModel.h"
#include <cmath>

GPRFlux GPRModel::computeFlux(const GPRState& s, double alpha2)
{
    GPRFlux f{};
    // mass flux (x-direction)
    f.mass = s.rho * s.u[0];

    // momentum flux tensor contracted in x-direction
    for (int i = 0; i < 3; ++i) {
        f.momentum[i] = s.rho * s.u[i] * s.u[0];
    }
    f.momentum[0] += s.p;

    // entropy (chi) flux
    f.entropy = s.rho * s.u[0] * s.chi + alpha2 * s.j[0];

    // thermal impulse flux
    for (int i = 0; i < 3; ++i) {
        f.impulse[i] = s.rho * s.j[i] * s.u[0];
    }
    f.impulse[0] += s.T;

    // energy flux
    f.energy = (s.rho * s.e + s.p) * s.u[0] + alpha2 * s.T * s.j[0];

    return f;
}

GPRState GPRModel::computeSource(const GPRState& s, double alpha2,
                                 double tau, double rho0, double T0)
{
    GPRState src{};
    double theta = tau * alpha2 * s.rho / rho0 * T0 / s.T;
    double j2 = s.j[0]*s.j[0] + s.j[1]*s.j[1] + s.j[2]*s.j[2];
    src.chi = - s.rho * alpha2 * j2 / (theta * s.T);
    for (int i = 0; i < 3; ++i) {
        src.j[i] = - alpha2 * s.rho * s.j[i] / theta;
    }
    // other components are zero
    return src;
}

std::array<double,8> GPRModel::eigenvalues(const GPRState& s,
                                           const std::array<double,3>& n,
                                           double cs, double ch,
                                           double cp, double cv)
{
    double un = s.u[0]*n[0] + s.u[1]*n[1] + s.u[2]*n[2];
    double disc = std::pow(cs*cs + ch*ch, 2) - 4.0 * cs*cs * ch*ch * cv / cp;
    disc = std::max(disc, 0.0);
    double root = std::sqrt(disc);
    double lamPlus = std::sqrt(0.5 * (cs*cs + ch*ch + root));
    double lamMinus = std::sqrt(0.5 * (cs*cs + ch*ch - root));

    std::array<double,8> lam{};
    lam[0] = un - lamPlus;
    lam[1] = un - lamMinus;
    for (int i=2;i<6;++i) lam[i] = un;
    lam[6] = un + lamMinus;
    lam[7] = un + lamPlus;
    return lam;
}

