#include "../src/GPR/GPRModel.h"
#include <iostream>
#include <array>

struct Conserv {
    double rho;
    std::array<double,3> mom;
    double chi;
    std::array<double,3> J;
    double rhoE;
};

static Conserv toConservative(const GPRState& s)
{
    Conserv c{};
    c.rho = s.rho;
    for(int i=0;i<3;++i) c.mom[i] = s.rho*s.u[i];
    c.chi = s.rho*s.chi;
    for(int i=0;i<3;++i) c.J[i] = s.rho*s.j[i];
    c.rhoE = s.rho*s.e;
    return c;
}

int main()
{
    GPRState left{};
    left.rho = 0.6635;
    left.u = {0.0,0.0,0.0};
    left.chi = 0.0;
    left.j = {0.0,0.0,0.0};
    left.e = 2.5;
    left.p = 0.0313;
    left.T = 0.9;

    GPRState right{};
    right.rho = 0.0178;
    right.u = {0.0,0.0,0.0};
    right.chi = 0.0;
    right.j = {0.0,0.0,0.0};
    right.e = 2.0;
    right.p = 0.0126;
    right.T = 0.8;

    double alpha2 = 1.0;
    auto fL = GPRModel::computeFlux(left, alpha2);
    auto fR = GPRModel::computeFlux(right, alpha2);

    Conserv UL = toConservative(left);
    Conserv UR = toConservative(right);

    double sL = -1.0; // left-going wave speed guess
    double sR = 1.0;  // right-going wave speed guess

    Conserv flux{};
    flux.rho = (sR*fL.mass - sL*fR.mass + sL*sR*(UR.rho-UL.rho))/(sR - sL);
    for(int i=0;i<3;++i)
        flux.mom[i] = (sR*fL.momentum[i] - sL*fR.momentum[i] + sL*sR*(UR.mom[i]-UL.mom[i]))/(sR - sL);
    flux.chi = (sR*fL.entropy - sL*fR.entropy + sL*sR*(UR.chi-UL.chi))/(sR - sL);
    for(int i=0;i<3;++i)
        flux.J[i] = (sR*fL.impulse[i] - sL*fR.impulse[i] + sL*sR*(UR.J[i]-UL.J[i]))/(sR - sL);
    flux.rhoE = (sR*fL.energy - sL*fR.energy + sL*sR*(UR.rhoE-UL.rhoE))/(sR - sL);

    std::cout << "HLL mass flux: " << flux.rho << "\n";
    std::cout << "HLL momentum flux x: " << flux.mom[0] << "\n";
    std::cout << "HLL energy flux: " << flux.rhoE << "\n";
}

