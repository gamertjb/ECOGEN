#include "../src/GPR/GPRModel.h"
#include <iostream>

int main()
{
    GPRState state{};
    state.rho = 1.0;
    state.u = {1.0, 0.0, 0.0};
    state.chi = 0.0;
    state.j = {0.1, 0.0, 0.0};
    state.e = 2.5;
    state.p = 1.0;
    state.T = 1.0;

    double alpha2 = 1.0;
    double tau = 0.1;
    double rho0 = 1.0;
    double T0 = 1.0;

    GPRFlux flux = GPRModel::computeFlux(state, alpha2);
    GPRState source = GPRModel::computeSource(state, alpha2, tau, rho0, T0);

    std::array<double,3> n{1.0,0.0,0.0};
    auto lam = GPRModel::eigenvalues(state, n, 1.0, 0.5, 1.4, 1.0);

    std::cout << "Mass flux: " << flux.mass << "\n";
    std::cout << "Source chi: " << source.chi << "\n";
    std::cout << "First eigenvalue: " << lam[0] << "\n";
    return 0;
}

