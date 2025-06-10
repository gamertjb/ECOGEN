#include "../src/GPR/GPRSolver1D.h"
#include <iostream>

int main()
{
    const int cells = 200;
    const double length = 1.0;
    const double dx = length / cells;
    const double dt = 0.0005;
    double alpha2 = 1.0;
    double tau = 0.1;
    double rho0 = 1.0;
    double T0 = 1.0;
    double cv = 1.0;
    double cp = 1.4;

    GPRSolver1D solver(cells, dx, dt, alpha2, tau, rho0, T0, cv, cp);

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

    solver.setInitialCondition([&](double x){
        if (x < 0.5) return left; else return right;
    });

    solver.run(0.02); // integrate for short time

    const auto &cellsVec = solver.cells();
    for (int i=0;i<cells;i+=20) {
        std::cout << i*dx << " " << cellsVec[i].rho << " " << cellsVec[i].u[0] << " " << cellsVec[i].p << "\n";
    }
}
