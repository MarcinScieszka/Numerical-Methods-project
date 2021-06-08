#include <iostream>
#include "calerf.h"
#include "cmath"

/*
 *  Application of differential methods to an approximate solution of the diffusion equation with initial condition and boundary conditions
 *
 *  Algorithms used:
 *  - FTCS Explicit Method
 *  - Laasonen method using LU decomposition
 *  - Laasonen method using Thomas algorithm
 *
 *  Initial condition: U(x,0) = 1
 *  First boundary condition: U(r,t) = 0
 *  Second boundary condition: U(r+a, t) = 1.0 - (r/(r+a))*calerf::ERFCL(a/(2.0*sqrt(D*t)));
 *
 * */

double const r = 1.0; // radius of the sphere
double const a = 10.0;
double const D = 1.0;// diffusion coefficient
double const T_0 = 0.0; // starting moment
double const T_MAX = 2.0;

//x - spatial coordinate
double const X_MIN = r;
double const X_MAX = (r + a);

//lambda = D * (dt/h*h);
const double lambda_kmb = 0.4; // lambda value used for direct method, NOTE: method stable only when lambda < 0.5
const double lambda_ml = 1.0; // lambda value used for Laasonen method

double analytical_solution(double x, double t);

int main() {
    double dt; // delta t
    double h; // step h on the spatial grid
    double tk; // current time on the time grid
    double xi; // current point on the spatial grid
    double result; // approximate solution at node xi at time level tk

    dt = 5e-3;
    h = sqrt((D*dt)/lambda_kmb);
    tk = T_MAX;

    for (xi=X_MIN; xi<X_MAX; xi+=h)
    {
        result = analytical_solution(xi, tk);
        std::cout << result << '\n';
    }

    return 0;
}

double analytical_solution(double x, double t)
{
    return 1.0 - (r/x)*calerf::ERFCL((x-r)/(2.0*sqrt(D*t)));
}
