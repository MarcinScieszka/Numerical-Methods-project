#include <iostream>
#include "calerf.h"
#include "cmath"

double const r = 1.0;
double const a = 10.0;

#define D 1

#define T_0 0
#define T_MAX 2

//x belongs to the interval from r to r+a
#define X_MIN r
#define X_MAX (r + a)

//lambda = D * (dt/h*h);
double lambda_kmb = 0.4;
double lambda_ml = 1.0;

// initial condition: U(x,0) = 1
// first boundary condition: U(r,t) = 0
// second boundary condition: U(r+a, t) = 1.0 - (r/(r+a))*calerf::ERFCL(a/(2.0*sqrt(D*t)));

double analytical_solution(double x, double t);

int main() {

    double result;

    result = analytical_solution(2.0, 1.0);

    std::cout << result << '\n';

    return 0;
}

double analytical_solution(double x, double t)
{
    return 1.0 - (r/x)*calerf::ERFCL((x-r)/(2.0*sqrt(D*t)));
}
