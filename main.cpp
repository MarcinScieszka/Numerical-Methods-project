#include <iostream>
#include <iomanip>
#include "calerf.h"
#include "cmath"

/**
 *  Application of differential methods to an approximate solution of the diffusion equation with initial condition and boundary conditions.
 *
 *  Algorithms used:
 *  - FTCS Explicit Method
 *  - Laasonen method using LU decomposition
 *  - Laasonen method using Thomas algorithm
 *
 *  Initial condition: U(x,0) = 1
 *  First boundary condition: U(r,t) = 0
 *  Second boundary condition: U(r+a, t) = 1.0 - (r/(r+a))*calerf::ERFCL(a/(2.0*sqrt(D*t)));
 * */

double const D = 1.0; // diffusion coefficient
double const T_0 = 0.0; // starting moment
double const T_MAX = 2.0; // max time
double const r = 1.0; // radius of the sphere
double const a = 10.0; // domain length
double const X_MIN = r; // x - spatial coordinate
double const X_MAX = (r + a);

//lambda = D*(dt/h*h) => dt=
const double lambda_dm = 0.4; // lambda value used for direct methods, NOTE: beware of numerical stability limitations - method stable only when lambda < 0.5
const double lambda_im = 1.0; // lambda value used for for indirect methods

//double analytical_solution(double x, double t);

double **allocate_matrix(int rows, int cols);
void print_matrix(double *const *matrix, int rows, int cols);
void free_matrix(double *const *matrix, int rows);

void print_array(double *array, int size);
void fill_array(double *array, int size, double value);

void fill_matrix(double **matrix, int rows, int cols, double value);

int main()
{
    double dt; // delta t
    double h; // step h on the spatial grid, delta x
    double tk; // current time on the time grid
    double xi; // current point on the spatial grid

//    ----- FTCS_explicit_method ----
    int nodes_x = 15; // nr of nodes on the time grid
    int nodes_t = 4; // nr of nodes on the spacial grid

    auto* ic_array = new double[nodes_x];  // initial condition array
    fill_array(ic_array, nodes_x, 1.0);
    print_array(ic_array, nodes_x);

    double **u = allocate_matrix(nodes_t, nodes_x); // matrix of approximate values
    fill_matrix(u, nodes_t, nodes_x, 0.0);
    print_matrix(u, nodes_t, nodes_x);
    int i, k;
    for (i=0; i<nodes_x; i++) u[0][i] = ic_array[i];
    print_matrix(u, nodes_t, nodes_x);

    for (k=0; k<nodes_t-1; k++)
    {
        for (i=1; i<nodes_x-1; i++)
        {
            u[k+1][i] = lambda_dm*u[k][i-1] + (1-lambda_dm)*u[k][i] + lambda_dm*u[k][i+1]; // approximate solution at node x_i at time level t_k+1
        }
    }
    print_matrix(u, nodes_t, nodes_x);


//    h = a/(nodes_x-1);
//    dt = T_MAX/(nodes_t-1);


//    double H_MIN = 1e-2;
//    int counter = 0;
//    while (h>H_MIN)
//    {
//        std::cout << "\nh: ";
//        std::cout << h << ", ";
//        std::cout << "\ndt: ";
//        std::cout << dt << ", ";
//        nodes_x *= 2;
//        nodes_t *= 2;
//        printf(" - nodes: %d \n", nodes_x);
//
////        u[i][k] = (h, tk);
//
//        counter++;
//        h = a/(nodes_x-1);
//        dt = T_MAX/(nodes_t-1);
//    }
//    std::cout << counter << "\n";


//    for (xi=X_MIN; xi<X_MAX; xi+=h)
//    {
//        continue;
//    }

    /**
     * TODO: dependence of the maximum absolute value of the error observed for T_MAX as a function of the spatial step h
     **/

//    tk = T_MAX - dt;

    delete[] ic_array;
    free_matrix(u, nodes_t);

    return 0;
}

double analytical_solution(double xi, double tk)
{
    return 1.0 - (r/xi)*double(calerf::ERFCL((xi-r)/(2.0*sqrt(D*tk))));
}

void fill_array(double *array, int size, double value)
{
    for (int i=0; i<size; i++)
    {
        array[i] = value;
    }
}

void print_array(double *array, int size)
{
    std::cout << "[ ";
    for(int i=0; i<size; i++)
    {
        std::cout << array[i] << " ";
    }
    std::cout << "] \n";
}

double **allocate_matrix(int rows, int cols)
{
    auto** matrix = new double*[rows];
    for (int i=0; i<rows; ++i)
        matrix[i] = new double[cols];
    return matrix;
}

void fill_matrix(double **matrix, int rows, int cols, double value)
{
    for (int i=0; i<rows; i++)
    {
        for (int j=0; j<cols; j++)
        {
            matrix[i][j] = value;
        }
    }
}

void print_matrix(double *const *matrix, const int rows, const int cols)
{
    for (int i=rows-1; i>=0; i--)
    {
        for (int j=0; j<cols; j++)
        {
            std::cout << matrix[i][j] << std::setw(7);
        }
        std::cout  << "\n";
    }
    std::cout  << "\n";
}

void free_matrix(double *const *matrix, int rows)
{
    for (int i=0; i<rows; ++i)
    {
        delete [] matrix[i];
    }
    delete[] matrix;
}
