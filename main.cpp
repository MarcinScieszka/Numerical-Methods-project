#include <iostream>
#include <iomanip>
#include "calerf.h"
#include <cmath>
#include <fstream>

/**
 *  Application of differential methods to an approximate solution of the diffusion equation with initial condition and boundary conditions.
 *
 *  Algorithms used:
 *  - FTCS Method
 *  - Laasonen method using LU decomposition
 *  - Laasonen method using Thomas algorithm
 *
 *  Initial condition: U(x,0) = 1
 *  First boundary condition: U(r,t) = 0
 *  Second boundary condition: U(r+a, t) = 1.0 - (r/(r+a))*erfc(a/(2.0*sqrt(D*t)))
 * */

double const D = 1.0; // diffusion coefficient
double const T_0 = 0.0; // starting moment
double const T_MAX = 2.0; // max time
double const r = 1.0; // radius of the sphere
double const a = 10.0; // domain length
double const X_MIN = r; // x - spatial coordinate
double const X_MAX = r + a;
double const lambda_dm = 0.4; // lambda value used for direct methods, NOTE: beware of numerical stability limitations - method stable only when lambda < 0.5
double const lambda_im = 1.0; // lambda value used for for indirect methods, NOTE: lambda = D*(dt/h*h)

double **allocate_matrix(int rows, int cols);
void print_matrix(double *const *matrix, int rows, int cols);
void free_matrix(double *const *matrix, int rows);
void print_array(double *array, int size);
void fill_array(double *array, int size, double value);
void fill_matrix(double **matrix, int rows, int cols, double value);
double analytical_solution(double tk, double xi);
double find_max_error(double tk, double h, int nodes_x, int k, double *const *u);
void ftcs_method();

int main()
{
    ftcs_method();

    return 0;
}

void ftcs_method()
{
    double dt; // delta t
    double h; // step h on the spatial grid, delta x
    double tk; // current time on the time grid, NOTE: dt = T_MAX/(nodes_t-1)
    double xi; // current point on the spatial grid
    int i=0, k=0;

    std::ofstream result_file("results_FTCS.txt");
    if(!result_file)
    {
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    int nodes_x; // nr of nodes on the time grid, NOTE: for 101 nodes h=0.1
    int nodes_t; // nr of nodes on the spacial grid, NOTE: for 21 nodes dt=0.1

    for(nodes_x = 21; nodes_x<101; nodes_x+=5)
    {
        auto *ic_array = new double[nodes_x];  // initial condition array
        fill_array(ic_array, nodes_x, 1.0); // filling values from initial condition

        h = a / (nodes_x - 1);
        dt = (lambda_dm * h * h) / D; //NOTE: to calculate dt based on current spacial step h we use formula: dt=(lambda*h*h)/D
        nodes_t = int(T_MAX / dt);
        std::cout << std::endl << "dt: " << dt << ", " << "nr of t nodes: " << nodes_t << std::endl;
        std::cout << "h: " << h << ", " << "nr of x nodes: " << nodes_x << std::endl;

        double **u = allocate_matrix(nodes_t, nodes_x); // matrix of approximate values
        fill_matrix(u, nodes_t, nodes_x, 0.0);

        for (i = 0; i < nodes_x; i++) u[0][i] = ic_array[i];

        for (k = 0; k < nodes_t - 1; k++) {
            for (i = 1; i < nodes_x -
                            2; i++) // we omit first and last i - values there are already determined from the initial condition
            {
                u[k + 1][i] = lambda_dm * u[k][i - 1] + (1 - lambda_dm) * u[k][i] +
                              lambda_dm * u[k][i + 1]; // approximate solution at node x_i at time level t_k+1
            }
            u[k + 1][0] = 0.0; // first boundary condition
            u[k + 1][nodes_x - 1] =
                    1.0 - (r / (r + a)) * calerf::ERFCL(a / (2.0 * sqrt(D * (k + 1)))); // second boundary condition
        }
        /**
         * dependence of the maximum absolute value of the error observed for T_MAX as a function of the spatial step h
         * TODO: save to file results for T_MAX and range of h's
         **/
        tk = T_MAX;
        k = nodes_t - 1;
        //(double tk, double h, int nodes_x, int k, double *const *u)
        double max_error = find_max_error(tk, h, nodes_x, k, u);
        std::cout << "max error: " << max_error << std::endl;

        delete[] ic_array;
        free_matrix(u, nodes_t);
    }

    result_file.close();
}

double find_max_error(double tk, double h, int nodes_x, int k, double *const *u)
{
    double max_error = 0.0;
    double error;
    double xi = X_MIN + h;

    for (int i=1; i<nodes_x-2; i++)
    {
        error = fabs(u[k][i]- analytical_solution(tk, xi));
//        std::cout << "err: i=" << i << " : "<< error << " \n";
        if (error > max_error)
        {
            max_error = error;
//            std::cout << "max:" << max_error << " \n";
        }
        xi += h;
    }

    return max_error;
}

double analytical_solution(double tk, double xi)
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
            std::cout << matrix[i][j] << std::setw(8);
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
