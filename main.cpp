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
void free_matrix(double *const *matrix, int rows);
void print_array(double *array, int size, int print_width=6);
void fill_array(double *array, int size, double value);
void fill_matrix(double **matrix, int rows, int cols, double value);
double analytical_solution(double tk, double xi);
double second_boundary_condition(double tk, double dt);
double first_boundary_condition();
double find_max_error(double tk, double h, int nodes_x, int k, double *const *u);
void ftcs_method();

int main()
{
    ftcs_method();

    return 0;
}

void ftcs_method()
{
    int nodes_x; // nr of nodes on the time grid, NOTE: for 201 nodes h=0.05
    int nodes_t; // nr of nodes on the spacial grid, NOTE: for 201 nodes dt=0.01
    double h; // step h on the spatial grid, delta x
    double dt; // delta t
    double tk; // current time on the time grid
    double xi; // current position on the spatial grid
    double res_u; // approximate result at node x_i at time level t_k+1
    double res_a; // analytical result
    int i, k;

    std::cout.setf(std::ios::fixed);
    std::ofstream FTCS_h_errors("FTCS_h-dependent_errors.txt");
    std::ofstream FTCS_t_errors("FTCS_t-dependent_errors.txt");
    std::ofstream FTCS_results("FTCS_results.txt");
    if(!FTCS_h_errors || !FTCS_results || !FTCS_t_errors)
    {
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    FTCS_h_errors << "h max_err\n";
    FTCS_results << "h approx_res analytical_res\n";

//    for(nodes_x = 31; nodes_x<1501; nodes_x+=10)
    for(nodes_x = 21; nodes_x<22; nodes_x+=50) // 1001 nodes -> h=0.01
    {
        h = a / (nodes_x - 1);
        nodes_t = static_cast<int>((T_MAX-T_0) / ((lambda_dm * h * h) / D)) + 1;
        dt = T_MAX/(nodes_t-1);

        double **u = allocate_matrix(nodes_t, nodes_x); // matrix of approximate values
        fill_array(u[0], nodes_x, 1.0); // filling values from initial condition

        tk = T_0;
        for(k = 0; k < nodes_t - 1; k++)
        {
            xi = X_MIN;
            u[k + 1][0] = first_boundary_condition();
            xi += h;
            for(i = 1; i < nodes_x - 1; i++) // first and last i omitted - values already determined from the initial condition
            {
//                res_u = u[k + 1][i] = lambda_dm * u[k][i - 1] + (1.0 - 2.0*lambda_dm) * u[k][i] + lambda_dm * u[k][i + 1];
                u[k + 1][i] = u[k][i] + lambda_dm * (u[k][i - 1] - 2.0*u[k][i] + u[k][i + 1]);
                res_u = u[k + 1][i];
                res_a = analytical_solution(tk+dt, xi);

                /**
                 * data to plot numerical and analytical solutions for a few selected values of time t from the whole interval t - plotted using 1001 x nodes
                 * */
//                if (k == 0) {FTCS_results << xi << " " << res_u << " " << res_a << "\n";} //results for T_0
                if (k == ((nodes_t - 1)/2)) {FTCS_results << xi << " " << res_u << " " << res_a << "\n";} // results in the middle of the time interval
//                if (k == nodes_t - 2) {FTCS_results << xi << " " << res_u << " " << res_a << "\n";} //results for T_MAX
//                if (k == 12500) {FTCS_results << xi << " " << res_u << " " << res_a << "\n";} //results for T=0.5
//                if (k == 37500) {FTCS_results << xi << " " << res_u << " " << res_a << "\n";} //results for T=1.5

                xi += h;
            }
            u[k + 1][nodes_x - 1] = second_boundary_condition(tk, dt);
            tk += dt;

            /**
             * results saved to file for plotting dependence of the maximum absolute value of the error observed for optimal h as a function of the time
             */
            FTCS_t_errors << tk << " " << find_max_error(tk, h, nodes_x, k+1, u) << "\n";
        }

        /**
         * results saved to file for plotting dependence of the maximum absolute value of the error observed for T_MAX as a function of the spatial step h
         **/
        FTCS_h_errors << log10(h) << " " << log10(find_max_error(tk, h, nodes_x, k, u)) << "\n";

        free_matrix(u, nodes_t);
    }

    FTCS_h_errors.close();
    FTCS_t_errors.close();
    FTCS_results.close();
}

double find_max_error(double tk, double h, int nodes_x, int k, double *const *u)
{
    double max_error = 0.0;
    double error;
    double xi = X_MIN + h;
    double res_u, res_a;

    for(int i=1; i<nodes_x-1; i++)
    {
//        error = fabs(u[k][i]- analytical_solution(tk, xi));
        res_u = u[k][i];
        res_a = analytical_solution(tk, xi);
        error = fabs(res_u - res_a);
//        std::cout << "err: i=" << i << " : "<< error << " \n";
        if(error > max_error)
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
    return 1.0 - (r/xi)*static_cast<double>(calerf::ERFCL((xi-r)/(2.0*sqrt(D*tk))));
}

double first_boundary_condition() { return 0.0; }

double second_boundary_condition(double tk, double dt) { return 1.0 - (r / (r + a)) * static_cast<double>(calerf::ERFCL(a / (2.0 * sqrt(D * (tk + dt))))); }

void fill_array(double *array, int size, double value)
{
    for(int i=0; i<size; i++)
    {
        array[i] = value;
    }
}

void print_array(double *array, int size, int print_width)
{
    std::cout << "[ ";
    for(int i=0; i<size; i++)
    {
        std::cout << array[i]  << std::setw(print_width);
    }
    std::cout << "] \n";
}

double **allocate_matrix(int rows, int cols)
{
    auto** matrix = new double*[rows];
    for(int i=0; i<rows; ++i)
        matrix[i] = new double[cols];
    return matrix;
}

void fill_matrix(double **matrix, int rows, int cols, double value)
{
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            matrix[i][j] = value;
        }
    }
}

void print_matrix(double *const *matrix, const int rows, const int cols, int print_width)
{
    for(int i=rows-1; i>=0; i--)
    {
        for(int j=0; j<cols; j++)
        {
            std::cout << matrix[i][j] << std::setw(print_width);
        }
        std::cout  << "\n";
    }
    std::cout  << "\n";
}

void free_matrix(double *const *matrix, int rows)
{
    for(int i=0; i<rows; ++i)
    {
        delete [] matrix[i];
    }
    delete[] matrix;
}
