#include <iostream>
#include <iomanip>
#include "calerf.h"
#include <cmath>
#include <fstream>

/**
 *  Application of differential methods to an approximate solution of the diffusion equation with initial condition and boundary conditions.
 *
 *  Algorithms used:
 *  - FTCS Explicit Method
 *  - Laasonen method (Backward Time, Centered Space - BTCS Scheme) using LU decomposition
 *  - Laasonen method using Thomas algorithm (TDMA)"
 *
 *  Initial condition: U(x,0) = 1
 *  First boundary condition: U(rad,t) = 0
 *  Second boundary condition: U(rad+a, t) = 1.0 - (rad/(rad+a))*erfc(a/(2.0*sqrt(D*t)))
 * */

double const D = 1.0; // diffusion coefficient
double const T_0 = 0.0; // starting moment
double const T_MAX = 2.0; // max time
double const rad = 1.0; // radius of the sphere
double const a = 10.0; // domain length
double const X_MIN = rad; // x - spatial coordinate
double const X_MAX = rad + a;
double const lambda_dm = 0.4; // lambda value used for direct methods, NOTE: beware of numerical stability limitations - method stable only when lambda < 0.5
double const lambda_im = 1.0; // lambda value used for for indirect methods, NOTE: lambda = D*(dt/h*h)

void ftcs_method(); // Forward Time Centered Space Scheme
void lm_ta(); // Laasonen method using Thomas algorithm
double get_initial_condition();
double get_first_boundary_condition();
double get_second_boundary_condition(double tk, double dt);
double find_max_error(double tk, double h, int nodes_x, const double *vector);
double analytical_solution(double tk, double xi);
double get_alpha();
double get_beta();
double get_gamma();
double get_phi();
double get_psi();
double get_theta(double tk, double dt);
void determine_eta(double *eta, const double *u, const double *d, const double *l, int nodes_x);
void Thomas_algorithm(double *y, const double *eta, const double *u, const double *l, const double *b, int nodes_x);
double **allocate_matrix(int rows, int cols);
void free_matrix(double *const *matrix, int rows);
void print_array(double *array, int size, int print_width=10);
void fill_array(double *array, int size, double value);
void fill_matrix(double **matrix, int rows, int cols, double value);

int main()
{
    ftcs_method();
//    lm_ta();

    return 0;
}

void lm_ta()
{
    int nodes_x, nodes_t;
    double h, dt;
    double tk, xi;
    double alpha, beta, gamma;
    double phi, psi, theta;
    int i, k;

    std::cout.setf(std::ios::fixed);
    std::ofstream LM_TA_results("LM_TA_results.txt");
    std::ofstream LM_TA_h_errors("LM_TA_h-dependent_errors.txt");
    std::ofstream LM_TA_t_errors("LM_TA_t-dependent_errors.txt");
    if( !LM_TA_results || !LM_TA_h_errors || !LM_TA_t_errors)
    {
        std::cerr << "Error: file could not be opened" << std::endl;
        exit(1);
    }

    LM_TA_results << "h approx_res analytical_res\n";
    LM_TA_h_errors << "h max_err\n";
    LM_TA_t_errors << "t max_err\n";

//    for(nodes_x = 26; nodes_x < 600; nodes_x += 10)
    for(nodes_x = 1001; nodes_x < 1002; nodes_x += 10)
    {
        auto *u = new double[nodes_x - 1]; // upper diagonal of the tridiagonal  matrix
        auto *d = new double[nodes_x]; // principal (main) diagonal of the tridiagonal  matrix
        auto *l = new double[nodes_x - 1]; // lower diagonal of the tridiagonal  matrix
        auto *b = new double[nodes_x];
        auto *eta = new double[nodes_x];
        auto *y = new double[nodes_x]; // solution vector

        h = a / (nodes_x - 1);
        nodes_t = static_cast<int>((T_MAX - T_0) / ((lambda_im * h * h) / D)) + 1;
        dt = T_MAX / (nodes_t - 1);
        tk = T_0;
        fill_array(y, nodes_x, get_initial_condition());

        alpha = get_alpha();
        beta = get_beta();
        phi = get_phi();
        psi = get_psi();
        gamma = get_gamma();

        for (k = 0; k < nodes_t - 1; k++)
        {
            xi = X_MIN;
            theta = get_theta(tk+dt, dt);

            d[0] = -alpha/h + beta;
            u[0] = alpha/h;
            b[0] = -gamma;

            for (i = 1; i < nodes_x-1; i++)
            {
                xi += h;
                d[i] = -(1.0 + 2.0 * lambda_im);
                u[i] = lambda_im * (1.0 + (h/xi));
                l[i-1] = lambda_im * (1.0 - (h/xi));
                b[i] = -y[i];
            }

            d[nodes_x-1] = phi/h + psi;
            l[nodes_x-2] = -phi/h;
            b[nodes_x-1] = theta;

            determine_eta(eta, u, d, l, nodes_x); // determining the vector eta

            Thomas_algorithm(y, eta, u, l, b, nodes_x);

            /**
             * data to plot numerical and analytical solutions for a few selected values of time t from the whole interval t - plotted when nodes_x = 1001
             * */
             xi = X_MIN;
//            if (k == ((nodes_t - 1)/8)) // results for T=0.25
//            if (k == ((nodes_t - 1)/4)) // results for T=0.5
//            if (k == ((nodes_t - 1)/2)) // results in the middle of the time interval
            if (k == (((nodes_t - 1)/4)*3)) // //results for T=1.5
            {
                for (i=1; i<nodes_x-2; i++)
                {
                    xi += h;
                    LM_TA_results << xi << " " << y[i] << " " << analytical_solution(tk+dt, xi) << "\n";
                }
            }

            tk += dt;

            /**
             * results saved to file for plotting dependence of the maximum absolute value of the error observed for optimal h as a function of the time
             */
//            LM_TA_t_errors << tk << " " << find_max_error(tk, h, nodes_x, y) << "\n";
        }

        /**
        * results saved to file for plotting dependence of the maximum absolute value of the error observed for T_MAX as a function of the spatial step h
        **/
        LM_TA_h_errors << log10(h) << " " << log10(find_max_error(T_MAX, h, nodes_x, y)) << "\n";

        delete[] u; delete[] d; delete[] l; delete[] b; delete[] eta; delete[] y;
    }

    LM_TA_results.close();
    LM_TA_h_errors.close();
    LM_TA_t_errors.close();
}

void determine_eta(double *eta, const double *u, const double *d, const double *l, int nodes_x)
{
    eta[0] = d[0];
    for (int i = 1; i < nodes_x; i++)
    {
        eta[i] = d[i] - l[i-1] * (1.0 / eta[i-1]) * u[i-1];
    }
}

void Thomas_algorithm(double *y, const double *eta, const double *u, const double *l, const double *b, int nodes_x)
{
    int i;
    auto *r = new double[nodes_x];

    /**
     * determining the vector r
     * */
    r[0] = b[0];
    for(i=1; i<nodes_x; i++)
    {
        r[i] = b[i] - l[i - 1] * (1.0 / eta[i - 1]) * r[i - 1];
    }

    /**
     * determining the vector y
     * */
    y[nodes_x-1] = (1.0/eta[nodes_x-1]) * r[nodes_x-1];
    for(i=nodes_x-2; i>=0; i--)
    {
        y[i] = (1.0/eta[i]) * (r[i] - (u[i]*y[i+1]));
    }

    delete[] r;
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

    FTCS_results << "h approx_res analytical_res\n";
    FTCS_h_errors << "h max_err\n";
    FTCS_t_errors << "t max_err\n";

    for(nodes_x = 1001; nodes_x<1002; nodes_x+=50) // 1001 nodes -> h=0.01
//    for(nodes_x = 22; nodes_x<601; nodes_x+=10)
    {
        h = a / (nodes_x - 1);
        nodes_t = static_cast<int>((T_MAX-T_0) / ((lambda_dm * h * h) / D)) + 1;
        dt = T_MAX/(nodes_t-1);

        double **u = allocate_matrix(nodes_t, nodes_x); // matrix of approximate values
        fill_array(u[0], nodes_x, get_initial_condition()); // filling values from initial condition

        tk = T_0;
        for(k = 0; k < nodes_t - 1; k++)
        {
            xi = X_MIN + h;
            for(i = 1; i < nodes_x - 1; i++) // first and last i omitted - values already determined from the initial condition
            {
                res_u = u[k + 1][i] = lambda_dm * u[k][i - 1]*(1.0 - h/xi) + (1.0 - 2.0*lambda_dm) * u[k][i] + lambda_dm * u[k][i + 1]*(1.0 + h/xi);
                res_a = analytical_solution(tk + dt, xi);

                /**
                 * data to plot numerical and analytical solutions for a few selected values of time t from the whole interval t - plotted using 1001 x nodes
                 * */
                if (k == ((nodes_t - 1)/8)) // results for T=0.25
//                if (k == ((nodes_t - 1)/4)) // results for T=0.5
//                if (k == (((nodes_t - 1)/4)*3)) // //results for T=1.5
//                if (k == ((nodes_t - 1)/2)) // results in the middle of the time interval
                    {FTCS_results << xi << " " << res_u << " " << res_a << "\n";}

                xi += h;
            }
            u[k + 1][0] = get_first_boundary_condition();
            u[k + 1][nodes_x - 1] = get_second_boundary_condition(tk, dt);
            tk += dt;

            /**
             * results saved to file for plotting dependence of the maximum absolute value of the error observed for optimal h as a function of the time
             */
            FTCS_t_errors << tk << " " << find_max_error(tk, h, nodes_x, u[k+1]) << "\n";
        }

        /**
         * results saved to file for plotting dependence of the maximum absolute value of the error observed for T_MAX as a function of the spatial step h
         **/
        FTCS_h_errors << log10(h) << " " << log10(find_max_error(T_MAX, h, nodes_x, u[k])) << "\n";

        free_matrix(u, nodes_t);
    }

    FTCS_h_errors.close();
    FTCS_t_errors.close();
    FTCS_results.close();
}

double get_alpha() { return 0.0; }

double get_beta() { return 1.0; }

double get_gamma() { return 0.0; }

double get_phi() { return 0.0; }

double get_psi() { return 1.0; }

double get_theta(double tk, double dt) { return 1.0 - (rad / (rad + a)) * static_cast<double>(calerf::ERFCL(a / (2.0 * sqrt(D * tk)))); }

double find_max_error(double tk, double h, int nodes_x, const double *vector)
{
    double max_error = 0.0;
    double error;
    double xi = X_MIN;
    double res_v, res_a;

    for(int i = 1; i < nodes_x-1; i++)
    {
        xi += h;
        res_v = vector[i];
        res_a = analytical_solution(tk, xi);
        error = fabs(res_v - res_a);
        if(error > max_error) { max_error = error; }
    }

    return max_error;
}

double analytical_solution(double tk, double xi) { return 1.0 - (rad / xi) * static_cast<double>(calerf::ERFCL((xi - rad) / (2.0 * sqrt(D * tk)))); }

double get_initial_condition() { return 1.0; }

double get_first_boundary_condition() { return 0.0; }

double get_second_boundary_condition(double tk, double dt) { return 1.0 - (rad / (rad + a)) * static_cast<double>(calerf::ERFCL(a / (2.0 * sqrt(D * (tk + dt))))); }

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
