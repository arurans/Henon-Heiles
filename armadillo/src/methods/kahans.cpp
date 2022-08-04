#include "kahans.h"

/**
 * @brief 
 * Given an identity matrix fill in the rest of the values of A
 * 
 * @param y current computed values
 * @param h length of timestep
 * @param A matrix to be filled
 */
void create_A(const vec& y, const double& h, mat& A)
{
    A.unsafe_col(2)[0] = h*(y[3] + 0.5);
    A.unsafe_col(3)[0] = h*y[2];
    A.unsafe_col(2)[1] = h*y[2];
    A.unsafe_col(3)[1] = h*(0.5 - y[3]);
    A.unsafe_col(0)[2] = -0.5*h;
    A.unsafe_col(1)[3] = -0.5*h;
}

/**
 * @brief
 * Given an array with the same length as y, fill in for b
 * 
 * @param y current computed values
 * @param h length of timestep
 * @param b vector to be filled
 */
void create_b(const vec& y, const double& h, vec& b)
{
    b[0] = y[0] - 0.5*h*y[2];
    b[1] = y[1] - 0.5*h*y[3];
    b[2] = y[2] + 0.5*h*y[0];
    b[3] = y[3] + 0.5*h*y[1];
}

/**
 * @brief 
 * Perform iteration of Kahan's method
 * 
 * @param y current values
 * @param h timestep length
 * @param A matrix for solving linear system
 * @param b vector for solving linear system
 */
void kahans_iteration(const vec& y, const double& h, mat& A, vec& b)
{
    create_A(y, h, A);
    create_b(y, h, b);
}

/**
 * @brief 
 * Kahan's method (implicit method of order 2)
 * 
 * @param t_0 start time
 * @param t_end end time
 * @param y0 initial condition
 * @param h timestep length
 * @return mat solution for all time
 */
mat kahans(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    std::pair<uword, double> vals = create_H(t_0, t_end, h);
    uword n = vals.first;
    double last_step = vals.second;

    uword m = y0.n_elem; 

    //Init a matrix to be of the same dimension as the init-cond
    mat Y = zeros(m, n);
    Y.unsafe_col(0) = y0;

    //Matrix and vector used for solving linear system each step
    mat A(m, m, arma::fill::eye);
    vec b(m, arma::fill::zeros);

    //Our method requires four stored values for each step
    mat Y_vec = zeros(m, 4);

    //Compute the system forward in time
    for (uword i = 0; i < n - 2; i++)
    {
        kahans_iteration(Y.unsafe_col(i), h, A, b);
        Y.unsafe_col(i+1) = solve(A, b, arma::solve_opts::fast);
    }

    //Use last_step as step size to compute the last step
    kahans_iteration(Y.unsafe_col(n-2), last_step, A, b);
    Y.unsafe_col(n-1) = arma::solve(A, b, arma::solve_opts::fast);

    return Y;
}