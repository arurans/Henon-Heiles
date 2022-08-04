#include "rk4.h"

/**
 * @brief 
 * The Hénon Heiles system for the fourth order runge kutta method
 * 
 * @param y The values for computing the next time step
 * @param Y_vec A vector containig the given step of a fourth order computation
 */
void henon_heiles_rk(const vec& y, vec Y_vec)
{
    double y2 = y[2];
    double y3 = y[3];

    Y_vec[0] = -y2*(1 + 2*y3);
    Y_vec[1] = -(y3 + pow(y2, 2) - pow(y3, 2));
    Y_vec[2] = y[0];
    Y_vec[3] = y[1];

    return;
}

/**
 * @brief 
 * A function to perform the fourth order stages
 * for each iteration forward in time
 * 
 * @param Y Vector to store the final result in
 * @param y the current values of the system
 * @param Y_vec values to compute the four stages of our method 
 * @param h timestep length
 */
void kutta_iteration(vec Y, const vec& y, mat& Y_vec, const double& h)
{
    henon_heiles_rk(y, Y_vec.unsafe_col(0));
    henon_heiles_rk(y + 0.5*h*Y_vec.unsafe_col(0), Y_vec.unsafe_col(1));
    henon_heiles_rk(y + 0.5*h*Y_vec.unsafe_col(1), Y_vec.unsafe_col(2));
    henon_heiles_rk(y + h*Y_vec.unsafe_col(2), Y_vec.unsafe_col(3));

    Y = y + h/6 * (Y_vec.unsafe_col(0) + 2*Y_vec.unsafe_col(1) + 2*Y_vec.unsafe_col(2) + Y_vec.unsafe_col(3));

    return;
}

/**
 * @brief 
 * A fourth order Runge Kutta method (Kutta's method of order 4)
 * implemented for the Hénon Heiles system
 * 
 * @param t_0 Start time
 * @param t_end End time
 * @param y0 Initial conditions of the system
 * @param h Length of timestep between iterations
 * @return Y matrix
 */
mat kuttas_method(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    std::pair<uword, double> vals = create_H(t_0, t_end, h);
    uword n = vals.first;
    double last_step = vals.second;


    //Init a matrix to be of the same dimension as the init-cond
    mat Y = zeros(y0.n_elem, n);
    Y.unsafe_col(0) = y0;

    //Our method requires four stored values for each step
    mat Y_vec = zeros(y0.n_elem, 4);

    //Compute the system forward in time
    for (uword i = 0; i < n - 2; i++)
    {
        kutta_iteration(Y.unsafe_col(i+1), Y.unsafe_col(i), Y_vec, h);
    }

    //Use last_step as step size to compute the last step
    kutta_iteration(Y.unsafe_col(n-1), Y.unsafe_col(n-2), Y_vec, last_step);

    return Y;
}