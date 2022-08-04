#include "sb.h"


/**
 * @brief 
 * The Hénon Heiles system for the Shampine-Bogacki method
 * 
 * @param y The values for computing the next time step
 * @param Y_vec A vector containig the given step of a third order computation
 */
void henon_heiles_sb(const vec& y, vec Y_vec)
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
 * A function to perform the third order stages
 * for each iteration forward in time
 * 
 * @param Y Vector to store the final result in
 * @param y the current values of the system
 * @param Y_vec values to compute the three stages of our method 
 * @param h timestep length
 */
void sb_iteration(vec Y, const vec& y, mat& Y_vec, const double& h)
{
    henon_heiles_sb(y, Y_vec.unsafe_col(0));
    henon_heiles_sb(y + 0.5*h*Y_vec.unsafe_col(0), Y_vec.unsafe_col(1));
    henon_heiles_sb(y + 0.75*h*Y_vec.unsafe_col(1), Y_vec.unsafe_col(2));

    Y = y + h/9 * (2*Y_vec.col(0) + 3*Y_vec.col(1) + 4*Y_vec.col(2));

    return;
}

/**
 * @brief 
 * A third order Runge Kutta method (Shampine-Bogacki)
 * implemented for the Hénon Heiles system
 * 
 * @param t_0 Start time
 * @param t_end End time
 * @param y0 Initial conditions of the system
 * @param h Length of timestep between iterations
 * @return Y matrix
 */
mat shampine_bogacki(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    std::pair<uword, double> vals = create_H(t_0, t_end, h);
    uword n = vals.first;
    double last_step = vals.second;


    //Init a matrix to be of the same dimension as the init-cond
    mat Y = zeros(y0.n_elem, n);
    Y.unsafe_col(0) = y0;

    //Our method requires three stored values for each step
    mat Y_vec = zeros(y0.n_elem, 3);

    //Compute the system forward in time
    for (uword i = 0; i < n - 2; i++)
    {
        sb_iteration(Y.unsafe_col(i+1), Y.col(i), Y_vec, h);
    }

    //Use last_step as step size to compute the last step
    sb_iteration(Y.unsafe_col(n-1), Y.col(n-2), Y_vec, last_step);

    return Y;
}