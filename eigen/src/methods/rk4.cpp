#include "rk4.h"

/**
 * @brief 
 * The Hénon Heiles system for the fourth order runge kutta method
 * 
 * @param y The values for computing the next time step
 * @param Y_vec A vector containig the given step of a fourth order computation
 */
void henon_heiles_rk(const Ref<const Array<double, 4, 1>> y, Ref<Array<double, 4, 1>> Y_vec)
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
 * @param y_curr the current values of the system
 * @param Y_vec values to compute the four stages of our method 
 * @param h timestep length
 */
void kutta_iteration(Ref<Array<double, 4, 1>> y_curr, Ref<Matrix<double, 4, 4>> Y_vec, const double& h)
{
    // Convert matrix types to array to enable use of coefficient-wise addition
    henon_heiles_rk(y_curr, Y_vec.col(0).array());
    henon_heiles_rk(y_curr + 0.5*h*Y_vec.col(0).array(), Y_vec.col(1));
    henon_heiles_rk(y_curr + 0.5*h*Y_vec.col(1).array(), Y_vec.col(2));
    henon_heiles_rk(y_curr + h*Y_vec.col(2).array(), Y_vec.col(3));

    y_curr = y_curr + h/6 * (Y_vec.col(0) + 2*Y_vec.col(1) + 2*Y_vec.col(2) + Y_vec.col(3)).array();

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
Matrix<double, 4, Dynamic> kuttas_method(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h)
{
    std::tuple<int, int, int, double> vals = create_H(t_0, t_end, h);
    int n = std::get<0>(vals);
    int m = std::get<1>(vals);
    int skip_storage = std::get<2>(vals);
    double last_step = std::get<3>(vals);

    //Init a matrix to be of the same dimension as the init-cond
    Matrix<double, 4, Dynamic> Y = Matrix<double, 4, Dynamic>::Zero(4, m);
    Y.col(0) = y0;

    //Our method requires four stored values for each step
    Matrix<double, 4, 4> Y_vec = Matrix<double, 4, 4>::Zero(4, 4);

    //Since we're not necessarily storing every iteration in the matrix,
    //we need an array to store the current iteration in time
    Array<double, 4, 1> y_curr = y0;

    //Index to keep count of where to store in matrix
    int storage_index = 1;

    //Compute the system forward in time
    for (int i = 1; i < n - 1; i++)
    {
        kutta_iteration(y_curr, Y_vec, h);
        if (!(i % skip_storage))
        {
            Y.col(storage_index) = y_curr;
            storage_index++;
        }
    }

    //Use last_step as step size to compute the last step
    kutta_iteration(y_curr, Y_vec, last_step);
    Y.col(m-1) = y_curr;

    return Y;
}