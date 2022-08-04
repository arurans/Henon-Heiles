#include "sv.h"

/**
 * @brief 
 * Hénon Heiles system for the Störmer-Verlet method (order 2)
 * 
 * @param y_curr current values of the system
 * @param h length of timestep
 * @param q_next computing helpers
 */
void henon_heiles_sv(Ref<Array<double, 4, 1>> y_curr, const double& h, Ref<Array<double, 2, 1>> q_next)
{
    double p1_half = y_curr[0] + q_next[0];
    double q1_next = y_curr[2] + h * p1_half;
    double p2_half = y_curr[1] + q_next[1];
    double q2_next = y_curr[3] + h * p2_half;

    q_next[0] = 0.5 * h * (-q1_next*(1 + 2*q2_next));
    q_next[1] = 0.5 * h * (-q2_next - pow(q1_next, 2) + pow(q2_next, 2));

    y_curr[0] = p1_half + q_next[0];
    y_curr[1] = p2_half + q_next[1];
    y_curr[2] = q1_next;
    y_curr[3] = q2_next;

    return;
}

/**
 * @brief 
 * The Störmer-Verlet method
 * 
 * @param t_0 start time
 * @param t_end end time
 * @param y0 initial condition
 * @param h length of timestep
 * @return mat 
 */
Matrix<double, 4, Dynamic> stormer_verlet(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h)
{
    std::tuple<int, int, int, double> vals = create_H(t_0, t_end, h);
    int n = std::get<0>(vals);
    int m = std::get<1>(vals);
    int skip_storage = std::get<2>(vals);
    double last_step = std::get<3>(vals);

    //Init a matrix to be of the same dimension as the init-cond
    Matrix<double, 4, Dynamic> Y = Matrix<double, 4, Dynamic>::Zero(4, m);
    Y.col(0) = y0;

    Array<double, 2, 1> q_next(
        0.5 * h * (-y0[2]*(1 + 2*y0[3])), 
        0.5 * h * (-y0[3] - pow(y0[2], 2) + pow(y0[3], 2))
    );

    //Since we're not necessarily storing every iteration in the matrix,
    //we need an array to store the current iteration in time
    Array<double, 4, 1> y_curr = y0;

    //Index to keep count of where to store in matrix
    int storage_index = 1;

    //Compute the system forward in time
    for (int i = 1; i < n - 1; i++)
    {
        henon_heiles_sv(y_curr, h, q_next);
        if (!(i % skip_storage))
        {
            Y.col(storage_index) = y_curr;
            storage_index++;
        }
    }

    //Use last_step as step size to compute the last step
    q_next << 0.5 * last_step * (-y_curr[2]*(1 + 2*y_curr[3])), 
              0.5 * last_step * (-y_curr[3] - pow(y_curr[2], 2) + pow(y_curr[3], 2));

    henon_heiles_sv(y_curr, last_step, q_next);
    Y.col(m-1) = y_curr;

    return Y;
}