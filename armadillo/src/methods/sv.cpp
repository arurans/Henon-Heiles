#include "sv.h"

/**
 * @brief 
 * Hénon Heiles system for the Störmer-Verlet method (order 2)
 * 
 * @param y current values
 * @param h length of timestep
 * @param Y_vec new values
 * @param q_next computing helpers
 */
void henon_heiles_sv(const vec& y, const double& h, vec Y_vec, vec& q_next)
{
    double p1_half = y[0] + q_next[0];
    double q1_next = y[2] + h * p1_half;
    double p2_half = y[1] + q_next[1];
    double q2_next = y[3] + h * p2_half;

    q_next[0] = 0.5 * h * (-q1_next*(1 + 2*q2_next));
    q_next[1] = 0.5 * h * (-q2_next - pow(q1_next, 2) + pow(q2_next, 2));

    Y_vec[0] = p1_half + q_next[0];
    Y_vec[1] = p2_half + q_next[1];
    Y_vec[2] = q1_next;
    Y_vec[3] = q2_next;

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
mat stormer_verlet(const double& t_0, const double& t_end, const vec& y0, const double& h)
{
    std::pair<uword, double> vals = create_H(t_0, t_end, h);
    uword n = vals.first;
    double last_step = vals.second;

    //Init a matrix to be of the same dimension as the init-cond
    mat Y = zeros(y0.n_elem, n);
    Y.unsafe_col(0) = y0;

    vec q_next({0.5 * h * (-y0[2]*(1 + 2*y0[3])), 0.5 * h * (-y0[3] - pow(y0[2], 2) + pow(y0[3], 2))});

    for (uword i = 0; i < n - 2; i++)
        henon_heiles_sv(Y.unsafe_col(i), h, Y.unsafe_col(i+1), q_next);

    vec y = Y.col(n-2);
    q_next = {0.5 * last_step * (-y[2]*(1 + 2*y[3])), 0.5 * last_step * (-y[3] - pow(y[2], 2) + pow(y[3], 2))};
    henon_heiles_sv(Y.unsafe_col(n-2), last_step, Y.unsafe_col(n-1), q_next);

    return Y;
}