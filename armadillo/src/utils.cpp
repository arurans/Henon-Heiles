#include "utils.h"


/**
 * @brief 
 * Calculate the number of iterations to compute, with a given h
 * if remainder is not zero the last step is the remainder
 * otherwise the last step is h
 * 
 * @param t_0 Start time for the system
 * @param t_end End time for the system
 * @param h Timestep length for iterations
 * @return std::pair<int, double> return the number of steps and the length of the last time step
 */
std::pair<uword,double> create_H(const double& t_0, const double& t_end, const double& h)
{
    double ratio = (t_end - t_0)/h;
    double remainder = std::remainder(t_end - t_0, h);
    uword n = int(ratio) + 1; //Add one to include the initial condition

    //Check if the remainder is non-zero because of floating point error
    if (remainder >= 1e-10)
        return std::pair<uword, double>(n + 1, remainder); //Add 1 to n to include the last step

    return std::pair<uword, double>(n, h);
}


/**
 * @brief Create initial condition for the system
 * 
 * @param H_0 Initial energy in the system
 * @return the initial condition as a vector
 */
vec create_init_cond(const double& H_0)
{
    double p2 = 0;
    double q1 = 0;
    double q2 = 0.45;
    double p1 = std::sqrt(2.0/3.0*pow(q2, 3) + 2*H_0 - pow(p2, 2) - pow(q1, 2) - pow(q2, 2) - 2*pow(q1,2)*q2);

    return vec({p1, p2, q1, q2});
}


/**
 * @brief 
 * Create a vector containing the time steps
 * 
 * @param t_0 start time
 * @param t_end end time
 * @param h length of each time step
 * @return T Time vector
 */
vec create_T(const double& t_0, const double& t_end, const double& h)
{
    std::pair<uword, double> vals = create_H(t_0, t_end, h);
    uword n = vals.first;
    double last_step = vals.second;

    vec T = zeros(n);
    T[0] = t_0;
    T[n-1] = t_end;
    for (uword i = 0; i < n-2; i++)
    {
        T[i+1] = T[i] + h;
    }

    return T;
}

/**
 * @brief 
 * Purely to have an easier time naming files
 * based on step size for time
 * 
 * @param value length of time step size
 * @return std::string 
 */
std::string decimal_to_string(double h)
{
    std::stringstream ss;
    ss << "_" << h << ".csv";
    std::string decimal = ss.str();

    return decimal;
}