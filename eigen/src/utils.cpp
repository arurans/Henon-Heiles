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
 * @return std::tuple<int, int, int, double> number of total iterations, number of elements to store in matrix, SKIP_STORAGE, and last step size
 */
std::tuple<int, int, int, double> create_H(const double& t_0, const double& t_end, const double& h)
{
    double remainder = std::remainder(t_end - t_0, h);

    //If there is a remainder, the last step is equal to it
    double last_step = (remainder >= 1e-10) ? remainder : h;

    //The total number of iterations
    //Ceil to include the initial condition and +1 to include last step
    double ratio = (t_end - t_0)/h;
    int n = int(std::ceil(ratio)) + 1;

    //If all values are to be stored, no more calculations are necessary
    if (SKIP_STORAGE == 1)
        return std::tuple<int, int, int, double>(n, n, SKIP_STORAGE, last_step);

    //Compute the number of iterations to store, to allocate array memory
    int m = int(std::ceil(n/SKIP_STORAGE));

    //If the time interval is not divisible by time step size, increase by 1 to include last step
    if (remainder >= 1e-10)
    {
        m++;
    } else 
    //If time interval is divisible by time step size
    {
        //If the number of iterations (excluding last one) is divisible by skip_storage
        //Increase by one for last step. If it is not divisible, then add 2 to include 
        //step before last
        double storage_remainder = std::remainder(n-1, SKIP_STORAGE);
        m = (storage_remainder >= 1e-10) ? m + 2 : m + 1;
    }

    return std::tuple<int, int, int, double>(n, m, SKIP_STORAGE, last_step);
}

/**
 * @brief Create initial condition for the system
 * 
 * @param H_0 Initial energy in the system
 * @return the initial condition as a vector
 */
Array<double, 4, 1> create_init_cond(const double& H_0)
{
    double p2 = 0;
    double q1 = 0;
    double q2 = 0.45;
    double p1 = std::sqrt(2.0/3.0*pow(q2, 3) + 2*H_0 - pow(p2, 2) - pow(q1, 2) - pow(q2, 2) - 2*pow(q1,2)*q2);

    return Array<double, 4, 1>(p1, p2, q1, q2);
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
Array<double, Dynamic, 1> create_T(const double& t_0, const double& t_end, const double& h)
{
    std::tuple<int, int, int, double> vals = create_H(t_0, t_end, h);
    int n = std::get<0>(vals);
    int m = std::get<1>(vals);
    int skip_storage = std::get<2>(vals);

    Array<double, Dynamic, 1> T = Array<double, Dynamic, 1>::Zero(m);
    T[0] = t_0;

    //Current time step
    double curr_time = 0;
    int storage_index = 1;
    for (int i = 1; i < n-1; i++)
    {
        curr_time += h;
        if (!(i % skip_storage))
        {
            T[storage_index] = curr_time;
            storage_index++;
        }
    }
    T[m-1] = t_end;

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

/**
 * @brief 
 * Function to save Eigen matrices and arrays to CSV
 * 
 * @param filename file name
 * @param M matrix/array to be saved to CSV
 */
void matrix_to_CSV(std::string filename, const Ref<const Matrix<double, Dynamic, Dynamic>> Y)
{
    // Format the way matrix is saved
    Eigen::IOFormat CleanFmt(10, 0, ", ", "\n");

    // Save matrix to file
    std::ofstream file(filename);
    file << Y.format(CleanFmt);
    file.close();

    return;
}