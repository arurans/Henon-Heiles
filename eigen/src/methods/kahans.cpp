#include "kahans.h"

/**
 * @brief 
 * Given an identity matrix fill in the rest of the values of A
 * 
 * @param y current computed values
 * @param h length of timestep
 * @param A matrix to be filled
 */
void create_A(const Ref<const Array<double, 4, 1>> y, const double& h, Ref<Matrix<double, 4, 4>> A)
{
    A.col(2)[0] = h*(y[3] + 0.5);
    A.col(3)[0] = h*y[2];
    A.col(2)[1] = h*y[2];
    A.col(3)[1] = h*(0.5 - y[3]);
    A.col(0)[2] = -0.5*h;
    A.col(1)[3] = -0.5*h;
}

/**
 * @brief
 * Given an array with the same length as y, fill in for b
 * 
 * @param y current computed values
 * @param h length of timestep
 * @param b vector to be filled
 */
void create_b(const Ref<const Array<double, 4, 1>> y, const double& h, Ref<Matrix<double, 4, 1>> b)
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
 * @param y_curr the current values of the system
 * @param h timestep length
 * @param A matrix for solving linear system
 * @param b vector for solving linear system
 */
void kahans_iteration(const Ref<const Array<double, 4, 1>> y_curr, const double& h, Ref<Matrix<double, 4, 4>> A, Ref<Matrix<double, 4, 1>> b)
{
    create_A(y_curr, h, A);
    create_b(y_curr, h, b);
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
Matrix<double, 4, Dynamic> kahans(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h)
{
    std::tuple<int, int, int, double> vals = create_H(t_0, t_end, h);
    int n = std::get<0>(vals);
    int m = std::get<1>(vals);
    int skip_storage = std::get<2>(vals);
    double last_step = std::get<3>(vals);

    //Init a matrix to be of the same dimension as the init-cond
    Matrix<double, 4, Dynamic> Y = Matrix<double, 4, Dynamic>::Zero(4, m);
    Y.col(0) = y0;

    //Matrix and vector used for solving linear system each step
    Matrix<double, 4, 4> A = Matrix<double, 4, 4>::Identity(4, 4);
    Matrix<double, 4, 1> b = Array<double, 4, 1>::Zero(4);

    //Since we're not necessarily storing every iteration in the matrix,
    //we need an array to store the current iteration in time
    Matrix<double, 4, 1> y_curr = y0;

    //Index to keep count of where to store in matrix
    int storage_index = 1;

    //Compute the system forward in time
    for (int i = 1; i < n - 1; i++)
    {
        kahans_iteration(y_curr, h, A, b);
        y_curr = A.partialPivLu().solve(b);
        if (!(i % skip_storage))
        {
            Y.col(storage_index) = y_curr;
            storage_index++;
        }
    }

    //Use last_step as step size to compute the last step
    kahans_iteration(y_curr, last_step, A, b);
    Y.col(m-1) = A.partialPivLu().solve(b);

    return Y;
}