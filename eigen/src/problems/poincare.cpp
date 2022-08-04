#include "poincare.h"

/**
 * @brief 
 * Compute the Poincaré map of a given matrix Y
 * 
 * @param Y matrix as a result of an implemented method
 * @return mat the created Poincaré map
 */
Matrix<double, 2, Dynamic> poincare(const Ref<const Matrix<double, 4, Dynamic>> Y)
{
    // First count how many times this occurs, instead of using the time-consuming resize function
    // and to not use unessecary memory by creating an array with "maximum theoretical length"
    int n = 0;
    for (int i = 1; i < Y.cols(); i++)
    {
        if (Y.col(i)[0] > 0 && Y.col(i)[2] * Y.col(i-1)[2] < 0)
            n++;
    }


    // Then perform the actual interpolation, to compute the desired points
    Matrix<double, 2, Dynamic> p_mat = Matrix<double, 2, Dynamic>::Zero(2, n);
    n = 0;
    double lam = 0;
    for (int i = 1; i < Y.cols(); i++)
    {
        if (Y.col(i)[0] > 0 && Y.col(i)[2] * Y.col(i-1)[2] < 0)
        {
            lam = Y.col(i-1)[2]/(Y.col(i-1)[2] - Y.col(i)[2]);
            p_mat.col(n)[0] = lam * Y.col(i)[3] + (1 - lam) * Y.col(i-1)[3];
            p_mat.col(n)[1] = lam * Y.col(i)[1] + (1 - lam) * Y.col(i-1)[1];
            n++;
        }
    }

    return p_mat;
}