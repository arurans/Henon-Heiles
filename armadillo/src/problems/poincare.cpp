#include "poincare.h"

/**
 * @brief 
 * Compute the Poincaré map of a given matrix Y
 * 
 * @param Y matrix as a result of an implemented method
 * @return mat the created Poincaré map
 */
mat poincare(const mat& Y)
{
    // First count how many times this occurs, instead of using the time-consuming resize function
    // and to not use unessecary memory by creating an array with "maximum theoretical length"
    uword n = 0;
    for (uword i = 1; i < Y.n_cols; i++)
    {
        if (Y.unsafe_col(i)[0] > 0 && Y.unsafe_col(i)[2] * Y.unsafe_col(i-1)[2] < 0)
            n++;
    }


    // Then perform the actual interpolation, to compute the desired points
    mat p_mat = zeros(2, n);
    n = 0;
    double lam = 0;
    for (uword i = 1; i < Y.n_cols; i++)
    {
        if (Y.unsafe_col(i)[0] > 0 && Y.unsafe_col(i)[2] * Y.unsafe_col(i-1)[2] < 0)
        {
            lam = Y.unsafe_col(i-1)[2]/(Y.unsafe_col(i-1)[2] - Y.unsafe_col(i)[2]);
            p_mat.unsafe_col(n) = {
                lam * Y.unsafe_col(i)[3] + (1 - lam) * Y.unsafe_col(i-1)[3],
                lam * Y.unsafe_col(i)[1] + (1 - lam) * Y.unsafe_col(i-1)[1]
            };
            n++;
        }
    }

    return p_mat;
}