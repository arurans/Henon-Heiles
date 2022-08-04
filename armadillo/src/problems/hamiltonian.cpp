#include "hamiltonian.h"

/**
 * @brief The hamiltonian for this Hénon Heiles system
 * 
 * @param Y The computed matrix
 * @return vec The hamiltonian as a column vector
 */
vec hamiltonian(const mat& Y)
{
    arma::rowvec v = 0.5 * (arma::pow(Y.row(0), 2) + arma::pow(Y.row(1), 2))
        +  0.5 * (arma::pow(Y.row(2), 2) + arma::pow(Y.row(3), 2))
        +  Y.row(3) % arma::pow(Y.row(2), 2) - 1.0/3.0 * arma::pow(Y.row(3), 3);

    return v.t(); //Return as column vector
}