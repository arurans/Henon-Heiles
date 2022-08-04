#include "hamiltonian.h"

/**
 * @brief The hamiltonian for this Hénon Heiles system
 * 
 * @param Y The computed matrix
 * @return vec The hamiltonian as a column vector
 */
Array<double, Dynamic, 1> hamiltonian(const Ref<const Matrix<double, 4, Dynamic>> Y)
{
    return 0.5 * (Y.row(0).array().square() + Y.row(1).array().square())
        +  0.5 * (Y.row(2).array().square() + Y.row(3).array().square())
        +  Y.row(3).array() * Y.row(2).array().square() - 1.0/3.0 * Y.row(3).array().cube();
}