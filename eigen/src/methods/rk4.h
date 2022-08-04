#pragma once

//Kutta's method (fourth order Runge Kutta method)

#include "../utils.h"


// We know that the dimension of our problem is 4, and Eigen is much quicker when smaller matrices are
// defined with dimension, as it will create a normal C-array, as opposed to dynamically allocating memory
void henon_heiles_rk(const Ref<const Array<double, 4, 1>> y, Ref<Array<double, 4, 1>> Y_vec);
void kutta_iteration(Ref<Array<double, 4, 1>> y_curr, Ref<Matrix<double, 4, 4>> Y_vec, const double& h);
Matrix<double, 4, Dynamic> kuttas_method(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h);