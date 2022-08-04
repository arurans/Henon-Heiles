#pragma once

//Kahans method of order 2

#include "../utils.h"
#include <eigen3/Eigen/LU>

void create_A(const Ref<const Array<double, 4, 1>> y, const double& h, Ref<Matrix<double, 4, 4>> A);
void create_b(const Ref<const Array<double, 4, 1>> y, const double& h, Ref<Matrix<double, 4, 1>> b);
void kahans_iteration(const Ref<const Array<double, 4, 1>> y_curr, const double& h, Ref<Matrix<double, 4, 4>> A, Ref<Matrix<double, 4, 1>> b);
Matrix<double, 4, Dynamic> kahans(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h);