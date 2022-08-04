#pragma once

//Shampine-Bogacki method of order 3

#include "../utils.h"

void henon_heiles_sb(const Ref<const Array<double, 4, 1>> y, Ref<Array<double, 4, 1>> Y_vec);
void sb_iteration(Ref<Array<double, 4, 1>> y_curr, Ref<Matrix<double, 4, 3>> Y_vec, const double& h);
Matrix<double, 4, Dynamic> shampine_bogacki(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h);