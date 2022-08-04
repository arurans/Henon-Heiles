#pragma once

//Störmer-Verlet method of order 2

#include "../utils.h"

void henon_heiles_sv(Ref<Array<double, 4, 1>> y_curr, const double& h, Ref<Array<double, 2, 1>> q_next);
Matrix<double, 4, Dynamic> stormer_verlet(const double& t_0, const double& t_end, const Ref<const Array<double, 4, 1>> y0, const double& h);