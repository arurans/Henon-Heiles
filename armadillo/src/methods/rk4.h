#pragma once

//Kutta's method (fourth order Runge Kutta method)

#include "../utils.h"


void henon_heiles_rk(const vec& y, vec Y_vec);
void kutta_iteration(vec Y, const vec& y, mat& Y_vec, const double& h);
mat kuttas_method(const double& t_0, const double& t_end, const vec& y0, const double& h);