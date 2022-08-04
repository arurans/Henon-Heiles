#pragma once

//Shampine-Bogacki method of order 3

#include "../utils.h"

void henon_heiles_sb(const vec& y, vec Y_vec);
void sb_iteration(vec Y, const vec& y, mat& Y_vec, const double& h);
mat shampine_bogacki(const double& t_0, const double& t_end, const vec& y0, const double& h);