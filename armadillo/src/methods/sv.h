#pragma once

//Störmer-Verlet method of order 2

#include "../utils.h"

void henon_heiles_sv(const vec& y, const double& h, vec Y_vec, vec& q_next);
mat stormer_verlet(const double& t_0, const double& t_end, const vec& y0, const double& h);