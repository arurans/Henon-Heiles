#pragma once

#include "../file_names.h"
#include "hamiltonian.h"
#include "poincare.h"

//Compute the hamiltonians/poincaré maps

void compute_hamiltonians(const double& t_0, const double& t_end, const vec& y0, const double& h);
void compute_poincare_maps(const double& t_0, const double& t_end, const vec& y0, const double& h);
void compute_both(const double& t_0, const double& t_end, const vec& y0, const double& h);