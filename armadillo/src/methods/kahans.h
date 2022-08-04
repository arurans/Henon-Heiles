#pragma once

//Kahans method of order 2

#include "../utils.h"


void create_A(const vec& y, const double& h, mat& A);
void create_b(const vec& y, const double& h, vec& b);
void kahans_iteration(const vec& y, const double& h, mat& A, vec& b);
mat kahans(const double& t_0, const double& t_end, const vec& y0, const double& h);