#pragma once

#include <eigen3/Eigen/Core>
#include <fstream>
#include <utility>

#include "storage_info.h"

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Matrix;
using Eigen::Ref;

std::tuple<int, int, int, double> create_H(const double& t_0, const double& t_end, const double& h);
Array<double, 4, 1> create_init_cond(const double& H_0);
Array<double, Dynamic, 1> create_T(const double& t_0, const double& t_end, const double& h);
std::string decimal_to_string(double h);
void matrix_to_CSV(std::string filename, const Ref<const Matrix<double, Dynamic, Dynamic>> M);