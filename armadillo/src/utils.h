#pragma once

#include <armadillo>
#include <cmath>
#include <iostream>
#include <fstream>
//#include <string>
#include <utility>

using arma::csv_ascii;
using arma::eye;
using arma::mat;
using arma::regspace;
using arma::solve;
using arma::uword;
using arma::vec;
using arma::zeros;

std::pair<uword,double> create_H(const double& t_0, const double& t_end, const double& h);
vec create_init_cond(const double& H_0);
vec create_T(const double& t_0, const double& t_end, const double& h);
std::string decimal_to_string(double h);