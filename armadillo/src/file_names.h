#pragma once

#include <string>

// File with the paths to store computed data

// Path to csv file to store the computed hamiltonians
extern const std::string hamiltonians_file;

//Path to csv file to store the computed Poincaré for Kutta's method
extern const std::string poincare_file_rk4;

//Path to csv file to store the computed Poincaré for Shampine-Bogacki
extern const std::string poincare_file_sb;

//Path to csv file to store the computed Poincaré for Kahans method
extern const std::string poincare_file_kahans;

//Path to csv file to store the computed Poincaré for Störmer-Verlet
extern const std::string poincare_file_sv;