#include "storage_info.h"

// Remember that this file resides in the /src folder

// store every skip_storage-th value in the matrix
const int SKIP_STORAGE = 1;

// Path to csv file to store the computed hamiltonians
const std::string hamiltonians_file = "../output/hamiltonians";

//Path to csv file to store the computed Poincaré for Kutta's method
const std::string poincare_file_rk4 = "../output/poincare_rk4";

//Path to csv file to store the computed Poincaré for Shampine-Bogacki
const std::string poincare_file_sb = "../output/poincare_sb";

//Path to csv file to store the computed Poincaré for Kahans method
const std::string poincare_file_kahans = "../output/poincare_kahans";

//Path to csv file to store the computed Poincaré for Störmer-Verlet
const std::string poincare_file_sv = "../output/poincare_sv";