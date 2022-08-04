#pragma once

// #include "../methods/kahans.h"
#include "../methods/rk4.h"
// #include "../methods/sb.h"
// #include "../methods/sv.h"

//Compute the hamiltonian of a Hénon Heiles system

Array<double, Dynamic, 1> hamiltonian(const Ref<const Matrix<double, 4, Dynamic>> Y);