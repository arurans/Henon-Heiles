#pragma once

#include "../methods/kahans.h"
#include "../methods/rk4.h"
#include "../methods/sb.h"
#include "../methods/sv.h"

//Find the Poincaré map of a Hénon Heiles system

Matrix<double, 2, Dynamic> poincare(const Ref<const Matrix<double, 4, Dynamic>> Y);