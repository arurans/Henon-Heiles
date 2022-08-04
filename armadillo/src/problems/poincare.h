#pragma once

#include "../methods/kahans.h"
#include "../methods/rk4.h"
#include "../methods/sb.h"
#include "../methods/sv.h"

//Find the Poincaré map of a Hénon Heiles system

mat poincare(const mat& Y);