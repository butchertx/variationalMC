#pragma once
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif // !_USE_MATH_DEFINES
#include "mkl.h"
#include <vector>

std::vector<double> solve_linear_system(std::vector<std::vector<double>> A, std::vector<double> b);