#pragma once

#include <random>

extern std::mt19937 rng;
extern std::uniform_real_distribution<double> random01;
extern std::uniform_int_distribution<> randomInt;

namespace lehmer {
  double random(double *seed, double coef);
  double random(double *seed);
}

