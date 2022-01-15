#include "rng.hpp"

using namespace std;

mt19937 rng;
uniform_real_distribution<double> random01(0.0,1.0);
uniform_int_distribution<> randomInt;

namespace lehmer {
  double random(double *seed, double coef) {
    double rd, rf;
    rd = 16807 * (*seed);
    rf = floor(rd / coef);
    *seed = rd - rf * coef;
    return (*seed / (coef + 1));
  }

  double random(double *seed) {
    double coef = 2048;
    coef *= 1024;
    coef *= 1024;
    coef -= 1;
    return random(seed,coef);
  }
}

