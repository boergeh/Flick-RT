#ifndef flick_uniform_random
#define flick_uniform_random

#include <random>
#include <chrono>

namespace flick {
  class uniform_random {
    mutable std::default_random_engine generator;
    mutable std::uniform_real_distribution<double> distribution;
  public:
    uniform_random()
    : distribution{std::uniform_real_distribution<double>(0,1)} {
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      generator = std::default_random_engine(seed);
    }
    double operator()() const {
      return distribution(generator);
    }
    double operator()(double min, double max) const {
      return min + distribution(generator)*(max-min);
    }
  };
}
#endif
