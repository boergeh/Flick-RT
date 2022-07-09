#ifndef flick_polarization_algorithm
#define flick_polarization_algorithm

#include "stokes.hpp"
#include "mueller.hpp"

namespace flick {
  stokes operator*(const mueller& m, const stokes& s)
  // sparse matrix and vector multiplication
  {
    std::vector<double> s1 = {s.I(),s.Q(),s.U(),s.V()};
    std::vector<double> s2 = {0,0,0,0};
    for (size_t n = 0; n < m.size(); ++n) {
      s2[m(n).row] += m(n).value * s1[m(n).col];
    }
    return stokes(s2);
  }
}

#endif
