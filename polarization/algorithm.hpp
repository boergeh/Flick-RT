#ifndef flick_polarization_algorithm
#define flick_polarization_algorithm

#include "stokes.hpp"
#include "basic_mueller.hpp"

namespace flick {
  stokes operator*(const basic_mueller& bm, const stokes& s) {
    std::vector<double> s1 = {s.I(),s.Q(),s.U(),s.V()};
    std::vector<double> s2 = {0,0,0,0};
    for (size_t i = 0; i<4; ++i) {
      for (size_t j = 0; j<4; ++j) {
	s2[i] += bm.element(i,j)*s1[j];
      }
    }
    return stokes(s2);
  }
}

#endif
