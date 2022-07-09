#ifndef flick_isotropic_mueller
#define flick_isotropic_mueller

#include "mueller.hpp"
#include "../numeric/constants.hpp"

namespace flick {
  mueller isotropic_mueller() {
    mueller m;
    m.add(0,0, 1/(4*constants::pi));
    return m;
  }
}

#endif
