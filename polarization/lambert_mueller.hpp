#ifndef flick_lambert_mueller
#define flick_lambert_mueller

#include "mueller.hpp"

namespace flick {
  mueller lambert_mueller()
  // Lambert surface with isotropic reflectance and transmittance
  {
    mueller m;
    m(0,0, 1/(2*constants::pi));
    return m;
  };
}

#endif
