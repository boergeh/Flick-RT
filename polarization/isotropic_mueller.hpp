#ifndef flick_isotropic_mueller
#define flick_isotropic_mueller

#include "basic_mueller.hpp"

namespace flick {
  class isotropic_mueller : public basic_mueller
  {
  public:
    isotropic_mueller() : basic_mueller()
    {
      m_[0][0] = 1/(4*constants::pi);
    }
  };
}
#endif
