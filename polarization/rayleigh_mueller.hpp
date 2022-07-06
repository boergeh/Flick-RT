#ifndef flick_rayleigh_mueller
#define flick_rayleigh_mueller

#include "basic_mueller.hpp"

namespace flick {
  class rayleigh_mueller : public basic_mueller
  // See e.g. book by Knut and Jakob Stamnes, 2015, page 21, and book
  // by Bohren and Huffman, 1983, page 156.
  {
  public:
    rayleigh_mueller(double angle, double depolarization=0)
    {
      double c = cos(angle);
      double s = sin(angle);
      double f = (1-depolarization)/(1+depolarization);
      double k = 3/(3+f)*4./3*3/(16*constants::pi);
      m_[0][0] = k*(1+f*c*c);
      m_[0][1] = -k*f*c*c;
      m_[1][0] = -k*f*s*s;
      m_[1][1] = k*f*(1+c*c);
      m_[2][2] = k*2*f*c;
      m_[3][3] = k*(3*f-1)*c;
    }
  };
}
#endif
