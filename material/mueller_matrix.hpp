#ifndef flick_mueller_matrix
#define flick_mueller_matrix

#include "../numeric/vector.hpp"

namespace flick {
  class mueller_matrix {
  protected:
    std::vector<std::vector<double>> mm_{4, std::vector<double>(4)};
  public:
    double element(size_t row, size_t col) {
      return mm_.at(row).at(col);
    }
    void element(size_t row, size_t col, double value) {
      mm_.at(row).at(col) = value;
    }
  };

  class rayleigh_mueller_matrix : public mueller_matrix
  // See e.g. book by Knut and Jakob Stamnes, 2015, page 21, and book
  // by Bohren and Huffman, 1983, page 156.
  {
  public:
    rayleigh_mueller_matrix(double angle, double depolarization=0)
    {
      double c = cos(angle);
      double s = sin(angle);
      double f = (1-depolarization)/(1+depolarization);
      double k = 3/(3+f)*4./3*3/(16*pi);
      mm_[0][0] = k*(1+f*c*c);
      mm_[0][1] = -k*f*c*c;
      mm_[1][0] = -k*f*s*s;
      mm_[1][1] = k*f*(1+c*c);
      mm_[2][2] = k*2*f*c;
      mm_[3][3] = k*(3*f-1)*c;
    }
  };
}

#endif
