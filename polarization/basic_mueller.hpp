#ifndef flick_basic_mueller
#define flick_basic_mueller

#include "../numeric/vector.hpp"

namespace flick {
  class basic_mueller
  // See e.g. Wikipedia Mueller calculus
  {
  protected:
    std::vector<std::vector<double>> m_{4, std::vector<double>(4)};
  public:
    double element(size_t row, size_t col) const {
      return m_.at(row).at(col);
    }
    void element(size_t row, size_t col, double value) {
      m_.at(row).at(col) = value;
    }
    friend std::ostream& operator<<(std::ostream &os, const basic_mueller& bm) {
      for (size_t i = 0; i < 4; ++i) {
	for (size_t j = 0; j < 4; ++j) {
	  os << bm.m_[i][j] << " ";
	}
	os << '\n';
      }
      return os;
    }
  };
}

#endif
