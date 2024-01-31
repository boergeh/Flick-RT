#ifndef flick_flist
#define flick_flist
 
#include "function.hpp"
#include "linalg/matrix.hpp"

namespace flick {
  template<class Interpolation>
  class flist {
    std::string header_;
    std::vector<function<Interpolation>> functions_;
    linalg::matrix all_values_;
  public:
    std::string header() const {
      return header_;
    }
    const function<Interpolation>& operator()(size_t n) const {
      return functions_.at(n);
    }
    const linalg::matrix& matrix() {
      return all_values_;
    }
  private:
    friend std::ostream& operator<<(std::ostream& os, const flist<Interpolation>& fl) {
      using namespace linalg;
      os << fl.header_;
      os << fl.all_values_ << std::endl;
      return os;
    }
    friend std::istream& operator>>(std::istream& is, flist<Interpolation>& fl) {
      using namespace linalg;
      fl.functions_.clear();
      fl.header_ = read_header(is);
      is >> fl.all_values_;
      linalg::matrix m = t(fl.all_values_);
      for (size_t i=1; i<m.size(); i++) {
	fl.functions_.push_back(function<Interpolation>(m.at(0),m.at(i)));
      }
      return is;
    }
  };

  using pl_flist = flist<piecewise_linear>;
  using pp_flist = flist<piecewise_power>;
  using pe_flist = flist<piecewise_exponential>;
}

#endif
