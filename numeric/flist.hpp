#ifndef flick_flist
#define flick_flist
 
#include "function.hpp"
#include "linalg/matrix.hpp"

namespace flick {
  template<class Interpolation>
  class flist {
    std::string header_;
    std::vector<function<Interpolation>> functions_;
  public:
    std::string header() const {
      return header_;
    }
    const function<Interpolation>& operator()(size_t n) const {
      return functions_.at(n);
    }
  private:
    friend std::ostream& operator<<(std::ostream &os, const flist<Interpolation>& fl) {
      os << fl.header_;
      for (size_t i=0; i<fl.functions_.size(); i++) {
	os << "Function " << i << ":\n";
	os << fl(i) << std::endl;
      }
      return os;
    }
    friend std::istream& operator>>(std::istream &is, flist<Interpolation>& fl) {
      using namespace linalg;
      fl.functions_.clear();
      fl.header_ = read_header(is);
      matrix m;
      is >> m;
      m = t(m);
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
