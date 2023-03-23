#ifndef flick_range
#define flick_range

#include <vector>
#include <valarray>

namespace flick {
  
  template<class T>
  class basic_range {
    double low;
    double high;
    size_t n_points;
  public:
    basic_range(double l, double h, size_t n)
      : low{l}, high{h}, n_points{n} {}
    T linspace() const {
      double step = (high-low)/(n_points-1);
      if (n_points < 2 || high < low) {
	T v;
	v.resize(1);
	v[0] = low;
	return v;
      }
      T r(n_points);
      r[0] = low;
      for (size_t i = 1; i < r.size(); ++i) {
	r[i] = r[i-1]+step;
      }
      return r;      
    }
    T logspace() const {
      return exp(basic_range<T>(log(low),log(high),n_points).linspace());  
    }
  private:
    T exp(T v) const {
      T v_out(v.size());
      for (size_t i = 0; i < v.size(); ++i) {
	v_out[i] = std::exp(v[i]);
      }
      return v_out;
    }
  };

  using range = basic_range<std::vector<double>>;
  using valarray_range = basic_range<std::valarray<double>>;

  std::ostream& operator<<(std::ostream& os, const std::valarray<double>& va) {
    for (size_t i = 0; i < va.size(); ++i)
      os << va[i] << " ";
    return os;
  }

  /*
  std::ostream& operator<<(std::ostream& os, const std::vector<double>& v) {
    for (size_t i = 0; i < v.size(); ++i)
      os << v[i] << " ";
    return os;
  }
  */
}
#endif
