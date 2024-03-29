#ifndef flick_bounded_type
#define flick_bounded_type

#include <iostream>
#include <sstream>
#include <ratio>

namespace flick {
  using negative_one = std::ratio<-1>;
  using zero = std::ratio<0>;
  using one = std::ratio<1>;
  using two = std::ratio<2>;
  using three = std::ratio<3>;
  using four = std::ratio<4>;
  using five = std::ratio<5>;
  using fifteen = std::ratio<15>;
  using hundred = std::ratio<100>;
  using one_pi = std::ratio<3141592653589793238,1000000000000000000>;
  using two_pi = std::ratio_multiply<std::ratio<2>,one_pi>;
  using four_pi = std::ratio_multiply<std::ratio<4>,one_pi>;
  using pi_half = std::ratio<3141592653589793238,2000000000000000000>;

  template<class T,  class lower_bound, class upper_bound>
  class bounded_type {
    T value_;
    void check_validity() {
      bool too_low = static_cast<T>(value_ * lower_bound::den)
	< lower_bound::num; 
      bool too_high = static_cast<T>(value_ * upper_bound::den)
	> upper_bound::num;
      if (too_low || too_high) {
	std::stringstream ss;
	ss << "bounded type value " << value_
	   << " is not within the required interval "
	   << double(lower_bound::num) / lower_bound::den
	   << " to " << double(upper_bound::num) / upper_bound::den << " ";
	throw std::invalid_argument(ss.str());
      }
    }
  public:
    bounded_type() = default;
    explicit bounded_type(T value) : value_{value} {
      check_validity();
    }
    T operator()() const {return value_;}
    T value() const {return value_;}
    bounded_type& operator=(T new_value) {
      value_ = new_value;
      check_validity();
      return *this;
    }
    
    friend std::ostream& operator<<(std::ostream &os,
    				    const bounded_type<T,lower_bound,
				    upper_bound>& bt) {
      os << bt();
      return os;
    }
  }; 
}

#endif
