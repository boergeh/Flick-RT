#ifndef flick_named_bounded_types
#define flick_named_bounded_types

#include "../numeric/bounded_type.hpp"

namespace flick {
  using scattering_coefficient = bounded_type<double, zero, std::exa>;
  using absorption_coefficient = bounded_type<double, zero, std::exa>;
  using asymmetry_factor = bounded_type<double, negative_one, one>;
  using albedo = bounded_type<double, zero, one>;
  using zenith_angle = bounded_type<double, zero, pi_half>;
  using thickness = bounded_type<double, zero, std::exa>;
  using number_of_packages = bounded_type<size_t, zero, std::exa>;
  using bottom_albedo = albedo;
  using incidence_angle = zenith_angle;
  /*
  class double_bounded_type {
  protected:
    double bt_=0;
  public:
    double operator()() const {return bt_;}
  };
  
  struct scattering_coefficient : public double_bounded_type {
    scattering_coefficient(double v) {
      bt_ = bounded_type<double, zero, std::exa>(v)();
    }
  };
  
  struct source_zenith_angle : public double_bounded_type {
    source_zenith_angle(double v) {
      bt_ = bounded_type<double, zero, pi_half>(v)();
    }
  };
  */
}

#endif
