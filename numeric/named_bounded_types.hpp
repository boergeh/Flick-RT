#ifndef flick_named_bounded_types
#define flick_named_bounded_types

#include "bounded_type.hpp"
#include <complex>

namespace flick {
  using unit_interval = bounded_type<double, zero, one>;
  using wavelength = bounded_type<double, zero, std::exa>;
  using scattering_coefficient = bounded_type<double, zero, std::exa>;
  using absorption_coefficient = bounded_type<double, zero, std::exa>;
  //char af_str[]="asymmetry factor";
  //using asymmetry_factor = bounded_type<double, negative_one, one,af_str>;
  using asymmetry_factor = bounded_type<double, negative_one, one>;
  using zenith_angle = bounded_type<double, zero, pi_half>;
  using polar_angle = bounded_type<double, zero, one_pi>;
  using azimuth_angle = bounded_type<double, zero, two_pi>;
  using thickness = bounded_type<double, zero, std::exa>;
  using number_of_packages = bounded_type<size_t, zero, std::exa>;
  using albedo = unit_interval;
  using bottom_albedo = albedo;
  using incidence_angle = zenith_angle;
  using sampling_asymmetry_factor = asymmetry_factor;
  using relative_refractive_index = std::complex<double>;
  //using real_refractive_index = bounded_type<double, zero, std::exa>;
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
