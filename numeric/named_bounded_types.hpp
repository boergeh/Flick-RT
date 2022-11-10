#ifndef flick_named_bounded_types
#define flick_named_bounded_types

#include "bounded_type.hpp"
#include <complex>

namespace flick {
  using unit_interval = bounded_type<double, zero, one>;
  using wavelength = bounded_type<double, zero, std::exa>;
  using scattering_coefficient = bounded_type<double, zero, std::exa>;
  using absorption_coefficient = bounded_type<double, zero, std::exa>;
  using asymmetry_factor = bounded_type<double, negative_one, one>;

  using zenith_angle = bounded_type<double, zero, pi_half>;
  using polar_angle = bounded_type<double, zero, one_pi>;
  using azimuth_angle = bounded_type<double, zero, two_pi>;
  using angle = bounded_type<double, zero, two_pi>;
  using vertex_angle = bounded_type<double, zero, one_pi>;

  using number = bounded_type<size_t, zero, std::exa>;
  using percentage = bounded_type<double, zero, hundred>;
  using thickness = bounded_type<double, zero, std::exa>;
  
  using albedo = unit_interval;
  using sampling_asymmetry_factor = asymmetry_factor;
  using relative_refractive_index = std::complex<double>;
  using real_refractive_index = bounded_type<double, zero, std::exa>;
}

#endif
