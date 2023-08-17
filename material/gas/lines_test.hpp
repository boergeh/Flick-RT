#include "lines.hpp"

namespace flick {
  begin_test_case(lines_test) {  
    // Test against Bohren and Clothiaux, chapter 2, figure 2.19, with
    // continuum
    auto h = lines("h2o");
    h.wing_cutoff(200);
    double hPa = 1e2;
    std::vector<double> pressure = {1013*hPa, 10*hPa, 1013*hPa, 10*hPa};
    std::vector<double> wn = {1085.7, 1085.7, 1084.7, 1084.7};
    std::vector<double> bench_mu2 = {2e-16, 4e-19, 1.8e-16, 1e-19};
    for (size_t i=0; i<wn.size(); ++i) {
      double T = 296;
      h.temperature(T);
      h.total_pressure(pressure.at(i));
      double p_self = 1e-3*pressure.at(i);
      h.partial_pressure(p_self);
      double wavelength = 1/wn.at(i)*1e-2;
      const double mu2_per_m2 = 1e12;
      double cross_section_mu2 = h.absorption_coefficient(wavelength) /
	lines::molecules_per_volume(p_self,T)*mu2_per_m2;
      double ratio = cross_section_mu2/bench_mu2.at(i); 
      check(ratio > 0.1 && ratio < 3);
    }
  } end_test_case()
}
