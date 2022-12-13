#ifndef flick_material_atmosphere
#define flick_material_atmosphere

#include "../material.hpp"
#include "../iop_profile.hpp"
#include "../../polarization/rayleigh_mueller.hpp"
#include "atmosphere_state.hpp"
#include "cross_section.hpp"
#include "../z_profile.hpp"

namespace flick {
namespace material {  
  class atmosphere : public z_profile {
    std::vector<lines> l_;
    o3_cross_section o3_;
  public:
    atmosphere(const atmosphere_state& atm) {
      auto gas_names = atm.gas_names();
      for (size_t i=0; i<gas_names.size(); ++i) {
	l_.push_back(lines(gas_names.at(i)));
      }
      std::vector<double> h = atm.height_grid();
      std::vector<double> a(h.size(),0.0);
      std::vector<double> s(h.size(),0.0);
      for (size_t i=0; i<h.size(); ++i) {
	for (size_t j=0; j<atm.gas_names().size(); ++j) {
	  if (gas_names.at(j)=="o3" && wavelength() < o3_.longest())
	    a[i] += o3_.value(wavelength(),atm.temperature(h[i]));
	  else {
	    double partial_pressure =
	      atm.gas_concentration(gas_names.at(j),h[i]) /
	      atm.air_concentration(h[i]);
	    l_[j].wing_cutoff(200);
	    l_[j].temperature(atm.temperature(h[i]));
	    l_[j].total_pressure(atm.pressure(h[i]));
	    l_[j].partial_pressure(partial_pressure);
	    a[i] += l_[j].absorption_coefficient(wavelength());
	  }
	}
	s[i] = scattering_cross_section(wavelength()) *
	  atm.air_concentration(h[i]);
      }
      a_profile_ = iop_z_profile(pe_function(h,a));
      s_profile_ = iop_z_profile(pe_function(h,s));
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) {
      return rayleigh_mueller(angle(scattering_direction),0.0279);
    }
    double real_refractive_index() {
      return 1;
    }
  private:
    double scattering_cross_section(double wavelength) const
    /* Thomas, G.E. and Stamnes, K., 2002. Radiative transfer in the
       atmosphere and ocean. Cambridge University Press. */
    {
      double to_m2 = 1e-32;
      double wl_mu = wavelength * 1e6;
      double a[] = {3.9729066, 4.6547659e-2, 4.5055995e-4, 2.3229848e-5};
      double s = 0;
      for (int i = 0; i < 4; i++) 
	s += pow(wl_mu, -4) * a[i] * pow(wl_mu, -2*i) * to_m2;
      return s;
    }
  };
}
}

#endif
