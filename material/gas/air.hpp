#ifndef flick_material_air
#define flick_material_air

#include "../material.hpp"
#include "../iop_profile.hpp"
#include "../../polarization/rayleigh_mueller.hpp"
#include "atmospheric_state.hpp"
#include "cross_section.hpp"
#include "lines.hpp"
#include "../z_profile.hpp"

namespace flick {
namespace material {  
  class air : public z_profile {
    std::vector<lines> l_;
    o3_cross_section o3_;
    atmospheric_state atm_;
  public:
    air(const atmospheric_state& atm)
      : atm_{atm} {
      const auto& h = atm_.height_grid();
      const auto& gas_names = atm_.gas_names();
      for (size_t i=0; i<gas_names.size(); ++i) {
	l_.push_back(lines(gas_names.at(i)));
      }
      for (size_t i=0; i<h.size(); ++i) {
	for (size_t j=0; j<gas_names.size(); ++j) {
	  double partial_pressure =
	    atm_.gas_concentration(gas_names.at(j),h[i]) /
	    atm_.air_concentration(h[i]);
	  l_[j].wing_cutoff(200);
	  l_[j].temperature(atm_.temperature(h[i]));
	  l_[j].total_pressure(atm_.pressure(h[i]));
	  l_[j].partial_pressure(partial_pressure);
	}
      }
      set_iop_profiles();
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      return rayleigh_mueller(angle(scattering_direction),0.0279);
    }
    double real_refractive_index() const override {
      return 1;
    }
    void set_wavelength(double wl) override {
      base::set_wavelength(wl);
      set_iop_profiles();
    }
  private:
    void set_iop_profiles() {
      const auto& h = atm_.height_grid();
      const auto& gas_names = atm_.gas_names();
      std::vector<double> a(h.size(),0.0);
      std::vector<double> s(h.size(),0.0);
      for (size_t i=0; i<h.size(); ++i) {
	for (size_t j=0; j<gas_names.size(); ++j) {
	  if (gas_names.at(j)=="o3" && wavelength() < o3_.longest()) {
	    a[i] += o3_.value(wavelength(),atm_.temperature(h[i]))
	       *  atm_.gas_concentration("o3",h[i]);
	  } else {
	    a[i] += l_[j].absorption_coefficient(wavelength());
	  }
	}
	s[i] = scattering_cross_section(wavelength()) *
	  atm_.air_concentration(h[i]);
      }
      a_profile_ = iop_z_profile(pe_function(h,a));
      s_profile_ = iop_z_profile(pe_function(h,s));
    }
    double scattering_cross_section(double wavelength) const
    /* See Thomas, G.E. and Stamnes, K., 2002, radiative transfer in
       the atmosphere and ocean, Cambridge University Press. */
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
