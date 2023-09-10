#ifndef flick_material_air
#define flick_material_air

#include "../material.hpp"
#include "../iop_profile.hpp"
#include "../../polarization/rayleigh_mueller.hpp"
#include "../../environment/input_output.hpp"
#include "atmospheric_state.hpp"
#include "cross_section.hpp"
#include "lines.hpp"
#include "../z_profile.hpp"

namespace flick {
namespace material {  
  class basic_air : public z_profile {
  protected:
    atmospheric_state atm_;
  public:
    basic_air(const atmospheric_state& atm)
      : atm_{atm} {}
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      return rayleigh_mueller(angle(scattering_direction),0.0279);
    }
    double real_refractive_index() const override {
      return 1;
    }
    void set_wavelength(double wl) override {
      base::set_wavelength(wl);
      make_iop_profiles();
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
  protected:
    virtual double absorption_coefficient(size_t gas_number, double height) = 0;
    void make_iop_profiles() {
      const auto& h = atm_.height_grid();
      std::vector<double> a(h.size(),0.0);
      std::vector<double> s(h.size(),0.0);
      for (size_t i=0; i<h.size(); ++i) {
	for (size_t j=0; j<atm_.gas_names().size(); ++j) {
	  a[i] += absorption_coefficient(j,h[i]);
	}
	s[i] = scattering_cross_section(wavelength())*atm_.air_concentration(h[i]);
      }
      a_profile_ = iop_z_profile(pe_function(h,a));
      s_profile_ = iop_z_profile(pe_function(h,s));
    }
  };
  
  class hitran_air : public basic_air {
    std::vector<lines> l_;
    o3_cross_section o3_;
    no2_cross_section no2_;
  public:
    hitran_air(const atmospheric_state& atm)
      : basic_air(atm) {
      const auto& gas_names = atm_.gas_names();
      for (size_t i=0; i<gas_names.size(); ++i) {
	l_.push_back(lines(gas_names.at(i)));
      }
      make_iop_profiles();
    }
  private:
    double absorption_coefficient(size_t gas_number, double height) override {
      double h = height;
      size_t j = gas_number;
      std::string gas_name = atm_.gas_names().at(j);
      if (gas_name=="o3" && wavelength() < o3_.longest()) {
	return o3_.value(wavelength(),atm_.temperature(h)) * atm_.gas_concentration("o3",h);
      }
      else if (gas_name=="no2" && wavelength() < no2_.longest()) {
	return no2_.value(wavelength()) * atm_.gas_concentration("no2",h);
      }
      else {
	double partial_pressure = atm_.gas_concentration(gas_name,h) /
	  atm_.air_concentration(h) * atm_.pressure(h);
	l_[j].wing_cutoff(200);
	l_[j].temperature(atm_.temperature(h));
	l_[j].total_pressure(atm_.pressure(h));
	l_[j].partial_pressure(partial_pressure);
	return l_[j].absorption_coefficient(wavelength());
      }
    }
  };

  class smooth_air : public basic_air {
    const double T0_ = constants::T_ntp;
    const double P0_ = constants::P_ntp;
    const double step_ = -0.1;
    std::vector<std::vector<pp_function>> cross_section_;
  public:
    smooth_air(const atmospheric_state& atm, const std::string& spectrum_name)
      : basic_air(atm), cross_section_(3) {
      const auto& gas_names = atm_.gas_names();
      std::vector<std::string> ext = {"T0_P0","T1_P0","T0_P1"};
      for (size_t i=0; i<ext.size(); ++i) {
	for (size_t j=0; j<gas_names.size(); ++j) {
	  std::string file_name = gas_names.at(j)+"_"+spectrum_name+"_"+
	    ext[i]+".txt";
	  cross_section_[i].push_back(read<pp_function>("material/gas/smooth_input/"+file_name));
	}
      }
      make_iop_profiles();
    }
  private:
    double absorption_coefficient(size_t gas_number, double height) override {
      double wl = wavelength();
      double C_00 = cross_section_[0][gas_number].value(wl);
      double C_10 = cross_section_[1][gas_number].value(wl);
      double C_01 = cross_section_[2][gas_number].value(wl);
      double dC_dT = (C_10-C_00)/(T0_*step_);
      double dC_dP = (C_01-C_00)/(P0_*step_);
      double T = atm_.temperature(0); 
      double P = atm_.pressure(0);
      double N = atm_.gas_concentration(atm_.gas_names()[gas_number],height);
      return (C_00 + dC_dT*(T-T0_) + dC_dP*(P-P0_))*N;
    }
  };

  class uv_air : public smooth_air {
  public:
    uv_air(const atmospheric_state& atm) : smooth_air(atm,"uv_air") {}
  };
}
}

#endif
