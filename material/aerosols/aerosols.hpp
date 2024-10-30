#ifndef flick_material_aerosols
#define flick_material_aerosols

#include "../../numeric/function.hpp"
#include "../../numeric/table.hpp"
#include "../../environment/input_output.hpp"
#include "../iop_profile.hpp"
#include "../z_profile.hpp"

namespace flick {
namespace material {
  enum aerosol_name{urban,rural};
  template<aerosol_name Name>
  class aerosols : public z_profile<pe_function>
  /* Implementation adapted from: Ruiz-Arias, J.A., Dudhia, J. and
     Gueymard, C.A., 2014. A simple parameterization of the short-wave
     aerosol optical properties for surface direct and diffuse
     irradiances assessment in a numerical weather
     model. Geoscientific Model Development, 7(3), pp.1159-1174. */
  {
    std::vector<pl_function> alpha_functions_;
    pl_table omega_;
    pl_table asymmetry_factor_;
    double rh_;
    double od550_;
    double layer_thickness_;
    const std::vector<double> humidity_{0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99};
  public:
    aerosols(double layer_thickness,
	     double optical_depth_at_550nm, double relative_humidity)
      : layer_thickness_{layer_thickness}, od550_{optical_depth_at_550nm},
	rh_{relative_humidity} {
      std::string name;
      if(Name==rural) {
	name = "rural";
	make_alpha({{1.036, 1.035, 1.030, 0.999, 0.946, 0.906, 0.818, 0.753},
		    {1.433, 1.430, 1.421, 1.382, 1.371, 1.357, 1.221, 1.152}});
      } else if (Name==urban) {
	name = "urban";
	make_alpha({{0.915, 0.919, 0.929, 0.921, 0.875, 0.803, 0.682, 0.588},
		    {1.198, 1.202, 1.202, 1.254, 1.265, 1.243, 1.164, 1.082}});
      }
      const std::string path = "material/aerosols/";
      omega_ = read<pl_table>(path+name+"_ssalbedo.txt");
      asymmetry_factor_ = read<pl_table>(path+name+"_asymmetry.txt");
      make_iop_profiles();
    }
    void set_wavelength(double wl) override {
      base::set_wavelength(wl);
      make_iop_profiles();
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      mueller m;
      double theta = angle(scattering_direction);
      double frequency = constants::c/wavelength();
      double g = asymmetry_factor_.value(rh_,frequency);
      m.add(0,0,flick::henyey_greenstein(g).phase_function(theta));
      return m;
    }    
    double real_refractive_index() const override {
      return 1;
    }
    double optical_thickness() const {
      size_t i = 0;
      if (wavelength() > 550e-9)
	i = 1;
      double alpha = alpha_functions_[i].value(rh_);
      return od550_ * pow(wavelength()/550e-9, -alpha);	
    }
  private:
    double height_profile() const {
      double s = 1e-5; // Gives almost constant concentration with height
      double h = layer_thickness_;
      double z = pose().position().z();
      double k = s/(1-exp(-s*h));
      return k * exp(-s*z);
    }
    double get_absorption_coefficient() const {
      return get_attenuation_coefficient() - get_scattering_coefficient();
    }
    double get_scattering_coefficient() const {
      double frequency = constants::c/wavelength();
      return omega_.value(rh_,frequency) * get_attenuation_coefficient();
    }
    void make_iop_profiles() {
      stdvector h = {0, layer_thickness_};
      stdvector a(2), s(2);
      base::set_position({0,0,0});
      a[0] = get_absorption_coefficient();
      s[0] = get_scattering_coefficient();
      base::set_position({0,0,layer_thickness_});
      a[1] = get_absorption_coefficient();
      s[1] = get_scattering_coefficient();
      base::set_position({0,0,0});
      a_profile_ = iop_z_profile<pe_function>(pe_function(h,a));
      s_profile_ = iop_z_profile<pe_function>(pe_function(h,s));
    }
    void make_alpha(const std::vector<std::vector<double>>& alpha) {
      alpha_functions_.resize(2);
      for (size_t i=0; i<alpha.size(); ++i)
	alpha_functions_[i] = pl_function(humidity_, alpha[i]);
    }
    double get_attenuation_coefficient() const {
      return optical_thickness() * height_profile();
    }
  };
  using rural_aerosols = aerosols<rural>;
  using urban_aerosols = aerosols<urban>;
}
}

#endif
