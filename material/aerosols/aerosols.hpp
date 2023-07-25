#ifndef flick_material_aerosols
#define flick_material_aerosols

#include "../../numeric/function.hpp"
#include "../../numeric/table.hpp"
#include "../../environment/input_output.hpp"
#include "../iop_profile.hpp"
#include "../z_profile.hpp"

namespace flick {
namespace material {
  class aerosols : public z_profile
  /* Implementation adapted from: Ruiz-Arias, J.A., Dudhia, J. and
     Gueymard, C.A., 2014. A simple parameterization of the short-wave
     aerosol optical properties for surface direct and diffuse
     irradiances assessment in a numerical weather
     model. Geoscientific Model Development, 7(3), pp.1159-1174. */
  {
  protected:
    std::vector<double> humidity_{0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99};
    std::vector<std::vector<double>> alpha_;
    std::vector<pl_function> alpha_functions_;
    pl_table asymmetry_factor_;
    double rh_{0};
    double od550_{0.1};

    void make_alpha() {
      alpha_functions_.resize(2);
      for (size_t i=0; i<alpha_.size(); ++i)
	alpha_functions_[i] = pl_function(humidity_, alpha_[i]);
    }
  public:
    aerosols() = default;
    aerosols(double optical_depth_at_550nm, double relative_humidity)
      : od550_{optical_depth_at_550nm}, rh_{relative_humidity} {
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double theta = angle(scattering_direction);
      double frequency = constants::c/wavelength();
      double g = asymmetry_factor_.value(rh_,frequency);
      m.add(0,0,flick::henyey_greenstein(g).phase_function(theta));
      return m;
    }
    double toa_optical_depth() const {
      size_t i = 0;
      if (wavelength() > 550e-9)
	i = 1;
      double alpha = alpha_functions_[i].value(rh_);
      return od550_ * pow(wavelength()/550e-9, -alpha);	
    }
    double real_refractive_index() const {
      return 1;
    }
  protected:
    void make_iop_z_profile(const std::string& name) {
      std::string path = "material/aerosols/";
      pl_table omega = read<pl_table>(path+name+"_ssalbedo.txt");
      pe_function d = read<pe_function>(path+"height_distribution.txt");
      d = d.normalize();
      std::vector<double> h = d.x();
      std::vector<double> a, s;
      for (size_t i=0; i<h.size(); ++i) {
	double att_coef = od550_ * d.value(h[i]);
	double frequency = constants::c/wavelength();
	s.push_back(omega.value(rh_,frequency) * att_coef);
	a.push_back(att_coef - s.back());
      }
      a_profile_ = iop_z_profile(pe_function(h,a));
      s_profile_ = iop_z_profile(pe_function(h,s));
    }
  };
  
  class rural_aerosols : public aerosols {
  public:
    rural_aerosols():rural_aerosols(0.1,0.5) {
    }
    rural_aerosols(double optical_depth_at_550nm, double relative_humidity)
      : aerosols::aerosols(optical_depth_at_550nm, relative_humidity)
    {
      alpha_ = {{1.036, 1.035, 1.030, 0.999, 0.946, 0.906, 0.818, 0.753},
		{1.433, 1.430, 1.421, 1.382, 1.371, 1.357, 1.221, 1.152}};
      make_alpha();
      asymmetry_factor_ = read<pl_table>("material/aerosols/rural_asymmetry.txt");
      make_iop_z_profile("rural");
    }
  };
  
  class urban_aerosols : public aerosols {
  public:
    urban_aerosols() : urban_aerosols(0.1,0.5) {
    }
    urban_aerosols(double optical_depth_at_550nm, double relative_humidity)
      : aerosols::aerosols(optical_depth_at_550nm, relative_humidity)
    {
      alpha_ = {{0.915, 0.919, 0.929, 0.921, 0.875, 0.803, 0.682, 0.588},
		{1.198, 1.202, 1.202, 1.254, 1.265, 1.243, 1.164, 1.082}};
      make_alpha();
      asymmetry_factor_ = read<pl_table>("material/aerosols/urban_asymmetry.txt");
      make_iop_z_profile("urban");
    }    
  };
}
}

#endif
