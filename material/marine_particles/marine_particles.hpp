#ifndef flick_material_marine_particles
#define flick_material_marine_particles

#include "../material.hpp"
#include "../../numeric/flist.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class marine_particles : public base {
    pe_function p_;
    pl_function a_star_total_; // [m^2/g]
    pl_function a_star_bleached_; // [m^2/g]
    pl_function b_star_; // [m^2/g]
    double mass_concentration_; // [kg/m^3]
    double bleaching_factor_;
    const double to_g_ = 1e3; 
    const double to_nm_ = 1e9;
    const size_t percentile_50th = 2;
  public:
    marine_particles(const std::string& name, double mass_concentration=1e-3,
		     double scattering_scaling_factor=1, double bleaching_factor=0)
      : mass_concentration_{mass_concentration},
	bleaching_factor_{bleaching_factor} {
      std::string p = "/material/marine_particles/iop_tables";
      pe_flist pf = read<pe_flist>(name+"_pf.txt", p);
      p_ = pf(percentile_50th);
      pl_flist ab = read<pl_flist>(name+"_ap_bp.txt", p);
      a_star_total_ = ab(0);
      a_star_bleached_ = ab(1);
      b_star_ = ab(2);
      b_star_.scale_y(scattering_scaling_factor);
      p_.add_constant_extrapolation();
      a_star_total_.add_constant_extrapolation();
      a_star_bleached_.add_constant_extrapolation();
      b_star_.add_constant_extrapolation();
    }
    void mass_concentration(double c) {
      mass_concentration_ = c;
    } 
    double absorption_coefficient() const {
      double wl = wavelength()*to_nm_;
      double a_total = a_star_total_.value(wl) * mass_concentration_ * to_g_;
      double a_bleached = a_star_bleached_.value(wl) * mass_concentration_ * to_g_;
      double a = a_total * (1-bleaching_factor_) + a_bleached * bleaching_factor_;
      if (bleaching_factor_ > 1)
	return a_bleached/bleaching_factor_;
      return a;
    }
    double scattering_coefficient() const {
      bool use_slope = false;
      if (use_slope) {
	double wl = wavelength()*to_nm_;
	double wl0 = 515;
	double slope = -2.0;
	return b_star_.value(wl0)*pow(wl/wl0,slope)*mass_concentration_*to_g_; 
      }
      return b_star_.value(wavelength()*to_nm_) * mass_concentration_*to_g_; 
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double cos_ang = std::clamp<double>(scattering_direction.mu(),-1,1);
      double ang_deg = 180/constants::pi * acos(cos_ang);
      return m.add(0,0,p_.value(ang_deg));
    }
  };
}
}

#endif
