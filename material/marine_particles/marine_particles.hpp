#ifndef flick_material_marine_particles
#define flick_material_marine_particles

#include "../material.hpp"
#include "../../numeric/flist.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class marine_particles : public base {
    pe_function p_;
    pl_function a_star_; // [m^2/g]
    pl_function b_star_; // [m^2/g]
    double mass_concentration_ = 1e-3; // [kg/m^3]
    const double to_mg_ = 1e3; 
    const double to_nm_ = 1e9; 
  public:
    marine_particles(const std::string& name, double mass_concentration=1e-3)
      : mass_concentration_{mass_concentration} {
      const std::string path{"/material/marine_particles/iop_tables"};
      pe_flist PF = read<pe_flist>(path+"/"+name+"_PFs.txt");
      p_ = PF(5);
      pl_flist ab = read<pl_flist>(path+"/"+name+"_ap_bp.txt");
      a_star_ = ab(0);
      b_star_ = ab(2);
    }
    void mass_concentration(double c) {
      mass_concentration_ = c;
    } 
    double absorption_coefficient() const {
      return a_star_.value(wavelength()*to_nm_) * mass_concentration_*to_mg_;
    }
    double scattering_coefficient() const {
      return b_star_.value(wavelength()*to_nm_) * mass_concentration_*to_mg_; 
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double cos_ang = std::clamp<double>(scattering_direction.mu(),-1,1);
      double ang_deg = 180/constants::pi * acos(cos_ang);
      return m.add(0,0,p_.value(ang_deg));
    }
  };
  template<int n>
  struct listable_marine_particles : public marine_particles {
    using marine_particles::marine_particles;
  };
}
}

#endif
