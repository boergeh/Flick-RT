#ifndef flick_material_z_profile
#define flick_material_z_profile

#include "material.hpp"
#include "iop_profile.hpp"
#include "../numeric/function.hpp"

namespace flick {
namespace material { 
  class z_profile : public base {
  protected:
    iop_z_profile a_profile_;
    iop_z_profile s_profile_;
    double real_refractive_index_{1};
  public:
    z_profile() = default;
    z_profile(base& mat, const stdvector& z) {
      stdvector a(z.size());
      stdvector s(z.size());
      for (size_t i=0; i<z.size(); i++) {
	vector r = mat.pose().position();
	r = {r.x(), r.y(), z[i]};
	mat.set(mat.pose().move_to(r));
	a[i] = mat.absorption_coefficient();
	s[i] = mat.scattering_coefficient();
      }
      a_profile_ = iop_z_profile(pe_function(z,a));
      s_profile_ = iop_z_profile(pe_function(z,s));
    }
    const iop_z_profile& a_profile() const {
      return a_profile_;
    }
    const iop_z_profile& s_profile() const {
      return s_profile_;
    }
    const stdvector& height_grid() const {
      return a_profile_.height_grid();
    }
    virtual mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double theta = angle(scattering_direction);
      double g = 0;
      m.add(0,0,flick::henyey_greenstein(g).phase_function(theta));
      return m;
    }
    virtual double real_refractive_index() const {
      return real_refractive_index_;
    } 
    double absorption_coefficient() const {
      return a_profile_.value(pose().position().z());
    }
    double scattering_coefficient() const {
      return s_profile_.value(pose().position().z());
    }
    double absorption_optical_depth(double distance) const {
      return a_profile_.optical_depth(pose(),distance);
    }
    double scattering_optical_depth(double distance) const {
      return s_profile_.optical_depth(pose(),distance);
    }
    double absorption_distance(double absorption_optical_depth) const {
      return s_profile_.distance(pose(),absorption_optical_depth);
    }
    double scattering_distance(double scattering_optical_depth) const {
      return s_profile_.distance(pose(),scattering_optical_depth);
    }
  private:
    friend std::ostream& operator<<(std::ostream &os, const z_profile& p) {
      os << "\n";
      os << "Absorption coefficient profile" << "\n";
      os << p.a_profile_;
      os << "Scattering coefficient profile" << "\n";
      os << p.s_profile_;
      return os;
    }
  };
}
}

#endif
