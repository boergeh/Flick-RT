#ifndef flick_material_z_profile
#define flick_material_z_profile

#include "material.hpp"
#include "iop_profile.hpp"

namespace flick {
namespace material { 
  class z_profile : public base {
  protected:
    iop_z_profile a_profile_;
    iop_z_profile s_profile_;
    double real_refractive_index_{1};
  public:
    const iop_z_profile& a_profile() const {
      return a_profile_;
    }
    const iop_z_profile& s_profile() const {
      return s_profile_;
    }
    double real_refractive_index() const {
      return real_refractive_index_;
    } 
    double absorption_coefficient() {
      return a_profile_.value(pose().position().z());
    }
    double scattering_coefficient() {
      return s_profile_.value(pose().position().z());
    }
    double absorption_optical_depth(double distance) {
      return a_profile_.optical_depth(pose(),distance);
    }
    double scattering_optical_depth(double distance) {
      return s_profile_.optical_depth(pose(),distance);
    }
    double absorption_distance(double absorption_optical_depth) {
      return s_profile_.distance(pose(),absorption_optical_depth);
    }
    double scattering_distance(double scattering_optical_depth) {
      return s_profile_.distance(pose(),scattering_optical_depth);
    }
  };

  class aggregate_z_profile : public z_profile {
    std::vector<double> heights_;
    //std::vector<angular_mueller> mueller_;
  public:
    aggregate_z_profile(const std::vector<double>& heights)
      : heights_{heights} {
    }
    void add(const z_profile& p) {
      a_profile_.add(p.a_profile(), heights_);
      s_profile_.add(p.s_profile(), heights_);
      if (p.real_refractive_index() > real_refractive_index_)
	real_refractive_index_ = p.real_refractive_index();
      
    }
    
    mueller mueller_matrix(const unit_vector& scattering_direction) {
      mueller m;
      //double theta = angle(scattering_direction);
      //m.add(0,0,all_angle_mueller.value(0,0,theta));
      return m;
    }
  };
}
}

#endif
