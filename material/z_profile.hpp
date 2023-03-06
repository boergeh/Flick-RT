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
    double height_{0};
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
    friend std::ostream& operator<<(std::ostream &os,
				    const z_profile& p) {
      os << std::endl;
      os << "Absorption coefficient profile" << std::endl;
      os << p.a_profile_;
      os << "Scattering coefficient profile" << std::endl;
      os << p.s_profile_;
      return os;
    }
  };

  class aggregate_z_profile : public z_profile {
    std::vector<angular_mueller> mueller_;
    std::vector<double> heights_;
    std::vector<double> angles_;
  public:
    aggregate_z_profile(const std::vector<double>& heights,
			const std::vector<double>& angles)
      : heights_{heights}, angles_{angles} {
      mueller_.resize(heights.size());
    }
    void add(const z_profile& p) {
      a_profile_.add(p.a_profile(), heights_);
      s_profile_.add(p.s_profile(), heights_);
      if (p.real_refractive_index() > real_refractive_index_) {
	real_refractive_index_ = p.real_refractive_index();
      }
      for (size_t i=0; i<heights_.size(); ++i) {
	double h = heights_[i];
	double s1 = s_profile_.value(h);
	double s2 = p.s_profile().value(h);
	double weight = s2/(s1+s2);
	angular_mueller am = fill_angular_mueller(p);
	mueller_[i].add(am,weight);
      }
    }    
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      double theta = angle(scattering_direction);
      size_t n = s_profile_.low_index_near(pose().position().z());
      for (size_t i=0; i<mueller_[n].size(); ++i) {
	size_t row = mueller_[n](i).row;
	size_t col = mueller_[n](i).col;
	std::vector<double> x{heights_[n],heights_[n+1]};
	std::vector<double> y{mueller_[n].value(row,col,theta),
	  mueller_[n+1].value(row,col,theta)};
	pe_function f{x,y};
	m.add(row,col,f.value(pose().position().z()));
      }
      return m;
    }
  private:
    angular_mueller fill_angular_mueller(const z_profile& p)
    // Can only do isotropic materials
    {
      pe_function f;
      for (size_t i=0; i<angles_.size(); ++i) {	
	mueller m = p.mueller_matrix(unit_vector{angles_[i],0});
	f.append({angles_[i],m.value(0,0)});
      }
      angular_mueller am{tabulated_phase_function{f}};
      for (size_t i=0; i<angles_.size(); ++i) {	
	mueller m = p.mueller_matrix(unit_vector{angles_[i],0});
	for (size_t j=1; j<m.size(); ++j) {
	  am.append(m(j).row,m(j).col,angles_[i],m.value(m(j).row,m(j).col));
	}
      }
      return am;
    }
  };
}
}

#endif
