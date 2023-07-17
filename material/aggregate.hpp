#ifndef flick_material_aggregate
#define flick_material_aggregate

#include "z_profile.hpp"

namespace flick {
namespace material { 
  class aggregate : public z_profile {
    stdvector angles_;
    stdvector heights_;
    std::vector<angular_mueller> mueller_;
  public:
    aggregate(const stdvector& angles, const stdvector& heights={-1,1})
      : angles_{angles}, heights_{heights} {
      mueller_.resize(heights.size());
    }
    void add(base& material, double z_low, double z_high) {
      add(material,pe_function{{z_low,z_high},{1,1}});
    }
    void add(base& material, pe_function scaling_profile) {
      double z_low = scaling_profile.x()[0];
      double z_high = scaling_profile.x().back();
      flick::pose initial_pose = material.pose();
      stdvector zs = z_low + (heights_ - heights_[0])
	/ (heights_.back()-heights_[0])*(z_high - z_low);
      //zs = insert_extra_points_for_extrapolation(zs);
      pe_function a;
      pe_function s;
      for (auto& z:zs) {
	material.set_position({0,0,z});
	double scaling_factor = scaling_profile.value(z);
	a.append({z,material.absorption_coefficient()*scaling_factor});
	s.append({z,material.scattering_coefficient()*scaling_factor});
      }
      a = add_extrapolation_points(a);
      s = add_extrapolation_points(s);
      material.set_position({0,0,z_low});

      double dz = z_high - z_low;
      a = scale_to_integral(a,material.absorption_optical_depth(dz));
      s = scale_to_integral(s,material.scattering_optical_depth(dz));
      a = multiply(a, scaling_profile, a.x());
      s = multiply(s, scaling_profile, s.x());
      a_profile_.add(iop_z_profile(a), heights_);
      s_profile_.add(iop_z_profile(s), heights_);

      if (material.real_refractive_index() > real_refractive_index_) {
	real_refractive_index_ = material.real_refractive_index();
      }
      for (size_t i=0; i<heights_.size(); ++i) {
	double z = heights_[i];
	double s1 = s_profile_.value(z);
	material.set_position({0,0,z});
	double s2 = material.scattering_coefficient();
	double weight = s2/(s1+s2);
	angular_mueller am = fill_angular_mueller(material);
	mueller_[i].add(am, weight);
      }
      material.set(initial_pose);
    }
    void add(base& material) {
      add(material, heights_[0], heights_.back());
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      mueller m;
      double theta = angle(scattering_direction);
      size_t n = s_profile_.low_index_near(pose().position().z());
      for (size_t i=0; i<mueller_[n].size(); ++i) {
	size_t row = mueller_[n](i).row;
	size_t col = mueller_[n](i).col;
	stdvector x{heights_[n],heights_[n+1]};
	stdvector y{mueller_[n].value(row,col,theta),
	  mueller_[n+1].value(row,col,theta)};
	pe_function f{x,y};
	m.add(row,col,f.value(pose().position().z()));
      }
      return m;
    }
  private:
    stdvector insert_extra_points_for_extrapolation(stdvector& z) const {
      double epsilon = 1e-4;
      double dz_begin = z[1]-z[0];
      double dz_end = z.end()[-1]-z.end()[-2];
      z.insert(z.begin()+1,z[0]+dz_begin*epsilon);
      z.insert(z.end()-1,z.back()-dz_end*epsilon);
      return z;
    }
    pe_function add_extrapolation_points(pe_function& f) const {
      double i0 = f.integral();
      f.add_extrapolation_points(0);
      return scale_to_integral(f,i0);
    }
    angular_mueller fill_angular_mueller(const base& material)
    // Can only do isotropic materials
    {
      pe_function f;
      for (size_t i=0; i<angles_.size(); ++i) {	
	mueller m = material.mueller_matrix(unit_vector{angles_[i],0});
	f.append({angles_[i],m.value(0,0)});
      }
      angular_mueller am{tabulated_phase_function{f}};
      for (size_t i=0; i<angles_.size(); ++i) {	
	mueller m = material.mueller_matrix(unit_vector{angles_[i],0});
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
