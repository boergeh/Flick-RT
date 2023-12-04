#ifndef flick_material_mixture
#define flick_material_mixture

#include "z_profile.hpp"
#include "../environment/configuration.hpp"
#include <map>

namespace flick {
namespace material {
  class mixture : public z_profile {
    bool should_update_iops_ = true;
    stdvector angles_;
    stdvector heights_;
    std::vector<angular_mueller> mueller_;
    std::map<std::string, std::shared_ptr<material::base>> materials_;
    std::map<std::string, std::vector<size_t>> range_;
    std::string name_extension_;
  public:
    mixture(const stdvector& angles, const stdvector& heights={-1,1})
      : angles_{angles}, heights_{heights} {
      mueller_.resize(heights.size());
    }
    void name_extension(const std::string& s) {
      name_extension_ = s;
    }
    const stdvector& heights() const {
      return heights_;
    }
    template<class Material, class... Args>
    void add_material(Args... a) {
      std::string key = id<Material>();
      if (exists(key))
      	throw std::runtime_error("material add: " + key + " already exists");
      materials_[key] = std::make_shared<Material>(a...);
      if (range_.find(key) == range_.end())
	range_[key]={0, heights_.size()-1};
      update_iops();
    }
    template<class Material>
    const Material& get_material() const {
      auto ptr = materials_.at(id<Material>()).get();
      return *static_cast<Material*>(ptr);
    }
    template<class Material>
    void set_range(size_t n_low, size_t n_high) {
      range_[id<Material>()] = {n_low, n_high};
      update_iops();
    }
    stdvector angle_range(size_t n) const {
      stdvector a = range(0,constants::pi,n+1).linspace();
      a.erase(a.begin());
      return a;
    }
    void should_update_iops(bool tf) {
      should_update_iops_ = tf;
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
    void set_wavelength(double wl) override {
      for (auto const& [key, material] : materials_) {
	material->set_wavelength(wl);
      }
      base::set_wavelength(wl);
      update_iops();
    }
    void update_iops() {
      should_update_iops_=true;
      if (should_update_iops_) {
	a_profile_.clear();
	s_profile_.clear();
	real_refractive_index_ = 1;
	for (auto const& [key, material] : materials_) {
	  std::vector<size_t> range = range_.at(key);
	  add(*material, range[0], range[1]);
	}
      }
    }
    std::vector<std::string> material_ids() const {
      std::vector<std::string> ids;
      for (auto& [name, val] : materials_) {
	ids.push_back(name);
      }
      return ids;
    }
  private:
    bool exists(const std::string& name) const {
      return (materials_.find(name) != materials_.end());
    }
    template<class Material>
    std::string id() const {
      return typeid(Material).name() + name_extension_;
    }
    void add(base& material, size_t n_low, size_t n_high) {
      if (n_low >= n_high or n_high >= heights_.size())
	throw std::runtime_error("mixtrue");
      flick::pose initial_pose = material.pose();
      stdvector zs = {heights_.begin()+n_low, heights_.begin()+n_high+1};
      pe_function a;
      pe_function s;
      for (auto& z:zs) {
	material.set_position({0,0,z});
	a.append({z,material.absorption_coefficient()});
	s.append({z,material.scattering_coefficient()});
      }
      a = add_extrapolation_points(a);
      s = add_extrapolation_points(s);
      material.set_position({0,0,heights_[n_low]});

      double dz = heights_[n_high] - heights_[n_low];
      a = scale_to_integral(a,material.absorption_optical_depth(dz));
      s = scale_to_integral(s,material.scattering_optical_depth(dz));
      a_profile_.add(iop_z_profile(a), heights_);
      s_profile_.add(iop_z_profile(s), heights_);
      if (material.real_refractive_index() > real_refractive_index_) {
	real_refractive_index_ = material.real_refractive_index();
      }
      for (size_t i=n_low; i<=n_high; ++i) {
	double z = heights_[i];
	double s1 = s_profile_.value(z);
	material.set_position({0,0,z});
	double s2 = material.scattering_coefficient();
	double weight = s2/(s1+s2);
	if (not std::isfinite(weight))
	  weight = 0;
	angular_mueller am = fill_angular_mueller(material);
	mueller_[i].add(am, weight);
      }
      material.set(initial_pose);
    }
    void add(base& material) {
      add(material, 0, heights_.size()-1);
    }
    pe_function add_extrapolation_points(pe_function f) const {
      double i0 = f.integral();
      f.add_extrapolation_points(0,1e-5);
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
