#ifndef flick_material_mixture
#define flick_material_mixture

#include "z_profile.hpp"
#include "../environment/configuration.hpp"
#include <map>

namespace flick {
namespace material {
  template<class Function>
  struct mixture : public z_profile<Function> {
    using z_profile<Function>::a_profile_;
    using z_profile<Function>::s_profile_;
    using z_profile<Function>::real_refractive_index_;
    struct configuration : basic_configuration {
      configuration() {
	add<size_t>("n_angles", 300, R"(Number of grid points used to sample the volume scattering function)");
	add<size_t>("n_heights", 8, R"(Number of grid points used sample vertical atmospheric gas profiles)");
      }
    };
  private:
    stdvector angles_;
    stdvector heights_;
    bool auto_update_iops_ = true;
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
      auto m = std::make_shared<Material>(a...);
      add_material(m, id<Material>());

      /*
      std::string key = id<Material>();
      if (exists(key))
      	throw std::runtime_error("material add: " + key + " already exists");
      materials_[key] = std::make_shared<Material>(a...);
      if (range_.find(key) == range_.end())
      	range_[key]={0, heights_.size()-1};
      if (auto_update_iops_)
	update_iops();
      */
    }
  
    void add_material(const std::shared_ptr<base>& m, const std::string& unique_name) {
      std::string key = unique_name;
      if (exists(key))
      	throw std::runtime_error("material add: " + key + " already exists");
      materials_[key] = m;
      if (range_.find(key) == range_.end())
      	range_[key]={0, heights_.size()-1};
      if (auto_update_iops_)
	update_iops();
    }
    /*
    template<class Material>
    void add_material(const std::shared_ptr<Material>& m, const std::string& unique_name="") {
      std::string key = id<Material>();
      if (not unique_name.empty())
	key = unique_name;
      if (exists(key))
      	throw std::runtime_error("material add: " + key + " already exists");
      materials_[key] = m;
      if (range_.find(key) == range_.end())
      	range_[key]={0, heights_.size()-1};
      if (auto_update_iops_)
	update_iops();
    }
    */
    
    template<class Material>
    const Material& get_material() const {
      auto ptr = materials_.at(id<Material>()).get();
      return *static_cast<Material*>(ptr);
    }
    
    template<class Material>
    void set_range(size_t n_low, size_t n_high) {
      set_range(n_low, n_high, id<Material>());
      /*
      range_[id<Material>()] = {n_low, n_high};
      if (auto_update_iops_)
	update_iops();
      */
    }
    
    void set_range(size_t n_low, size_t n_high, const std::string& name) {
      if (not exists(name))
      	throw std::runtime_error("material set_range: " + name + " does not exists");
      range_[name] = {n_low, n_high};
      if (auto_update_iops_)
	update_iops();
    }
    
    stdvector angle_range(size_t n) const {
      stdvector a = range(0,constants::pi,n+1).linspace();
      a.erase(a.begin());
      return a;
    }
    void auto_update_iops(bool b) {
      auto_update_iops_ = b;
      if (auto_update_iops_)
	update_iops();	
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      mueller m;
      double theta = z_profile<Function>::angle(scattering_direction);
      size_t n = s_profile_.low_index_near(z_profile<Function>::pose().position().z());
      for (size_t i=0; i<mueller_[n].size(); ++i) {
	size_t row = mueller_[n](i).row;
	size_t col = mueller_[n](i).col;
	stdvector x{heights_[n],heights_[n+1]};
	stdvector y{mueller_[n].value(row,col,theta),
	  mueller_[n+1].value(row,col,theta)};
	pe_function f{x,y};
	m.add(row,col,f.value(z_profile<Function>::pose().position().z()));
      }
      return m;
    }
    void set_wavelength(double wl) override {
      for (auto const& [key, material] : materials_) {
	material->set_wavelength(wl);
      }
      base::set_wavelength(wl);
      if (auto_update_iops_)
	update_iops();
    }
    void update_iops() {
      a_profile_.clear();
      s_profile_.clear();
      real_refractive_index_ = 1;
      mueller_.clear();
      mueller_.resize(heights_.size());
      for (auto const& [key, material] : materials_) {
	std::vector<size_t> range = range_.at(key);
	add(*material, range[0], range[1]);
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
      add_mueller(material,n_low,n_high);
      add_absorption_and_scattering(material, n_low, n_high);
      if (material.real_refractive_index() > real_refractive_index_) {
	real_refractive_index_ = material.real_refractive_index();
      }
      material.set(initial_pose);
    }
    void add_absorption_and_scattering(base& material, size_t n_low, size_t n_high) {
      stdvector zs = {heights_.begin()+n_low, heights_.begin()+n_high+1};
      Function a;
      Function s;
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
      a_profile_.add(iop_z_profile<Function>(a), heights_);
      s_profile_.add(iop_z_profile<Function>(s), heights_);
    }
    void add_mueller(base& material, size_t n_low, size_t n_high)
    // Must be run before adding scattering coefficients
    {
      for (size_t i=n_low; i<=n_high; ++i) {
	double z = heights_[i];
	double s1 = 0;
	double weight = 0;
	if (s_profile_.size() > 0)
	  s1 = s_profile_.value(z);
	material.set_position({0,0,z});
	double s2 = material.scattering_coefficient();
	if (s1+s2 > 0)
	  weight = s2/(s1+s2);
	if (not std::isfinite(weight))
	  weight = 0;
	angular_mueller am = fill_angular_mueller(material);
	mueller_[i].add(am, weight);
      }
    }
    void add(base& material) {
      add(material, 0, heights_.size()-1);
    }
    Function add_extrapolation_points(Function f) const {
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
