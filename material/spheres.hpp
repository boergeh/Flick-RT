#ifndef flick_material_spheres
#define flick_material_spheres

#include "material.hpp"
#include "water/pure_water.hpp"
#include "../mie/mie.hpp"
#include "../polarization/mueller.hpp"

namespace flick {
namespace material {
  template<class Size_distribution, class Host_material,
	   class Sphere_material, class Monodisperesed_mie>
  class spheres : public base {
    double volume_fraction_;
    Size_distribution size_distribution_;
    Host_material host_material_;
    Sphere_material sphere_material_;

    const std::vector<int> row_{0,0,1,1,2,2,3,3};
    const std::vector<int> col_{0,1,0,1,2,3,2,3};
    stdvector angles_ = range(0,constants::pi,300).linspace();
    int precision_{3};
    using pm = polydispersed_mie<Monodisperesed_mie,Size_distribution>;

    mutable std::shared_ptr<pm> poly_mie_;
    mutable bool has_changed_{true};
    mutable std::vector<pl_function> scattering_matrix_elements_;
    
    
    void update_mie() const {
      if (has_changed_) {
	flick::refractive_index m_host = host_material_.refractive_index();
	flick::refractive_index m_sphere = sphere_material_.refractive_index();	
	Monodisperesed_mie mono_mie(m_host,m_sphere,flick::wavelength{wavelength()});
	mono_mie.precision(precision_);
	mono_mie.angles(angles_);
	poly_mie_ = std::make_shared<pm>(mono_mie, size_distribution_);
	
	for (size_t i=0; i < row_.size(); ++i) {
	  scattering_matrix_elements_[i] =
	    pl_function(angles_,poly_mie_->scattering_matrix_element(row_[i],col_[i]));
	}
      }
      has_changed_ = false;
    }
  public:
    spheres(double volume_fraction,
	    const Size_distribution& sd,
	    const Host_material& hm,
	    const Sphere_material& sm)
      : volume_fraction_{volume_fraction}, size_distribution_{sd},
	host_material_{hm}, sphere_material_{sm} {
      scattering_matrix_elements_.resize(row_.size());
    }
    void set(const flick::wavelength& wl) {
      base::set(wl);
      host_material_.set(wl);
      sphere_material_.set(wl);
      has_changed_ = true;
    }
    void precision(int n) {
      precision_ = n;
      has_changed_ = true;
    }
    void angles(const stdvector& angles) {
      angles_ = angles;
      has_changed_ = true;
    }
    double absorption_coefficient() const {
      update_mie();
      return poly_mie_->absorption_cross_section()
	* size_distribution_.particles_per_volume(volume_fraction_);
    }
    double scattering_coefficient() const {
      update_mie();
      return poly_mie_->scattering_cross_section()
	      * size_distribution_.particles_per_volume(volume_fraction_);
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      update_mie();
      mueller m;
      for (size_t i=0; i<scattering_matrix_elements_.size(); ++i)
	m.add(row_[i],col_[i],
	      scattering_matrix_elements_[i].value(angle(scattering_direction)));
      return m;  
    }
    double real_refractive_index() const {
      return 1;
    }   
  };
  
}
}

#endif
