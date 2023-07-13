#ifndef flick_material_spheres
#define flick_material_spheres

#include "material.hpp"
#include "water/pure_water.hpp"
#include "ice/pure_ice.hpp"
#include "../mie/polydispersed_mie.hpp"
#include "../mie/monodispersed_mie.hpp"
#include "../mie/parameterized_monodispersed_mie.hpp"
#include "../polarization/mueller.hpp"

namespace flick {
namespace material {
  template<class Size_distribution, class Host_material,
	   class Sphere_material, class Monodispersed_mie>
  class spheres : public base {
    double volume_fraction_;
    Size_distribution size_distribution_;
    Host_material host_material_;
    Sphere_material sphere_material_;

    const std::vector<int> row_{0,0,1,1,2,2,3,3};
    const std::vector<int> col_{0,1,0,1,2,3,2,3};
    int precision_{3};
    using pm = polydispersed_mie<Monodispersed_mie,Size_distribution>;

    mutable stdvector tabulated_angles_{0};
    mutable std::shared_ptr<pm> poly_mie_;
    mutable bool has_changed_{true};
    mutable std::vector<pl_function> scattering_matrix_elements_;
    
    void update_mie() const {
      if (has_changed_) {
	stdcomplex m_host = host_material_.refractive_index();
	stdcomplex m_sphere = sphere_material_.refractive_index();	
	Monodispersed_mie mono_mie(m_host,m_sphere,wavelength());
	mono_mie.angles(tabulated_angles_);
	poly_mie_ = std::make_shared<pm>(mono_mie, size_distribution_);
	for (size_t i=0; i < row_.size(); ++i) {
	  scattering_matrix_elements_[i] =
	    pl_function(tabulated_angles_,poly_mie_->
			scattering_matrix_element(row_[i],col_[i]));
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
    void set_wavelength(double wl) {
      base::set_wavelength(wl);
      host_material_.set_wavelength(wl);
      sphere_material_.set_wavelength(wl);
      has_changed_ = true;
    }
    void precision(int n) {
      precision_ = n;
      has_changed_ = true;
    }
    void tabulated_angles(const stdvector& angles) {
      tabulated_angles_ = angles;
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
      double theta = angle(scattering_direction);
      if (tabulated_angles_.size() <= 1)
	has_changed_ = true;
	tabulated_angles_ = stdvector{theta};
      update_mie();
      mueller m;
      for (size_t i=0; i<scattering_matrix_elements_.size(); ++i) {
	double s = scattering_matrix_elements_[i].value(theta);
	m.add(row_[i],col_[i],s/poly_mie_->scattering_cross_section());
      }
      return m;  
    }
    double real_refractive_index() const {
      return 1;
    }   
  };


  template<class Monodispersed_mie>
  struct bubbles_in_ice : public spheres<log_normal_distribution,
					 material::pure_ice,
					 material::vacuum,
					 Monodispersed_mie> {
    bubbles_in_ice(double volume_fraction,double mu, double sigma) :
      spheres<log_normal_distribution,
	      material::pure_ice,
	      material::vacuum,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::pure_ice(),
				 material::vacuum()) {}
  };
  
  template<class Monodispersed_mie>
  struct brines_in_ice : public spheres<log_normal_distribution,
					 material::pure_ice,
					 material::pure_water,
					 Monodispersed_mie> {
    brines_in_ice(double volume_fraction, double mu, double sigma,
		   double salinity) :
      spheres<log_normal_distribution,
	      material::pure_ice,
	      material::pure_water,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::pure_ice(),
				 material::pure_water(salinity,273.15)) {}
  };
  
  template<class Monodispersed_mie>
  struct water_cloud : public spheres<log_normal_distribution,
				      material::vacuum,
				      material::pure_water,
				      Monodispersed_mie> {
    water_cloud(double volume_fraction, double mu, double sigma) :
      spheres<log_normal_distribution,
	      material::vacuum,
	      material::pure_water,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::vacuum(),
				 material::pure_water()) {}
  };
  
  template<class Monodispersed_mie>
  struct ice_cloud : public spheres<log_normal_distribution,
				      material::vacuum,
				      material::pure_ice,
				      Monodispersed_mie> {
    ice_cloud(double volume_fraction, double mu, double sigma) :
      spheres<log_normal_distribution,
	      material::vacuum,
	      material::pure_ice,
	      Monodispersed_mie>(volume_fraction,
				 log_normal_distribution(mu,sigma),
				 material::vacuum(),
				 material::pure_ice()) {}
  };
  
}
}

#endif
