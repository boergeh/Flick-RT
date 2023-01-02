#ifndef flick_coating
#define flick_coating

#include <complex>
#include "../environment/input_output.hpp"
#include "../numeric/constants.hpp"
#include "../numeric/named_bounded_types.hpp"
#include "../numeric/rotation.hpp"
#include "../polarization/mueller.hpp"
#include "../polarization/fresnel.hpp"

namespace flick {
namespace coating {
  /*
    A coating may be added to the outside of a volume surface.
  */
  class base {
  protected:
    wavelength wavelength_{500e-9};
    rotation incidence_{quaternion{1,0,0,0}};
    std::complex<double> relative_refractive_index_{1,0};
    double ui_polar_{0};
    double ui_azimuth_{0};
    unit_vector facing_surface_normal_{0,0,1};
    const double pi = constants::pi;
  public:
    virtual mueller reflection_mueller_matrix() = 0;
    virtual mueller transmission_mueller_matrix() = 0;
    virtual double unpolarized_reflectance() = 0;
    virtual double unpolarized_transmittance() = 0;
    virtual rotation reflection_rotation() = 0;
    virtual rotation transmission_rotation() = 0;
    
    void set_direction_parameters(const unit_interval& ui_polar,
				  const unit_interval& ui_azimuth) {
      ui_polar_ = ui_polar();
      ui_azimuth_ = ui_azimuth();
    }
    void set_incidence(const rotation& r) {
      incidence_ = r;
    }
    void set(const unit_vector& facing_surface_normal) {
      facing_surface_normal_ = facing_surface_normal;
    }
    void set(const wavelength& wl) {
      wavelength_ = wl;
    }
    void set(const std::complex<double>& relative_refractive_index) {
      relative_refractive_index_ = relative_refractive_index;
    }
    friend std::ostream& operator<<(std::ostream &os, base& b) {
      os << " at wavelength " <<  b.wavelength_
	 << ", polar_ui " << b.ui_polar_ << ", and azimuth_ui "
	 << b.ui_azimuth_ <<": ";
      os << "reflectance "<< b.unpolarized_reflectance();
      os << ", transmittance "<< b.unpolarized_transmittance();
      return os;
    }
  };
  
  class grey_lambert : public base {
    double reflectivity_{1};
    double transmissivity_{0};
  public:
    grey_lambert(double reflectivity, double transmissivity)
      : reflectivity_{reflectivity}, transmissivity_{transmissivity}
    {}
    mueller reflection_mueller_matrix() {
      mueller m;
      m.add(0,0, reflectivity_);
      return m;
    }
    mueller transmission_mueller_matrix() {
      mueller m;
      m.add(0,0, transmissivity_);
      return m;
    }
    double unpolarized_reflectance() {
      return reflectivity_;
    }
    double unpolarized_transmittance() {
      return transmissivity_;
    }
    rotation reflection_rotation() {
      double theta = asin(sqrt(ui_polar_));
      double phi = 2*pi*ui_azimuth_;
      rotation r{facing_surface_normal_};
      r.rotate_about_local_z(phi);
      r.rotate_about_local_x(theta);
      return r;
    }
    rotation transmission_rotation() {
      rotation r{reflection_rotation()};
      r.rotate_about_local_x(pi);
      return r;
    }
  };

  class white_lambert : public grey_lambert {
  public:
    white_lambert() : grey_lambert::grey_lambert(1,0) {}
  };

  class fresnel : public base {   
    flick::fresnel f_;
  public:
    mueller reflection_mueller_matrix() {
      update_fresnel();
      return reflection_mueller(f_);
    }  
    mueller transmission_mueller_matrix() {
      update_fresnel();
      return transmission_mueller(f_);
    }
    double unpolarized_reflectance() {
      update_fresnel();
      return f_.R();
    }
    double unpolarized_transmittance() {
      update_fresnel();
      return f_.T();
    }
    rotation reflection_rotation() {
      update_fresnel();
      rotation r{incidence_};
      r.rotate_about(facing_surface_normal_,pi);
      r.rotate_about_local_x(pi);
      return r;
    }
    rotation transmission_rotation() {
      update_fresnel();
      rotation r{-facing_surface_normal_};
      assert(fabs(dot(r.x_direction(),facing_surface_normal_))<1e-6);
      r.rotate_about_local_x(real(f_.transmission_angle()));
      return r;
    }
  private:
    void update_fresnel() {
      double cos_theta = dot(-incidence_.z_direction(),facing_surface_normal_);
      f_ = flick::fresnel(relative_refractive_index_, cos_theta);
    } 
  };
}
}

#endif
