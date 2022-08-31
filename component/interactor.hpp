#ifndef flick_interactor
#define flick_interactor
/*
#include "../numeric/direction_generator.hpp"
#include "../numeric/physics_function.hpp"
#include "../polarization/fresnel.hpp"
#include "radiation_package.hpp"

namespace flick {
  class basic_interactor {
  protected:
    radiation_package rp_;
    mueller m_;
    uniform_random rnd_;
    direction_generator dg_;
  public:
    void radiation_package(const radiation_package& rp) {
      rp_=rp;
    }
  }
  *;
  /*
  class surface_interactor : public basic_interactor {
  protected:
    unit_vector normal_{0,0,1};
  public:
    void normal(unit_vector normal) {
      normal_ = normal;
    } 
    unit_vector relative_normal() const {
      double reldir = dot(rp_.pose().z_direction(), normal_);
      if (reldir < 0)
	return normal_;
      return -normal_;
    }
  };
  
  class lambert_surface_interactor : public surface_interactor
  {
    double reflectance_{1};
    double transmittance_{0};
  public:
    lambert_surface_interactor(double reflectance, double transmittance)
      : reflectance_{reflectance}, transmittance_{transmittance} {}
    auto leaving_radiation() const {
      auto rp = rp_;
      if (reflectance_ > rnd_(0,1)) { 
	rp.scale_intensity(reflectance_);
	rp.traveling_direction(dg_.lambertian(relative_normal()));
      } else { 
	rp.scale_intensity(transmittance_);
	rp.traveling_direction(dg_.lambertian(-relative_normal()));
      }
      return rp;
    }
  };
  
  class fresnel_surface_interactor : public surface_interactor {
    fresnel f_;
  public:
    fresnel_surface_interactor(const fresnel& f) : f_{f} {}
    auto leaving_radiation() const {
      auto rp = rp_;
      double interaction_plane_angle=acos(dot(relative_normal(), rp.pose().y_direction()));
      if (f_.R() > rnd_(0,1)) { 
	rp.change_polarization(reflection_mueller(f_),interaction_plane_angle);
	rp.scale_intensity(f_.R());
	rp.rotate_about_local_x(2*(constants::pi/2-f_.reflection_angle()));
      } else {
	rp.change_polarization(transmission_mueller(f_),interaction_plane_angle);
	rp.scale_intensity(f_.T());
	rp.rotate_about_local_x(-f_.reflection_angle()+real(f_.transmission_angle()));
      }
      return rp;
    }
  };

  class isotropic_particle_interactor : public basic_interactor {
    double asymmetry_factor_{0};
  public:
    auto leaving_radiation() const {
      auto rp = rp_;
      double polar_angle = henyey_greenstein{0}.inverted_accumulated_angle(rnd_(0,1));
      double azimuth_angle = rnd_(0,2*constants::pi);
      rp.change_polarization(m_,azimuth_angle);
      rp.rotate_about_local_x(polar_angle);
      return rp;
    }
  };
  */
}

#endif
