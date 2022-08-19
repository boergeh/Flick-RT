#ifndef flick_mc_basic
#define flick_mc_basic

#include "../geometry/volume.hpp"
#include "../component/content.hpp"
#include "../component/emitter.hpp"

namespace flick {
  using volume = geometry::volume<content>;

   class wall_interactor {
    geometry::navigator<content>& nav_;
    std::optional<pose> next_wall_intersection_;
    volume* next_volume_;
    const coating::base* coating_;
    radiation_package& rp_;
    uniform_random& rnd_;
    coating::lambert_angle_generator& lag_;
    bool is_reflected_{false};
    bool is_transmitted_{true};
  public:
    wall_interactor(geometry::navigator<content>& nav,
		    radiation_package& rp,
		    uniform_random& ur,
		    coating::lambert_angle_generator& lag)
      : nav_{nav}, rp_{rp}, rnd_{ur}, lag_{lag} {   
      next_wall_intersection_ = nav_.next_intersection(rp_.pose());
      next_volume_ = &nav_.next_volume(rp_.pose());
      volume* v = next_volume_;
      if (is_moving_outward())
      	v = &nav_.current_volume();
      if (v->content().has_coating() && next_wall_intersection_.has_value()) {
	coating_ = &v->content().coating();
	lag_ = coating::lambert_angle_generator{rnd_(0,1),rnd_(0,1)};
	double r = rnd_(0,1);
	is_reflected_ = (r < coating_->reflectivity());
	is_transmitted_ = (1-r < coating_->transmissivity());
      }
    }  
    bool will_intersect() const {
      if (next_wall_intersection_.has_value())
	return true;
      return false;
    }
    void move_to_intersection() const {
      rp_.move_to((*next_wall_intersection_).position());
    }    
    void interact() {
      unit_vector n = facing_surface_normal();
      if (is_reflected_) {
	rp_.move(-nav_.current_volume().small_step());
	change_orientation(n);
	change_radiation();
      } else if (is_transmitted_) {
	if (nav_.current_volume().content().has_coating() &&
	    next_wall_intersection_.has_value()) {
	  change_orientation(-n);
	  change_radiation();
	}
	rp_.move(nav_.current_volume().small_step());
	if (next_wall_intersection_.has_value())
	  store_radiation_if_receivers_are_activated();
	nav_.go_to(*next_volume_);
      } else {
	absorb_radiation_package();
      }
    }
    void store_radiation_if_receivers_are_activated() {
      if (is_moving_inward())
	next_volume_->content().inward_receiver().receive(rp_);
      else if (is_moving_outward())
	nav_.current_volume().content().outward_receiver().receive(rp_);
    }
    auto& current_volume() const {
      return nav_.current_volume();
    }  
    unit_vector facing_surface_normal() const {
      unit_vector n = (*next_wall_intersection_).z_direction();
      if (is_moving_outward())
	n = -n;
      return n;
    }
    void absorb_radiation_package() {
      rp_.scale_intensity(0);
    }
    bool is_moving_inward() const {
      if (!next_wall_intersection_.has_value())
	return false;
      if (nav_.is_moving_inward(*next_wall_intersection_,rp_.pose()))
	return true;
      return false;
    }    
    bool is_moving_outward() const {
      return !is_moving_inward();
    }   
    void change_orientation(const unit_vector& facing_surface_normal) {
      rp_.x_direction_parallel_with_plane_of_incidence(facing_surface_normal);
      rp_.traveling_direction(facing_surface_normal);
      rp_.rotate_about_local_z(lag_.azimuth_angle());
      rp_.rotate_about_local_x(lag_.polar_angle());
    }  
    void change_radiation() {
      double f = 1;
      if (is_reflected_) {
	double brdf = coating_->reflection_mueller_matrix().value(0,0);
	f = brdf/lag_.brdf();
      } else if (is_transmitted_) {
	double brdf = coating_->transmission_mueller_matrix().value(0,0);
	f = brdf/lag_.brdf();
      }
      rp_.scale_intensity(f);
    }
    radiation_package& leaving_radiation_package() {
      return rp_;
    } 
  };
  
  class mc_basic {
    volume* outer_volume_;
    geometry::navigator<content> nav_;
    radiation_package rp_;   
    uniform_random rnd_;
    
    bool lost_in_space() {
      if (!nav_.current_volume().has_outer_volume()
	  && rp_.weighted_traveling_length() > 0)
      	return true;
      return false;
    }
  public:
    mc_basic(geometry::volume<content>& v) {
      outer_volume_ = &v;
      nav_ =  geometry::navigator<content>(*outer_volume_);
    }
    void run(emitter& em, geometry::volume<content>& emitter_volume) {
      uniform_random ur;
      coating::lambert_angle_generator lag;
      while (!em.is_empty()) {
	nav_.go_to(emitter_volume);
	rp_ = em.emit();
	while (!rp_.is_empty() && !lost_in_space()) {
	  wall_interactor wi(nav_,rp_,ur,lag);
	  //particle_interactor pi(nav_,rp_);
	  //if (wi.will_intersect()) {
	    wi.move_to_intersection();	   	  
	    wi.interact();	  	   
	    //rp_ = wi.leaving_radiation_package();
	    //} //else if (pi.will_scatter()) {
	  // pi.move_to_scattering_event()
	  //pi.scatter()
	  // }
	}
      }
    }
  };
  /*
  namespace geometry {
    namespace plane_parallel {
      class basic_slab {
	size_t n_packages_{0};
	double optical_depth_{0};
	double single_scattering_albedo_{0};
	double bottom_albedo_{0};
	double source_zenith_angle_{0};
	semi_infinite_box<content> geometry_;
      public:
	basic_slab()=default;

	void n_radiation_packages(size_t n) {
	  n_packages_ = n;
	}
	void single_scattering_albedo(double ssa) {
	  single_scattering_albedo_ = ssa;
	}
	void bottom_albedo(double ba) {
	  bottom_albedo_ = ba;
	}
	void source_zenith_angle(double sza) {
	  source_zenith_angle_ = sza;
	}
	double hemispherical_reflectance()
	// See Wikipedia reflectance for definition
	{
	  run(n_packages_);
	  auto nav_ =  geometry::navigator<content>(geometry_);
	  auto& v = nav_.find("layer");
	  auto& inward = v.content().inward_receiver();
	  auto& outward = v.content().outward_receiver();
	  return outward.radiant_flux()/inward.radiant_flux();
	}
	double directional_reflectance(double polar_angle,
				       double azimuth_angle)
	// See Wikipedia reflectance for definition
	{
	  //tbi
	  // should search for optimal angle. Gaussian aceptance?
	  // should gess direction directly towards detector
	  return 0;
	}
	double hemispherical_transmittance()
	// See Wikipedia transmittance for definition
	{
	  run(n_packages_);
	  auto nav_ =  geometry::navigator<content>(geometry_);
	  auto& vl = nav_.find("layer");
	  auto& incident = vl.content().inward_receiver();
	  auto& vb = nav_.find("bottom");
	  auto& transmitted = vb.content().inward_receiver();
	  return transmitted.radiant_flux()/incident.radiant_flux();
	}
	double directional_transmittance(double polar_angle,
					 double azimuth_angle)
	// See Wikipedia transmittance for definition
	{
	  //tbi
	  return 0;
	}
      private:
	void build_geometry() {
	  semi_infinite_box<content> layer;
	  semi_infinite_box<content> bottom;
	  geometry_.name("geometry");
	  layer.name("layer");
	  bottom.name("bottom");
	  layer.content().inward_receiver().activate();
	  layer.content().outward_receiver().activate();
	  bottom.content().inward_receiver().activate();
	  bottom.content().set_coating<coating::grey_lambert>(bottom_albedo_,1-bottom_albedo_);
	  geometry_.move_by({0,0,2});
	  layer.move_by({0,0,1});
	  layer.insert(bottom);
	  geometry_.clear();
	  geometry_.insert(layer);
	}
	void run(size_t n_packages) {
	  emitter emitter{{0,0,1.5},stokes{1,0,0,0},n_packages};
	  emitter.set_wavelength<monocromatic>(500e-9);
	  auto direction = unit_vector{constants::pi-source_zenith_angle_,0}; 
	  emitter.set_direction<unidirectional>(direction);
	  build_geometry();
	  mc_basic(geometry_).run(emitter, geometry_);
	}

      };

      double remote_sensing_reflectance(const basic_slab& bs) {
	//tbi
	return 0;
      }
      class henyey_greenstein_slab : public basic_slab {
	double asymmetr_factor_{0};
      public:	
      };
      class user_defined_slab : public basic_slab {
	
      public:
	void tabulated_phase_function(const std::string& filename) {
	}
      };
      
    }
  
  }
  */
}

#endif
