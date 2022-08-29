#ifndef flick_wall_interactor
#define flick_wall_interactor

#include "../geometry/volume.hpp"
#include "../component/content.hpp"
#include "../component/emitter.hpp"

namespace flick {
  using cube = geometry::cube<content>;
  using sphere = geometry::sphere<content>;
  using semi_infinite_box = geometry::semi_infinite_box<content>;

  class wall_interactor {
    geometry::navigator<content>& nav_;
    std::optional<pose> next_wall_intersection_;
    geometry::volume<content>* next_volume_;
    const coating::base* coating_;
    radiation_package& rp_;
    uniform_random& rnd_;
    coating::lambert_angle_generator& lag_;
    bool is_reflected_{false};
    bool is_transmitted_{true};
    bool is_moving_inward_;
    bool is_moving_outward_;
    unit_vector facing_surface_normal_;
  public:
    wall_interactor(geometry::navigator<content>& nav,
		    radiation_package& rp,
		    uniform_random& ur,
		    coating::lambert_angle_generator& lag)
      : nav_{nav}, rp_{rp}, rnd_{ur}, lag_{lag} {   
      next_wall_intersection_ = nav_.next_intersection(rp_.pose());
      next_volume_ = &nav_.next_volume(rp_.pose());
      is_moving_inward_ = is_moving_inward();
      is_moving_outward_ = is_moving_outward();
      
      geometry::volume<content>* v = next_volume_;
      if (is_moving_outward_) {
      	v = &nav_.current_volume();
      }
      if (v->content().has_coating() && next_wall_intersection_.has_value()) {
	coating_ = &v->content().coating();
	lag_ = coating::lambert_angle_generator{rnd_(0,1),rnd_(0,1)};
	double r = rnd_(0,1);
	is_reflected_ = (r < coating_->reflectivity());
	is_transmitted_ = (1-r < coating_->transmissivity());
      }
      facing_surface_normal_ = facing_surface_normal();
    }
    void move_to_wall() {
      if (next_wall_intersection_.has_value())
	rp_.move_to((*next_wall_intersection_).position());
    }
    void move_close_to_wall() {
      move_to_wall();
      //std::cout<<"pos " << nav_.current_volume().small_step()*facing_surface_normal_*1e7<<std::endl;
      rp_.move_by(nav_.current_volume().small_step()*facing_surface_normal_);
      //std::cout << rp_ << std::endl;
      //rp_.move(-nav_.current_volume().small_step()); 
    }
    void move_through_wall() {
      move_to_wall();
      rp_.move_by(-nav_.current_volume().small_step()*facing_surface_normal_);
      //rp_.move(nav_.current_volume().small_step());
    }    
    void interact_with_wall() {
      unit_vector& n = facing_surface_normal_;
      if (is_reflected_) {
	move_close_to_wall();
	reorient_traveling_direction(n);
	reshape_polarization(n);
	likelihood_scale_intensity();
      } else if (is_transmitted_) {
	move_through_wall();
	if (nav_.current_volume().content().has_coating() &&
	    next_wall_intersection_.has_value()) {
	  reorient_traveling_direction(-n);
	  reshape_polarization(-n);
	  likelihood_scale_intensity();
	}
	if (next_wall_intersection_.has_value()) {
	  increment_activated_receiver();
	}
	nav_.go_to(*next_volume_);
      } else {
	absorb_radiation_package();
      }
    }
    void increment_activated_receiver() {
      if (is_moving_inward_)
	next_volume_->content().inward_receiver().receive(rp_);
      else if (is_moving_outward_)
	nav_.current_volume().content().outward_receiver().receive(rp_);
    }
    auto& current_volume() const {
      return nav_.current_volume();
    }  
    void absorb_radiation_package() {
      rp_.scale_intensity(0);
    }
    void set_x_axis_parallel_with_plane_of_incidence(const unit_vector&
						     surface_normal) {
      pose p = rp_.pose();
      double ang = acos(dot(surface_normal, p.y_direction()));
      rp_.rotate_about_local_z(ang);
    }
    void reorient_traveling_direction(const unit_vector& facing_surface_normal) {      
      rp_.traveling_direction(facing_surface_normal);
      rp_.rotate_about_local_z(lag_.azimuth_angle());
      rp_.rotate_about_local_x(lag_.polar_angle());
    }
    void reshape_polarization(const unit_vector& facing_surface_normal) {
      set_x_axis_parallel_with_plane_of_incidence(facing_surface_normal);
      if (is_reflected_)
	rp_.reshape_polarization(coating_->reflection_mueller_matrix());
      else if (is_transmitted_)
	rp_.reshape_polarization(coating_->transmission_mueller_matrix());	
    }
    void likelihood_scale_intensity() {
      rp_.scale_intensity(1/lag_.brdf());
    }
    double distance_to_wall() {
      if (!next_wall_intersection_.has_value())
	return 0;
      return norm((*next_wall_intersection_).position()-rp_.pose().position()); 
    }
  private:
    unit_vector facing_surface_normal() const {
      unit_vector n = (*next_wall_intersection_).z_direction();
      if (dot(n,rp_.pose().z_direction()) < 0)
	return n;
      return -n;
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
  };
}

#endif
