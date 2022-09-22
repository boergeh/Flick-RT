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
    geometry::volume<content>* current_volume_;
    geometry::volume<content>* next_volume_;
    coating::base* coating_;
    radiation_package& rp_;
    uniform_random& rnd_;
    bool is_reflected_{false};
    bool is_transmitted_{true};
    bool is_moving_inward_;
    unit_vector facing_surface_normal_;
  public:
    wall_interactor(geometry::navigator<content>& nav,
		    radiation_package& rp,
		    uniform_random& ur)
      : nav_{nav}, rp_{rp}, rnd_{ur} {   
      next_wall_intersection_ = nav_.next_intersection(rp_.pose());
      current_volume_ = &nav_.current_volume();
      next_volume_ = &nav_.next_volume(rp_.pose());
      is_moving_inward_ = is_moving_inward();
      if (!current_volume_->content().has_coating()) {	  
	current_volume_->content().coat<coating::fresnel>();
      }
      if (!next_volume_->content().has_coating()) {	  
	next_volume_->content().coat<coating::fresnel>();
      }
      if (!next_volume_->content().has_material()) {	  
	next_volume_->content().fill<material::vacuum>();
      } 
     
      facing_surface_normal_ = facing_surface_normal();
      move_to_wall();
      align_rp_x_axis_with_wall();
      set_coating();
      double r = rnd_(0,1);
      if (coating_!=nullptr) {
	is_reflected_ = (r < coating_->unpolarized_reflectance());
	is_transmitted_ = (1-r < coating_->unpolarized_transmittance());
      }
    }
    void interact_with_wall() {
      if (!next_wall_intersection_.has_value()) {
	nav_.go_to(*next_volume_);
      }
      else if (is_reflected_) {
	rp_.interact_with_matter(coating_->reflection_mueller_matrix());
	rp_.scale_intensity(1/coating_->unpolarized_reflectance());
	rp_.rotate_to(coating_->reflection_rotation());
	receive_reflected_packages();
	step_back_from_wall();
      }
      else if (is_transmitted_) {
	receive_transmitted_packages();
	if (coating_!=nullptr) {
	  rp_.interact_with_matter(coating_->transmission_mueller_matrix());
	  rp_.scale_intensity(1/coating_->unpolarized_transmittance());
	  rp_.rotate_to(coating_->transmission_rotation());
	}
	step_through_wall();
	nav_.go_to(*next_volume_);
      }
      else {
	receive_transmitted_packages();
	absorb_radiation_package();
      }
    }
  private:
    void receive_reflected_packages() {
      if (is_moving_inward_)
	next_volume_->content().outward_receiver().receive(rp_);
      else
	current_volume_->content().inward_receiver().receive(rp_);
    }
    void receive_transmitted_packages() {
      if (is_moving_inward_)
	next_volume_->content().inward_receiver().receive(rp_);
      else
	current_volume_->content().outward_receiver().receive(rp_);
    }
    void absorb_radiation_package() {
      rp_.scale_intensity(0);
    }
    void set_coating() {
      if (is_moving_inward_)
	coating_ = &next_volume_->content().coating();
      else
	coating_ = &current_volume_->content().coating();
      if (coating_!=nullptr) {
	coating_->set_incidence(rotation{rp_.pose().rotation()});
	coating_->set(facing_surface_normal_);
	unit_interval ui_polar{rnd_(0,1)};
	unit_interval ui_azimuth{rnd_(0,1)};
	coating_->set_direction_parameters(ui_polar,ui_azimuth);
	coating_->set(relative_refractive_index());
      }
    }
    std::complex<double> relative_refractive_index() {
      std::complex<double> m1 = current_volume_->
	content().material().refractive_index();
      std::complex<double> m2 = next_volume_->
	content().material().refractive_index();
      return m2 / m1;
    }
    void move_to_wall() {
      if (next_wall_intersection_.has_value())
	rp_.move_to((*next_wall_intersection_).position());
    }
    void step_back_from_wall() {
      rp_.move_by(current_volume_->small_step()*facing_surface_normal_);
    }
    void step_through_wall() {
      rp_.move_by(-current_volume_->small_step()*facing_surface_normal_);
    }    
    void align_rp_x_axis_with_wall()
    // Not too proud over this..
    {
      pose p = rp_.pose();
      vector normal = cross(p.z_direction(),facing_surface_normal_);
      if (norm(normal) > std::numeric_limits<double>::epsilon()*10) {
	double d = dot(normalize(normal),p.x_direction());
	if (d <= 1 && d >= -1) // avoid rounding errors
	  rp_.rotate_about_local_z(acos(d));
	if (dot(normal,rp_.pose().x_direction()) < 0)
	  rp_.rotate_about_local_z(constants::pi);
      }
    }
    unit_vector facing_surface_normal() const {
      unit_vector n = (*next_wall_intersection_).z_direction();
      if (dot(n,rp_.pose().z_direction()) < 0)
	return n;
      return -n;
    }
    bool is_moving_inward() const {
      //if (!next_wall_intersection_.has_value())
      //	return false;
      return nav_.is_moving_inward(*next_wall_intersection_,rp_.pose());
    }    
  };
}

#endif
