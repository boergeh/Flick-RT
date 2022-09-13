#ifndef flick_material_interactor
#define flick_material_interactor

namespace flick {
namespace transporter {  
  class material_interactor {
    radiation_package& rp_;
    material::base& m_;
    uniform_random& rnd_;
    double scattering_optical_depth_;
    double g_;
    unit_vector scattering_direction_;
    double scattering_angle_;
    double distance_to_scattering_;
    double scattering_polar_angle_;
    double scattering_azimuth_angle_;
  public:
    material_interactor(radiation_package& rp,
			material::base& m,
			uniform_random& rnd,
			double scattering_optical_depth,
			double sampling_asymmetry_factor)
      : rp_{rp}, m_{m}, rnd_{rnd}, g_{sampling_asymmetry_factor},
	scattering_optical_depth_{scattering_optical_depth} {
      m_.set(rp_.pose());
      distance_to_scattering_ = m_.scattering_distance(scattering_optical_depth_);
    }
    void move_to_scattering_event() {
      rp_.move(distance_to_scattering_);
    }
    void deposite_energy_to_heat(double distance) {
      rp_.scale_intensity(exp(-m_.absorption_optical_depth(distance)));
    }
    double distance_to_scattering() {
      return distance_to_scattering_;
    }
    void find_scattering_direction() {
      scattering_polar_angle_ =
	henyey_greenstein{g_}.inverted_accumulated_angle(rnd_(0,1));
      scattering_azimuth_angle_ = rnd_(0,2*constants::pi);
      pose p = rp_.pose();
      p.rotate_about_local_z(scattering_azimuth_angle_);
      p.rotate_about_local_y(scattering_polar_angle_);
      scattering_direction_ =  p.z_direction();
    }
    void reorient_traveling_direction() {
      // unit_vector olddir = rp_.pose().z_direction();
      rp_.rotate_about_local_y(scattering_polar_angle_);
      //double mu = dot(rp_.pose().z_direction(),olddir);
      //std::cout << mu << std::endl;
    }
    void reshape_polarization() {
      rp_.interact_with_matter(m_.mueller_matrix(scattering_direction_));
    }
    void likelihood_scale_intensity() {
      double hg = henyey_greenstein{g_}.phase_function(scattering_polar_angle_);
      rp_.scale_intensity(1/hg);
    }
    void scatter() {
      move_to_scattering_event();
      find_scattering_direction();
      likelihood_scale_intensity();
      make_x_axis_parallel_with_scattering_plane();
      reshape_polarization();    
      reorient_traveling_direction();
    }
  private:
    void make_x_axis_parallel_with_scattering_plane() {
      //const pose& p = rp_.pose();
      //vector n = cross(p.z_direction(),scattering_direction_);
      //double ang = acos(dot(n, p.x_direction()));
      //rp_.rotate_about_local_z(ang);
      rp_.rotate_about_local_z(scattering_azimuth_angle_);
    }
  };
}
}

#endif
