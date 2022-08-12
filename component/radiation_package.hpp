#ifndef flick_radiation_package
#define flick_radiation_package

#include "../numeric/direction_generator.hpp"
#include "../polarization/stokes.hpp"
#include "../polarization/mueller.hpp"
#include "../polarization/algorithm.hpp"

namespace flick {
  class radiation_package
  // First Stokes parameter is intensity weight. All radiation
  // variables may change during transport. All stokes
  // parameters have units of 'per solid angle'.
  {
    pose pose_;
    double wavelength_{500e-9};
    stokes stokes_{1,0,0,0};
    double weighted_traveling_length_{0};
    //size_t weighted_scattering_events_{0};
    //bool do_not_scatter_?
  public:
    radiation_package() = default;
    radiation_package(const pose& p, const stokes& s)
      : pose_{p}, stokes_{s} {}
    void scale_intensity(double factor) {
      stokes_.scale(factor);
    }
    void move(double distance) {
      pose_.move_by(pose_.z_direction()*distance);
      weighted_traveling_length_ += stokes_.I()*distance;
    }
    void move_to(const vector& position) {
      double distance = norm(position-pose_.position());
      move(distance);
    }
    void wavelength(double wl) {
      wavelength_ = wl;
    }
    void traveling_direction(const unit_vector& direction) {
      pose_.rotate_to(direction);
    }
    void x_direction_parallel_with_plane_of_incidence(const unit_vector&
						  surface_normal) {
      double ang = acos(dot(surface_normal, pose_.y_direction()));
      pose_.rotate_about_local_z(ang);
      stokes_.rotate(-ang);
    }
    void x_direction_parallel_with_scattering_plane(const unit_vector&
						  scattering_direction) {
      vector n = cross(pose_.z_direction(),scattering_direction);
      double ang = acos(dot(n, pose_.x_direction()));
      pose_.rotate_about_local_z(ang);
      stokes_.rotate(-ang);
    }
    void change_polarization(const mueller& m)
    {
      stokes_ = m * stokes_;
    }
    void rotate_about_local_x(double angle) {
      pose_.rotate_about_local_x(angle);
    }
   
    void rotate_about_local_y(double angle) {
      pose_.rotate_about_local_y(angle);
    }
    void rotate_about_local_z(double angle) {
      pose_.rotate_about_local_z(angle);
    }
 
    bool is_empty() const {
      double epsilon = 1e-9;
      if (stokes_.I() < 1e-9)
	return true;
      return false;
    }
    auto weighted_traveling_length() const {
      return weighted_traveling_length_;
    }
    auto pose() const {
      return pose_;
    }
    auto wavelength() {
      return wavelength_;
    }
    friend std::ostream& operator<<(std::ostream &os, const radiation_package& rp) {
      os << "xyz: " << rp.pose_.position() << ", dir: "<<rp.pose_.z_direction()
	 << ", wl: " << rp.wavelength_
	 << ", " << rp.stokes_ << ", length: "
	 << rp.weighted_traveling_length_ << " ";
	//<< rp.weighted_scattering_events_;
      return os;
    }
  };
}

#endif
