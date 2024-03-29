#ifndef flick_radiation_package
#define flick_radiation_package

#include "../numeric/direction_generator.hpp"
#include "../polarization/stokes.hpp"
#include "../polarization/algorithm.hpp"

namespace flick {
  class radiation_package
  {
    flick::pose pose_;
    double wavelength_{500e-9};
    flick::stokes stokes_{1,0,0,0};
    double traveling_length_{0};
    unit_vector emission_direction_;
    //size_t scattering_events_{0};?
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
      traveling_length_ += distance;
    }
    void move_to(const vector& position) {
      double distance = norm(position-pose_.position());
      pose_.move_to(position);
      traveling_length_ += distance;
    }
    void move_by(const vector& v) {
      move_to(pose_.position()+v);
    }
    void wavelength(double wl) {
      wavelength_ = wl;
    }
    void rotate_to(const rotation& r) {
      pose_ = {pose_.position(),r()};
    }
    void interact_with_matter(const mueller& m) { 
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
      stokes_.rotate(-angle);
    }
    bool is_empty() const {
      return (stokes_.I() < 1e-9);
    }
    double traveling_length() const {
      return traveling_length_;
    }
    const unit_vector& emission_direction() const {
      return emission_direction_;
    }
    void emission_direction(const unit_vector& ed) {
      emission_direction_ = ed;
    }
    const flick::pose& pose() const {
      return pose_;
    }
    auto wavelength() const {
      return wavelength_;
    }
    const flick::stokes& stokes() const {
      return stokes_;
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const radiation_package& rp) {
      os << "xyz " << rp.pose_.position() << ", dir "<<rp.pose_.z_direction()
	 << ", wl " << rp.wavelength_
	 << ", " << rp.stokes_ << ", length "
	 << rp.traveling_length_ << ", emission_dir "<<rp.emission_direction_
	 << " ";
      return os;
    }
  };
}

#endif
