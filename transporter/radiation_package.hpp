#ifndef flick_radiation_package
#define flick_radiation_package

#include "../numeric/direction_generator.hpp"
#include "../polarization/stokes.hpp"

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
    size_t weighted_scattering_events_{0};
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
    void traveling_direction(const unit_vector& direction) {
      pose_.rotate_to(direction);
      weighted_scattering_events_ += stokes_.I(); 
    }
    void wavelength(double wl) {
      wavelength_ = wl;
    }
    void rotate_polarization_plane(double angle) {
      pose_.rotate_about_local_z(angle);
    }
    auto pose(){
      return pose_;
    }
    auto wavelength() {
      return wavelength_;
    }
    friend std::ostream& operator<<(std::ostream &os, const radiation_package& rp) {
      os << rp.pose_ << " " << rp.wavelength_ << " " << rp.stokes_ << " "
	 << rp.weighted_traveling_length_ << " "
	 << rp.weighted_scattering_events_;
      return os;
    }
  };
}

#endif
