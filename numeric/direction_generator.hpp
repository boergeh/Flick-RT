#ifndef flick_direction_generator
#define flick_direction_generator

#include "uniform_random.hpp"
#include "pose.hpp"

namespace flick {
  class direction_generator {
    uniform_random rnd;
  public:
    unit_vector isotropic() const {
      double phi = rnd(0,2*pi);
      double theta = acos(rnd(-1,1));
      return unit_vector{theta,phi};
    }
    unit_vector conic(double solid_angle, const unit_vector& cone_direction) const {
      double phi = rnd(0,2*pi);
      double theta = acos(rnd(1-solid_angle/(2*pi),1));
      pose p;
      p.rotate_to(cone_direction);
      p.rotate_about_local_z(phi);
      p.rotate_about_local_y(theta);
      return p.z_direction();
    }
    unit_vector lambertian(const unit_vector& surface_normal) {
      double phi = rnd(0,2*pi);
      double theta = asin(sqrt(rnd(0,1)));
      pose p;
      p.rotate_to(surface_normal);
      p.rotate_about_local_z(phi);
      p.rotate_about_local_y(theta);
      return p.z_direction();
    }
  };
}
#endif
