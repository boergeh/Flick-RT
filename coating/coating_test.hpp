#include "coating.hpp"

namespace flick {  
  using namespace constants;
  begin_test_case(coating_test) {
    const double pi = constants::pi;
    coating::white_lambert l;
    check_close(l.unpolarized_reflectance(),1,1e-12);
    coating::fresnel f;
    double n = 1.33;
    double k = 0;
    f.set(relative_refractive_index{n,k});
    f.set_incidence(rotation_about_y(pi));
    rotation r = f.reflection_rotation();
    rotation t = f.transmission_rotation();
    check_small(rms(r.z_direction(),{0,0,1}),1e-12);
    check_small(rms(t.z_direction(),{0,0,-1}),1e-12);

    double theta0 = 3*pi/4;
    r = rotation{rotation_about_x(-theta0)};
    f.set_incidence(r);
    r = f.reflection_rotation();
    t = f.transmission_rotation();
    double theta_t = pi-t.z_direction().theta(); 
    double theta_t_bench = asin(sin(pi-theta0)/n);
    check_close(theta_t,theta_t_bench,1e-12);

  } end_test_case()
}
