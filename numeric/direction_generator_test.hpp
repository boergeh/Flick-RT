#include "direction_generator.hpp"

namespace flick {
  begin_test_case(direction_generator_test) {
    direction_generator dg;

    double theta = dg.isotropic().theta(); 
    double phi = dg.isotropic().phi(); 
    check(theta >= 0 && theta <= pi);
    check(phi >= 0 && theta <= 2*pi);

    unit_vector cone_direction = {pi,pi/2};
    theta = dg.conic(2*pi,cone_direction).theta(); 
    phi = dg.conic(2*pi,cone_direction).phi(); 
    check(theta >= 0 && theta <= pi);
    check(phi >= 0 && theta <= pi);
        
    unit_vector surface_normal = {pi,pi/3};
    theta = dg.lambertian(surface_normal).theta(); 
    phi = dg.lambertian(surface_normal).phi(); 
    check(theta >= pi/2 && theta <= pi);
    check(phi >= 0 && theta <= 2*pi);

  } end_test_case()
}
