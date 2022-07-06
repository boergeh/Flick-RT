#ifndef flick_polarization_algorithm_test
#define flick_polarization_algorithm_test

#include "algorithm.hpp"
#include "isotropic_mueller.hpp"
#include "rayleigh_mueller.hpp"
//#include "isotropic_surface_mueller.hpp"
//#include "dielectric_surface_mueller.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"

namespace flick {
  begin_test_case(algorithm_test) {
    double pi = constants::pi;

    isotropic_mueller im;
    stokes s1(1,0,0,0);
    stokes s2 = im*s1;
    check_close(s2.I()*4*pi,1,1e-12);

    pl_function f;
    size_t n = 150;
    auto ang = range(0,pi,n).linspace();
    for (size_t i=0; i<n; ++i) {
      rayleigh_mueller rm(ang[i],0.2);
      f.append({ang[i],rm.element(0,0)*sin(ang[i])});
    }
    check_close(2*pi*f.integral(), 1, 0.01);

    rayleigh_mueller rm(pi/2);
    stokes s3 = rm*s1;
    check_small(s3.V(),1e-9);
    check_close(s3.degree_of_polarization(),1,1e-12);
    check_close(s3.Q()/s3.I(),1,1e-12);
    check_small(s3.rotation_angle(),1e-12);
    
    stokes s4 = s3.rotate(pi/2);
    check_close(s4.Q()/s4.I(),-1,1e-12);
    check_close(s4.rotation_angle(),pi/2,1e-12);
    s4.rotate(-pi/4);
    check_close(s4.U()/s4.I(),1,1e-12);
      
  } end_test_case()
}

#endif
