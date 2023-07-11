#include "algorithm.hpp"
#include "isotropic_mueller.hpp"
#include "rayleigh_mueller.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"

namespace flick {
  begin_test_case(algorithm_test) {
    using namespace constants;
    mueller m;
    m.add(0,0,1);
    m.add(0,1,1);
    m.add(1,0,1);
    stokes s1 = stokes{1,1,0,0};
    stokes s2 = m * s1;
    
    mueller im = isotropic_mueller();
    stokes s3(1,0,0,0);
    stokes s4 = im*s3;
    check_close(s4.I()*4*pi,1.0_pct);
    
    pl_function f;
    size_t n = 150;
    auto ang = range(0,pi,n).linspace();
    for (size_t i=0; i<n; ++i) {
      mueller rm = rayleigh_mueller(ang[i],0.2);
      f.append({ang[i], rm(0).value * sin(ang[i])});
    }
    check_close(2*pi*f.integral(), 1, 0.01_pct);
    
    mueller rm = rayleigh_mueller(pi/2,0);
    stokes s5{1,0,0,0};
    stokes s6 = rm * s5;
    check_small(s6.V());
    check_close(s6.degree_of_polarization(),1.0_pct);
    check_close(s6.Q()/s6.I(),-1);
    check_small(s6.eccentricity_angle());
    check_close(s6.rotation_angle(),pi/2);
    
    stokes s7 = s6.rotate(pi/2);
    check_close(s7.Q()/s7.I(),1.0_pct);
    check_small(s7.U()/s7.I());
    check_small(s7.rotation_angle());
    s7.rotate(-pi/4);   
    check_small(s7.Q()/s7.I());
    check_close(s7.U()/s7.I(),1.0_pct);
    check_small(s7.V()/s7.I());
    check_close(s7.rotation_angle(),pi/4);
    check_small(s7.eccentricity_angle());
  } end_test_case()
}

