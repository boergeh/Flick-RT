#include "spheres.hpp"
#include "water/pure_water.hpp"

namespace flick {
  begin_test_case(spheres_test) {
    using namespace flick;
    material::vacuum v;
    material::pure_water pw;
    double r = 0.01;
    log_normal_distribution sd(log(r),0.0001);
    double f = 0.1;
    material::spheres<log_normal_distribution,material::vacuum ,
    		material::pure_water,
    		parameterized_monodispersed_mie> s(f,sd,v,pw);
    s.set(wavelength{400e-9});
    check_close(s.scattering_coefficient(),3./2*f/r,0.1);

  } end_test_case()
}
