#include "spheres.hpp"
#include "water/pure_water.hpp"

namespace flick {
  begin_test_case(spheres_test) {
    using namespace flick;
    material::vacuum v;
    material::pure_water pw;
    log_normal_distribution sd(10e-6,0.0001);
    //material::spheres<log_normal_distribution,material::vacuum ,
    //			material::pure_water,
    //			parameterized_monodispersed_mie> s(0.15,sd,v,pw);
    material::parameterized_spheres s(0.15,sd,v,pw);
    s.set(wavelength{1500e-9});

    material::water_cloud<log_normal_distribution,
			  parameterized_monodispersed_mie> wc(1e-6,sd);

  } end_test_case()
}
