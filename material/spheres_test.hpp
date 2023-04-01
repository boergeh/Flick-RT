#include "spheres.hpp"
#include "water/pure_water.hpp"
namespace flick {
  begin_test_case(spheres_test_A) {
    using namespace flick;
    material::vacuum v;
    material::pure_water pw;
    double r = 10e-6;
    log_normal_distribution sd(log(r),0.0001);
    double f = 0.1;  
    material::spheres<log_normal_distribution,material::vacuum ,
    		material::pure_water,
    		parameterized_monodispersed_mie> s(f,sd,v,pw);
    s.set(wavelength{400e-9});
    check_close(s.scattering_coefficient(),3./2*f/r,0.001);
    
    material::bubbles_in_ice<monodispersed_mie> bi(1,log(1e-10),0.0001);
    auto m1 = bi.mueller_matrix(unit_vector{0,0});
    auto m2 = rayleigh_mueller(0,0);
    check_close(m1.value(0,0),m2.value(0,0),1e-4);
    check_close(m1.value(2,2),m2.value(2,2),1e-4);
    check_close(m1.value(3,3),m2.value(3,3),1e-4);

    double theta = constants::pi/2;
    m1 = bi.mueller_matrix(unit_vector{theta,0});
    m2 = rayleigh_mueller(theta,0);
    check_close(m1.value(0,0),m2.value(0,0),1e-4);
    check_close(m1.value(0,1),m2.value(0,1),1e-4);
    check_close(m1.value(1,0),m2.value(1,0),1e-4);
  } end_test_case()

  begin_test_case(spheres_test_B) {
    using namespace flick;
   
    material::bubbles_in_ice<parameterized_monodispersed_mie> bu(1,log(1e-6),0.0001);
    auto m = bu.mueller_matrix(unit_vector{1,0});
    material::brines_in_ice<parameterized_monodispersed_mie> br(1,log(1e-6),0.0001,100);
    m = br.mueller_matrix(unit_vector{1,0});
    material::water_cloud<parameterized_monodispersed_mie> cl(1,log(1e-6),0.0001);
    m = cl.mueller_matrix(unit_vector{1,0});
  } end_test_case()
}
