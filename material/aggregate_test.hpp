#include "aggregate.hpp"
#include "aerosols/aerosols.hpp"

namespace flick {
  using namespace constants;
  using namespace flick;
  
  begin_test_case(aggregate_test_hg) {
    material::white_isotropic m1{1.0};
    stdvector ang{1,2,3.14};
    material::aggregate ma{ang};
    ma.add(m1);
    check_close(ma.scattering_coefficient(),1);
    material::white_isotropic m2{2.0};
    ma.add(m2);
    ma.add(m2);
    check_close(ma.scattering_coefficient(),5);
    check_small(ma.absorption_coefficient());
    check_small(ma.absorption_optical_depth(100));

    stdvector h{0,1,2,10};
    material::aggregate ma2(ang,h);
    ma2.add(m1);
    check_close(ma2.scattering_coefficient(),1);
    ma2.add(m1,1,2);
    check_close(ma2.scattering_optical_depth(10),11);
    check_small(ma2.absorption_optical_depth(100));
  } end_test_case()
  
  begin_test_case(aggregate_test_aerosols) {
    material::rural_aerosols ra;
    material::aggregate ag({0,1,3.14}, {0,1e3,100e3});
    ag.set_direction({0,0,1});
    ag.add(ra);
    check_close(ag.scattering_optical_depth(100e3),
		ra.scattering_optical_depth(100e3));
    unit_vector u{0,0};
    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0));
    ag.add(ra);
    check_close(ag.absorption_optical_depth(100e3),
    		2*ra.absorption_optical_depth(100e3));
    ag.add(ag);
    check_close(ag.absorption_optical_depth(100e3),
    		4*ra.absorption_optical_depth(100e3));

    u = unit_vector{3.14,0};
    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0));
    material::urban_aerosols ua;
    ag.add(ua);
    check_close(ag.absorption_optical_depth(100e3),
    		4*ra.absorption_optical_depth(100e3)+
		ua.absorption_optical_depth(100e3));

    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0)*0.999,0.07_pct); 
  } end_test_case()

  begin_test_case(aggregate_test_aerosols_hg) {
    double cl_scat_coef = 1e-4;
    double zc_low = 410;
    double zc_high = 890;
    double toa = 120e3;
    material::urban_aerosols aero;
    material::white_isotropic cloud{cl_scat_coef};
    stdvector h = {0,100,700,800,2000,toa};
    stdvector ang = {0,3.14};
    material::aggregate atmosphere{ang,h};
    atmosphere.add(aero);
    atmosphere.add(cloud,zc_low,zc_high);
    double tau1 = atmosphere.optical_depth(toa);
    double tau2 = aero.optical_depth(toa);
    double tau3 = cl_scat_coef*(zc_high-zc_low);
    check_close(tau1,tau2+tau3);
    atmosphere.add(cloud,pe_function{{zc_low,zc_high},{0,0}});
    check_close(tau1,tau2+tau3);    
  } end_test_case()
}
