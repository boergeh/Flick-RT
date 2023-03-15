#include "mie.hpp"

namespace flick {
  begin_test_case(mie_test) {
    refractive_index m_host{1,0};
    refractive_index m_sphere{1.3,1e-6};
    
    parameterized_monodisperesed_mie mono_mie(m_host,
					    m_sphere,
					    wavelength(500e-9));
    double r = 10e-6;
    mono_mie.radius(r);
    mono_mie.precision(4);
    //std::cout << mono_mie.precision() << std::endl;

    check_small(mono_mie.absorption_cross_section(), 1e-12);
    //check_close(mono_mie.extinction_cross_section(), 2*3.14159*pow(r,2),0.1);
    check(mono_mie.scattering_function(0,0)[0] >
	  mono_mie.scattering_function(0,0)[1]);
 
    polydispersed_mie poly_mie(mono_mie,log_normal_distribution{log(r),1e-5});
    check_close(mono_mie.extinction_cross_section(),
    		poly_mie.extinction_cross_section(),0.01);
    check_close(mono_mie.absorption_cross_section(),
    		poly_mie.absorption_cross_section(),0.01);
    
  } end_test_case()
}
