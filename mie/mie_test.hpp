#include "mie.hpp"

namespace flick {
  begin_test_case(mie_test) {
    refractive_index m_host{1,0};
    refractive_index m_sphere{1.3,1e-4};
    
    parameterized_monodisperesed_mie mono_mie(m_host,
					    m_sphere,
					    wavelength(500e-9));
    double r = 10e-6;
    mono_mie.radius(r);
    mono_mie.precision(3);
    
    //std::cout << mono_mie.precision() << std::endl;
    //check_small(mono_mie.absorption_cross_section(), 1e-12);
    //check_close(mono_mie.extinction_cross_section(), 2*3.14159*pow(r,2),0.1);

    check(mono_mie.scattering_matrix_element(0,0)[0] >
	  mono_mie.scattering_matrix_element(0,0)[1]);

    //polydispersed_mie poly_mie(volume_fraction,
    //			       log_normal_distribution{log(r),1e-5},
    //			       mono_mie);
    polydispersed_mie poly_mie(mono_mie,
			       log_normal_distribution{log(r),1e-5});
    check_close(poly_mie.absorption_cross_section(),
    		mono_mie.absorption_cross_section(),0.02);
    check_close(poly_mie.scattering_cross_section(),
    		mono_mie.scattering_cross_section(),0.02);
    check_close(poly_mie.scattering_matrix_element(0,0)[2],
    		mono_mie.scattering_matrix_element(0,0)[2],0.02);
    
  } end_test_case()
}
