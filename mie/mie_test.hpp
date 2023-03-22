#include "mie.hpp"

namespace flick {
  begin_test_case(mie_test_A) {
    refractive_index m_host{1.3,0};
    refractive_index m_sphere{1.0,0};
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     wavelength(500e-9));
    const double pi = constants::pi;
    double r = 0.001;
    mono_mie.radius(r);
    check_close(mono_mie.scattering_cross_section(),2*pi*pow(r,2),0.1);
    check_small(mono_mie.absorption_cross_section(),1e-12);
  } end_test_case()

  begin_test_case(mie_test_B) {
    refractive_index m_host{1,0};
    refractive_index m_sphere{1.3,1e-4};
    
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     wavelength(500e-9));
    double r = 10e-6;
    mono_mie.radius(r);
    mono_mie.precision(4);

    double mu_v = log(r);
    double sigma_v = 1e-9;
    auto [mu, sigma] =
      log_normal_distribution::from_volume_distribution(mu_v,sigma_v);
    log_normal_distribution distribution{mu,sigma};

    check_close(mono_mie.scattering_cross_section()
		+mono_mie.absorption_cross_section(),
		2*3.14159*pow(r,2),0.1);

    check(mono_mie.scattering_matrix_element(0,0)[0] >
	  mono_mie.scattering_matrix_element(0,0)[1]);

    polydispersed_mie poly_mie(mono_mie,distribution);
    check_close(poly_mie.absorption_cross_section(),
    		mono_mie.absorption_cross_section(),0.02);
    check_close(poly_mie.scattering_cross_section(),
    		mono_mie.scattering_cross_section(),0.02);
    check_close(poly_mie.scattering_matrix_element(0,0)[0],
    		mono_mie.scattering_matrix_element(0,0)[0],0.02);

    check_close(distribution.particles_per_volume(1),
		1/(4./3*constants::pi*pow(r,3)),0.1);
    
  } end_test_case()
}
