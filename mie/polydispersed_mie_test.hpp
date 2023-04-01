#include "polydispersed_mie.hpp"
#include "parameterized_monodispersed_mie.hpp"
#include "monodispersed_mie.hpp"

namespace flick {
  begin_test_case(poly_mie_test_A) {
    stdcomplex m_host{1,0};
    stdcomplex m_sphere{1.3,1e-4};
    
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     500e-9);
    double r = 10e-6;
    mono_mie.radius(r);
    mono_mie.angles({0,1.5,3.14});    
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
    poly_mie.precision(4);
    check_close(poly_mie.absorption_cross_section(),
    		mono_mie.absorption_cross_section(),0.02);
    check_close(poly_mie.scattering_cross_section(),
    		mono_mie.scattering_cross_section(),0.02);
    check_close(poly_mie.scattering_matrix_element(0,0)[0],
    		mono_mie.scattering_matrix_element(0,0)[0],0.02);

    check_close(distribution.particles_per_volume(1),
		1/(4./3*constants::pi*pow(r,3)),0.1);
  } end_test_case()

   begin_test_case(poly_mie_test_B) {
    double pi = constants::pi;
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.33 + 1e-5i;
    double wl = 500e-9;
    double r = 10e-6;
    stdvector angles = range(0,pi,100).linspace();
    monodispersed_mie mono_mie(m_host,m_sphere,wl);
    mono_mie.angles(angles);
    log_normal_distribution sd{log(r),0.01};
    polydispersed_mie poly_mie(mono_mie,sd);
    //std::cout << poly_mie.scattering_cross_section();
    //std::cout << poly_mie.scattering_matrix_element(0,0);
  } end_test_case()

   begin_test_case(poly_mie_test_C) {
    stdcomplex m_host = {1,0};
    stdcomplex m_sphere = {1.3,0.0001};
    double r = 10e-6;
    monodispersed_mie mono_mie(m_host,m_sphere,500e-9);
    mono_mie.angles({0,3,3.14});

    log_normal_distribution sd{log(r),0.0};
    polydispersed_mie poly_mie(mono_mie,sd);
    poly_mie.precision(4);
    //double Csca = poly_mie.scattering_cross_section();
    //std::cout << Csca << std::endl;
    //check_close(Csca,5.673e-9,10);
    //std::cout << poly_mie.absorption_cross_section();
    //std::cout << poly_mie.scattering_matrix_element(0,0) << std::endl;

  } end_test_case()
}
