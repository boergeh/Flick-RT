#include "mie.hpp"
#include "monodispersed_mie.hpp"

namespace flick {
  begin_test_case(mie_test_A) {
    stdcomplex m_host{1.3,0};
    stdcomplex m_sphere{1.0,0};
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     500e-9);
    const double pi = constants::pi;
    double r = 0.001;
    mono_mie.radius(r);
    check_close(mono_mie.scattering_cross_section(),2*pi*pow(r,2),0.1);
    check_small(mono_mie.absorption_cross_section(),1e-12);
  } end_test_case()

  begin_test_case(mie_test_B) {
    stdcomplex m_host{1,0};
    stdcomplex m_sphere{1.3,1e-4};
    
    parameterized_monodispersed_mie mono_mie(m_host,
					     m_sphere,
					     500e-9);
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

  begin_test_case(mie_test_C) {
    //bessel_first_kind j(120.3+200.*1i,9);
    //std::cout << j.times_z_derivatives();
    stdcomplex m1 = 1.0+0.05i;
    stdcomplex m2 = 1.53;
    double wl = 2*constants::pi;
    double r = 10/(2*constants::pi)*wl;
    monodispersed_mie mie(m1,m2,wl);
    mie.radius(r);
    auto [a, b] = mie.ab_coefficients();
    double percent = 1e-10; 
    check_close(a[1].real(),0.82786371508743,percent);
    check_close(b[1].imag(),0.91474090929954,percent);
    check_close(a[10].imag(),0.28299838641514,percent);
    check_close(b[10].real(),-0.04440119340674,percent);
    m1 = 1.33 + 0.1i;
    m2 = 1.0;
    r = 2500;
    monodispersed_mie mie_b(m1,m2,wl);
    mie_b.radius(r);
    auto [a_b, b_b] = mie_b.ab_coefficients();
    check_close(a_b[1].real(),4.3914709187499176e216,percent);
  } end_test_case()
}
