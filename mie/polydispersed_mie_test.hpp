#include "polydispersed_mie.hpp"
#include "parameterized_monodispersed_mie.hpp"
#include "monodispersed_mie.hpp"
#include "../environment/input_output.hpp"

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
    double sigma_v = 1e-6;
    auto [mu, sigma] =
      log_normal_distribution::from_volume_distribution(mu_v,sigma_v);
    log_normal_distribution distribution{mu,sigma};

    check_close(mono_mie.scattering_cross_section()
		+mono_mie.absorption_cross_section(),
		2*3.14159*pow(r,2),0.1_pct);

    check(mono_mie.scattering_matrix_element(0,0)[0] >
	  mono_mie.scattering_matrix_element(0,0)[1]);

    polydispersed_mie poly_mie(mono_mie,distribution);
    double p = 0.02;
    poly_mie.percentage_accuracy(p);
    check_close(poly_mie.absorption_cross_section(),
    		mono_mie.absorption_cross_section(),p);
    check_close(poly_mie.scattering_cross_section(),
    		mono_mie.scattering_cross_section(),p);
    check_close(poly_mie.scattering_matrix_element(0,0)[0],
    		mono_mie.scattering_matrix_element(0,0)[0],p);

    check_close(distribution.particles_per_volume(1),
		1/(4./3*constants::pi*pow(r,3)),0.1_pct);
  } end_test_case()

  begin_test_case(poly_mie_test_B) {
    stdcomplex m_host = {1,0};
    stdcomplex m_sphere = {1.3,0.0};
    double r = 100e-6;
    parameterized_monodispersed_mie mono_mie(m_host,m_sphere,500e-9);
    log_normal_distribution sd{log(r),0.1};
    polydispersed_mie poly_mie(mono_mie,sd);
    double p = 0.01;
    poly_mie.percentage_accuracy(p);
    check_close(poly_mie.scattering_efficiency(),2,p);
  } end_test_case()

  begin_test_case(poly_mie_test_C) {
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.03+1e-9i;
    double wl = 355e-9;
    monodispersed_mie mono_mie(m_host,m_sphere,wl);  
    mono_mie.angles({0,constants::pi});
    log_normal_distribution sd{log(1e-6),0.2};
    polydispersed_mie poly_mie(mono_mie,sd);
     poly_mie.percentage_accuracy(0.1);
    stdvector F11 = poly_mie.scattering_matrix_element(0,0);
    check(F11[0]/F11[1] > 10);
  } end_test_case()

   begin_test_case(poly_mie_test_D) {
    stdcomplex m_host = 1.33;
    stdcomplex m_sphere = 1.5+1e-10i;
    double wl = 500e-9;
    monodispersed_mie mono_mie(m_host,m_sphere,wl);  
    mono_mie.angles({1.88});
    log_normal_distribution sd{log(8e-6),0.1};
    polydispersed_mie poly_mie(mono_mie,sd);
    double p = 5;
    poly_mie.percentage_accuracy(p);
    double F11 = poly_mie.scattering_matrix_element(0,0)[0];
    double bench = 6.131e-13;
    check_close(F11,bench,p);
    check_fast(2000 * cpu_duration());
  } end_test_case()

   begin_test_case(poly_mie_test_t_matrix) {
    stdcomplex m_host = 1.36;
    stdcomplex m_sphere = 1.580+0.3915e-3i;
    double wl = 368e-9;
    monodispersed_mie mono_mie(m_host,m_sphere,wl);  
    mono_mie.angles({0});
    log_normal_distribution sd{log(182e-9),log(1.11)};
    polydispersed_mie poly_mie(mono_mie,sd);
    poly_mie.percentage_accuracy(0.1);
    double t_matrix_C_scat = 0.9849e5*pow(1e-9,2); 
    double t_matrix_C_abs = 0.9894e5*pow(1e-9,2)-t_matrix_C_scat;
    double t_matrix_F11 = 20.77 * t_matrix_C_scat/(4*constants::pi);
    check_close(poly_mie.scattering_cross_section(), t_matrix_C_scat,0.4_pct);
    check_close(poly_mie.absorption_cross_section(), t_matrix_C_abs,0.8_pct);
    check_close(poly_mie.scattering_matrix_element(0,0)[0],t_matrix_F11,0.8_pct);
  } end_test_case()

   begin_test_case(poly_mie_test_no_absorption) {
    stdcomplex m_host = {1,0};
    stdcomplex m_sphere = {1.3,0};
    double r = 1e-6;
    parameterized_monodispersed_mie mono_mie(m_host,m_sphere,500e-9);
    log_normal_distribution sd{log(r),0.1};
    polydispersed_mie poly_mie(mono_mie,sd);
    double p = 0.01;
    poly_mie.percentage_accuracy(p);
    check_small(poly_mie.absorption_cross_section());
    check_fast(10 * cpu_duration());

  } end_test_case()
}
