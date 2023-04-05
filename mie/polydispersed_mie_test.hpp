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
		2*3.14159*pow(r,2),0.1,"a");

    check(mono_mie.scattering_matrix_element(0,0)[0] >
	  mono_mie.scattering_matrix_element(0,0)[1],"b");

    polydispersed_mie poly_mie(mono_mie,distribution);
    poly_mie.precision(3);
    check_close(poly_mie.absorption_cross_section(),
    		mono_mie.absorption_cross_section(),0.02,"c");
    check_close(poly_mie.scattering_cross_section(),
    		mono_mie.scattering_cross_section(),0.02,"d");
    check_close(poly_mie.scattering_matrix_element(0,0)[0],
    		mono_mie.scattering_matrix_element(0,0)[0],0.02,"e");

    check_close(distribution.particles_per_volume(1),
		1/(4./3*constants::pi*pow(r,3)),0.1,"f");
  } end_test_case()

  begin_test_case(poly_mie_test_B) {
    stdcomplex m_host = {1,0};
    stdcomplex m_sphere = {1.3,0.0};
    double r = 100e-6;
    parameterized_monodispersed_mie mono_mie(m_host,m_sphere,500e-9);
    log_normal_distribution sd{log(r),0.1};
    polydispersed_mie poly_mie(mono_mie,sd);
    poly_mie.precision(3);
    check_close(poly_mie.scattering_efficiency(),2,0.01);
  } end_test_case()

   begin_test_case(poly_mie_test_C) {
    double pi = constants::pi;
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.65 + 1e-5i;
    double wl = 355e-9;
    double r = 1e-6;
    monodispersed_mie mono_mie(m_host,m_sphere,wl);
    log_normal_distribution sd{log(r),0.5};
    polydispersed_mie poly_mie(mono_mie,sd);
    poly_mie.precision(2);
    double Cscat = poly_mie.scattering_cross_section();
    //std::cout << Cscat << std::endl;
    //std::cout << poly_mie.xy_points().integral()/Cscat << std::endl;
    write<pl_function>(poly_mie.xy_points(),"mie/xy_points.txt",9);
  } end_test_case()

}
