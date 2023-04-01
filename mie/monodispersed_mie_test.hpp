#include "monodispersed_mie.hpp"
#include "parameterized_monodispersed_mie.hpp"

namespace flick {
  begin_test_case(mono_mie_bessel_test_A) {
    stdcomplex z{1 , 0};
    spherical_bessel j(z,3);
    check_close(norm(j.terms()[2]),norm(stdcomplex{0.06203505201137386,
						     0}),1e-13);
  } end_test_case()
  
  begin_test_case(mono_mie_bessel_test_B) {
    double epsilon = 1e-9;
    stdcomplex z_a{2*constants::pi, 0};
    stdcomplex z_b{2*constants::pi+epsilon , 0};
    spherical_bessel j_a(z_a,2);
    spherical_bessel j_b(z_b,2);
    check_close(abs(vec::sum(j_a.terms())),abs(vec::sum(j_b.terms())),1e-4);
  } end_test_case()

  begin_test_case(mono_mie_test_A) {
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

  begin_test_case(mono_mie_test_B) {
    stdcomplex m1 = 1.0+0.05i;
    stdcomplex m2 = 1.53;
    double wl = 2*constants::pi;
    double r = 10/(2*constants::pi)*wl;
    monodispersed_mie mie(m1,m2,wl);
    mie.radius(r);
    auto [a, b] = mie.ab_coefficients();
    double percent = 1e-10; 
    check_close(a[1].real(),0.82786371508743,percent,"a");
    check_close(b[1].imag(),0.91474090929954,percent,"b");
    check_close(a[10].imag(),0.28299838641514,percent,"c");
    check_close(b[10].real(),-0.04440119340674,percent,"d");
    m1 = 1.33 + 0.1i;
    m2 = 1.0;
    r = 2500;
    monodispersed_mie mie2(m1,m2,wl);
    mie2.radius(r);
    auto [a2, b2] = mie2.ab_coefficients();
    check_close(a2[1].real(),4.3914709187499176e216,percent,"e");
  } end_test_case()

   begin_test_case(mono_mie_test_C) {
    double x = 50;
    double m_imag = 0.01;
    stdcomplex m_host = {1.3, m_imag};
    stdcomplex m_sphere = 1.3;
    double wl = 500e-9;
    double r = x * wl / (2*constants::pi);
    monodispersed_mie mie(m_host,m_sphere,wl);
    mie.radius(r);
    double area = constants::pi*pow(r,2);
    double Q_ext = (mie.absorption_cross_section()
		    +mie.scattering_cross_section())/area;
    check_close(Q_ext,-1.99948,1e-3);
  } end_test_case()
  
  begin_test_case(mono_mie_test_D) {
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.3;
    double wl = 500e-9;
    double r = 0.5e-6;
    monodispersed_mie mie(m_host,m_sphere,wl);
    mie.radius(r);
    auto [Cext,Csca] = mie.es_coefficients();
    double area = constants::pi*pow(r,2);
    double Qext = Cext/area;
    double Qsca = Csca/area;
    double Qabs = Qext-Qsca;
    check_small(Qabs,1e-11);
  } end_test_case()
  
  begin_test_case(mono_mie_test_E) {
    double pi = constants::pi;
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.33 + 1e-5i;
    double wl = 500e-9;
    double r = 100e-6;
    stdvector angles = range(0,pi,100).linspace();
    monodispersed_mie mie(m_host,m_sphere,wl);
    mie.radius(r);
    mie.angles(angles);
    parameterized_monodispersed_mie pmie(m_host,m_sphere,wl);
    pmie.radius(r);
    pmie.angles(angles);
    check_close(pmie.absorption_cross_section(),
    		mie.absorption_cross_section(),2.1,"abs");
    check_close(pmie.scattering_cross_section(),
    		mie.scattering_cross_section(),1,"scat");

    auto [S11, S22] = mie.s_functions();
    double k1 = 2*pi*m_host.real()/wl;
    double Cext = 4*pi/k1*imag(S11[0]);
    check_close(Cext,mie.scattering_cross_section()
		+mie.absorption_cross_section(),1e-12);
    
    stdvector f = mie.scattering_matrix_element(0,0)*vec::sin(angles);
    double Cscat = 2*pi*pl_function(angles,f).integral();
    check_close(Cscat, mie.scattering_cross_section(),50);

    f = pmie.scattering_matrix_element(0,0)*vec::sin(angles);
    Cscat = 2*pi*pl_function(angles,f).integral();
    check_close(Cscat, pmie.scattering_cross_section(),0.8);
  } end_test_case()
  
    begin_test_case(mono_mie_test_F) {
    double pi = constants::pi;
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.33 + 0i;
    double wl = 500e-9;
    double r = 1e-10;
    stdvector angles = {0};
    monodispersed_mie mie(m_host,m_sphere,wl);
    mie.radius(r);
    mie.angles(angles);
    check_close(mie.scattering_matrix_element(0,0)[0]
		/mie.scattering_cross_section(),0.119366,1e-3);
    check_close(mie.scattering_matrix_element(3,3)[0]
		/mie.scattering_cross_section(),0.119366,1e-3);
  } end_test_case()
  
    begin_test_case(mono_mie_test_G) {
    stdcomplex m_host = 1.33 + 1e-5i;
    stdcomplex m_sphere = 1;
    double wl = 500e-9;
    double r = 50e-6;
    monodispersed_mie mie(m_host,m_sphere,wl);
    mie.radius(r);
    parameterized_monodispersed_mie pmie(m_host,m_sphere,wl);
    pmie.radius(r);
    check_close(pmie.absorption_cross_section(),
    		mie.absorption_cross_section(),11,"abs");
    check_close(pmie.scattering_cross_section(),
    		mie.scattering_cross_section(),2,"scat");

  } end_test_case()

}
