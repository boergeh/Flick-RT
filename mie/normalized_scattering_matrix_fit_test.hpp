#include "normalized_scattering_matrix_fit.hpp"
#include "polydispersed_mie.hpp"
#include "../environment/input_output.hpp"

namespace flick {  
  begin_test_case(normalized_scattering_matrix_fit_test) {
    stdcomplex m_host = 1.0;
    stdcomplex m_sphere = 1.65 + 1e-5i;
    double wl = 355e-9;
    double r = 2e-7;
    monodispersed_mie mono_mie(m_host,m_sphere,wl);
    mono_mie.quadrature_angles(256);
    log_normal_distribution sd{log(r),0.0};
    polydispersed_mie poly_mie(mono_mie,sd);
    auto [a, b, x] = poly_mie.ab_functions();
    normalized_scattering_matrix_fit m_fit(a,b,x,8);
    double p = 1;
    for (size_t i=0; i<4; ++i) {
      write(pl_function{x,a[i]},"mie/ab_functions/a"+std::to_string(i)+".txt",15);
      write(pl_function{x,m_fit.scaled_a(i)},"mie/ab_functions/a_scaled"+std::to_string(i)+".txt",15);
      write(pl_function{x,m_fit.fitted_scaled_a(i)},"mie/ab_functions/a_scaled_fitted"+std::to_string(i)+".txt",15);
      check_small(vec::rms(m_fit.scaled_a(i)-m_fit.fitted_scaled_a(i)),p,"a"+std::to_string(i));
    }
    for (size_t i=0; i<2; ++i) {
      write(pl_function{x,b[i]},"mie/ab_functions/b"+std::to_string(i)+".txt",15);
      write(pl_function{x,m_fit.scaled_b(i)},"mie/ab_functions/b_scaled"+std::to_string(i)+".txt",15);
      write(pl_function{x,m_fit.fitted_scaled_b(i)},"mie/ab_functions/b_scaled_fitted"+std::to_string(i)+".txt",15);
      check_small(vec::rms(m_fit.scaled_b(i)-m_fit.fitted_scaled_b(i)),p,"b"+std::to_string(i));
    }
  } end_test_case()
}
