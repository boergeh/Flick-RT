#include "normalized_scattering_matrix_fit.hpp"
#include "../mie/polydispersed_mie.hpp"
#include "../environment/input_output.hpp"

namespace flick {  
  begin_test_case(normalized_scattering_matrix_fit_test) {
    stdcomplex m_host = 1.33;
    stdcomplex m_sphere = 1.38;
    double wl = 500e-9;
    double r = 0.2e-6;
    monodispersed_mie mono_mie(m_host,m_sphere,wl);
    mono_mie.angles(range(0,constants::pi,256).linspace());
    log_normal_distribution sd{log(r),0.2};
    polydispersed_mie poly_mie(mono_mie,sd);
    poly_mie.percentage_accuracy(1);
    auto [a, b, x] = poly_mie.ab_functions();
    size_t n_terms = 10;
    normalized_scattering_matrix_fit m(a,b,x,n_terms);
    write(n_terms,"material/ab_functions/n_terms.txt");
    write(m.scattering_scaling_factor(),"material/ab_functions/scaling_factor.txt");
    double p = 1;
    for (size_t i=0; i<4; ++i) {
      write(pl_function{x,a[i]},"material/ab_functions/a"+std::to_string(i)+".txt",15);
      //write(pl_function{x,a[i]},"material/ab_functions/a_scaled"+std::to_string(i)+".txt",15);
      write(pl_function{x,m.fitted_a(i)},"material/ab_functions/a_fitted"+std::to_string(i)+".txt",15);
      check_small(vec::rms(a[i]-m.fitted_a(i)),p,"a"+std::to_string(i));
    }
    for (size_t i=0; i<2; ++i) {
      write(pl_function{x,b[i]},"material/ab_functions/b"+std::to_string(i)+".txt",15);
      //write(pl_function{x,b[i]},"material/ab_functions/b_scaled"+std::to_string(i)+".txt",15);
      write(pl_function{x,m.fitted_b(i)},"material/ab_functions/b_fitted"+std::to_string(i)+".txt",15);
      check_small(vec::rms(b[i]-m.fitted_b(i)),p,"b"+std::to_string(i));
    }
  } end_test_case()
}