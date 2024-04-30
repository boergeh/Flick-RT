#ifndef flick_material_ab_functions
#define flick_material_ab_functions

#include "material.hpp"
#include "normalized_scattering_matrix_fit.hpp"
#include "../numeric/wigner/wigner_fit.hpp"

namespace flick {
namespace material {
  template<class Material>
  std::tuple<std::vector<std::vector<double>>,
	     std::vector<std::vector<double>>,
	     std::vector<double>> mueller_ab_functions(const Material& mat,
						       size_t n_angles)
  // a and b mueller matrix elements, where scattering angle is
  // taken relative to current traveling direction, keeping azimuth
  // angle equal to zero
  {
    std::vector<double> x = wigner_x_values(n_angles);
    std::vector<std::vector<double>> a(4,std::vector<double>(x.size()));
    std::vector<std::vector<double>> b(2,std::vector<double>(x.size()));
    for (size_t i=0; i<n_angles; ++i) {
      x[i] = std::clamp<double>(x[i],-1,1);
      mueller m = mat.mueller_matrix(unit_vector{acos(x[i]),0});
      a[0][i] = m.value(0,0);
      a[1][i] = m.value(1,1);
      a[2][i] = m.value(2,2);
      a[3][i] = m.value(3,3);
      b[0][i] = m.value(0,1);
      b[1][i] = m.value(2,3);
    }
    return {a,b,x};
  }

  template<class Material>
  std::tuple<std::vector<std::vector<double>>,
	     std::vector<std::vector<double>>> fitted_mueller_alpha_beta(const Material& mat, size_t n_terms, std::optional<size_t> n_points = std::optional<size_t>(), fit phase_function_fit = fit::relative)
  // Wigner d-function fits to the non-zero Mueller matrix a and b
  // functions. Returned coefficients are normalized such the 4*pi
  // solid angle integral over a[0] (phase function) equals one as
  // is common practice in Flick.
  {
    if (not n_points)
      n_points = wigner_n_sampling_points(n_terms);
    auto [a,b,x] = mueller_ab_functions(mat, *n_points);
    normalized_scattering_matrix_fit fit(a,b,x,n_terms,phase_function_fit);
    return {fit.alpha(), fit.beta()};
  }
  
  template<class Material>
  std::tuple<std::vector<std::vector<double>>,
	     std::vector<std::vector<double>>,
	     std::vector<double>> fitted_mueller_ab_functions(const Material& mat, size_t n_terms, std::optional<size_t> n_points = std::optional<size_t>(), fit phase_function_fit = fit::relative)
  {
    if (not n_points)
      n_points = wigner_n_sampling_points(n_terms);
    auto [a,b,x] = mueller_ab_functions(mat, *n_points);
    normalized_scattering_matrix_fit m(a,b,x,n_terms,phase_function_fit);
    std::vector<std::vector<double>> a_fitted(4);
    for (size_t i = 0; i < a_fitted.size(); ++i) {
      a_fitted[i] = m.fitted_a(i);
    }
    std::vector<std::vector<double>> b_fitted(2);
    for (size_t i = 0; i < b_fitted.size(); ++i)
      b_fitted[i] = m.fitted_b(i);
    return {a_fitted, b_fitted, x};
  }
}
}

#endif
