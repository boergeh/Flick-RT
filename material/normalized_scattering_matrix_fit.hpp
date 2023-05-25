#ifndef flick_normalized_scattering_matrix_fit
#define flick_normalized_scattering_matrix_fit

#include "../numeric/wigner/wigner_fit.hpp"
#include <armadillo>
 
namespace flick {
  class normalized_scattering_matrix_fit
  /* For definition of the normalized scattering matrix, see e.g.:
   Mishchenko, M.I. and Yang, P., 2018. Far-field Lorenz–Mie
   scattering in an absorbing host medium: theoretical formalism and
   FORTRAN program. Journal of Quantitative Spectroscopy and Radiative
   Transfer, 205, pp.241-252. */
  {
    std::vector<stdvector> a_;
    std::vector<stdvector> b_;
    stdvector x_;
    std::vector<stdvector> alpha_;
    std::vector<stdvector> beta_;
  public:
    normalized_scattering_matrix_fit(const std::vector<stdvector>& a,
				     const std::vector<stdvector>& b,
				     const stdvector& x,
				     size_t n_terms)
      : a_{a}, b_{b}, x_{x}, alpha_(4), beta_(2) { 
      assert(a_.size()==4 && b_.size()==2);
      pe_function a0(x_,a_[0]);
      wigner_fit wf_a0(a0,0,0,n_terms,fit::relative);
      alpha_[0] = wf_a0.coefficients();
      pe_function scaling_function{x_,1/a_[0]};
      stdvector alpha2p3 = wigner_fit(pl_function{x_,a_[1]+a_[2]},2,2,n_terms,fit::scaling,
				      scaling_function).coefficients();
      stdvector alpha2m3 = wigner_fit(pl_function{x_,a_[1]-a_[2]},2,-2,n_terms,fit::scaling,
				      scaling_function).coefficients();
      alpha_[1] = 0.5*(alpha2p3+alpha2m3);
      alpha_[2] = 0.5*(alpha2p3-alpha2m3);
      alpha_[3] = wigner_fit(pl_function{x_,a_[3]},0,0,n_terms,fit::scaling,
				      scaling_function).coefficients();  
      beta_[0] = -1*wigner_fit(pl_function{x_,b_[0]},0,2,n_terms,fit::scaling,
				      scaling_function).coefficients();
      beta_[1] = -1*wigner_fit(pl_function{x_,b_[1]},0,2,n_terms,fit::scaling,
				      scaling_function).coefficients();
    }
    const std::vector<stdvector>& alpha() const {
      return alpha_;
    }
    const std::vector<stdvector>& beta() const {
      return beta_;
    }    
    const stdvector& alpha(size_t i) const {
      return alpha_.at(i);
    }
    const stdvector& beta(size_t i) const {
      return beta_.at(i);
    } 
    stdvector fitted_a(size_t i) const {
      if (i==0)
	return wigner_evaluate(alpha_.at(0),x_,0,0);
      if (i==1 or i==2) {
	stdvector a2p3 = wigner_evaluate(alpha_[1]+alpha_[2],x_,2,2);
	stdvector a2m3 = wigner_evaluate(alpha_[1]-alpha_[2],x_,2,-2);
	if (i==1)
	  return 0.5*(a2p3+a2m3);
	if (i==2)
	  return 0.5*(a2p3-a2m3);
      }
      if (i==3)
	return wigner_evaluate(alpha_.at(i),x_,0,0);
      else
	return {0}; 
    }
    stdvector fitted_b(size_t i) const {
      return wigner_evaluate(-1*beta_.at(i),x_,0,2);
    }
    double scattering_scaling_factor() const {
      return alpha_[0][0];
    }
  };
}

#endif
