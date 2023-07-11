#ifndef flick_parameterized_monodispersed_mie
#define flick_parameterized_monodispersed_mie

#include "basic_monodispersed_mie.hpp"
#include "../numeric/physics_function.hpp"
#include "../environment/input_output.hpp"

namespace flick { 
  class parameterized_monodispersed_mie : public basic_monodispersed_mie
  // Approximate Mie-code solutions for large spheres, see: Stamnes,
  // K., Hamre, B., Stamnes, J.J., Ryzhikov, G., Biryulina, M.,
  // Mahoney, R., Hauss, B. and Sei, A., 2011. Modeling of radiation
  // transport in coupled atmosphere-snow-ice-ocean systems. Journal
  // of Quantitative Spectroscopy and Radiative Transfer, 112(4),
  // pp.714-726. Some modifications to be consistent with full Mie
  // theory, which gives negative absorption for transparent particles
  // in an absorbing host medium - gives correctly less total
  // absorption when e.g. bubbles are added in water.
  {
    double Qa_{0};
    double n_ = real(m_sphere_ / m_host_);
    pl_function g0_ = read<pl_function>("mie/g_parameterized.txt");
    
    void update_efficiency() {
      stdcomplex n = m_sphere_ / m_host_;
      stdcomplex arg = 1./n * (pow(n,3) - pow(pow(n,2)-1., 3./2));
      double Qa0 = 8./3*imag(m_sphere_-m_host_) * real(size_parameter_in_host())
      	* std::abs(arg);
      Qa_ = 0.94 * tanh(Qa0/0.94);
    }
    double geometrical_cross_section() const {
      return pi_ * pow(radius_,2);
    }
    double asymmetry_factor() const {
      return pow(g0_.value(n_), pow(1-Qa_, 0.6));
    }
    double extinction_cross_section() const {
      return 2 * geometrical_cross_section();
    }
  public:
    using basic_monodispersed_mie::basic_monodispersed_mie;
    void radius(double r) {
      radius_ = r;
      update_efficiency();
    }
    void angles(const stdvector& angles) {
      angles_ = angles;
    }
    double absorption_cross_section() const {
      return Qa_ * geometrical_cross_section();
    }
    double scattering_cross_section() const {
      return extinction_cross_section() - absorption_cross_section();
    }
    stdvector scattering_matrix_element(size_t row, size_t col) const
    // Note that integratinig element 0,0 over all 4*pi solid angles gives
    // the scattering cross section, where we count from 0 instead of one.
    {
      if (row == 0 && col == 0) {
	stdvector hg(angles_.size());
	double g = asymmetry_factor();
	for (size_t i=0; i<hg.size(); ++i) {
	  hg[i] = henyey_greenstein(g).phase_function(angles_[i]);
	}
	return hg * scattering_cross_section();
      }
      else
	return stdvector(angles_.size(),0);
    }    
  };
}

#endif
  

