#ifndef flick_basic_monodispersed
#define flick_basic_monodispersed

#include "../numeric/constants.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/function.hpp"

namespace flick {
  class basic_monodispersed_mie {
  protected:
    const double pi_{constants::pi};

    stdcomplex m_host_;
    stdcomplex m_sphere_;
    double vacuum_wl_;
    double radius_{1e-6};
    stdvector angles_{0};
    stdcomplex wavenumber_in_host_{2*pi_*m_host_/vacuum_wl_};
    stdcomplex m_sphere_at_r0_;
    
    stdcomplex size_parameter_in_host() {
      return wavenumber_in_host_ * radius_;
    }
    double size_parameter_in_vacuum() {
      return 2*pi_*radius_/vacuum_wl_;
    }
  public:
    basic_monodispersed_mie(const stdcomplex& m_host,
			     const stdcomplex& m_sphere,
			     double vacuum_wl)
      : m_host_{m_host}, m_sphere_{m_sphere}, vacuum_wl_{vacuum_wl},
	m_sphere_at_r0_{m_sphere} {
    }
    virtual void radius(double r) = 0;
    virtual void angles(const stdvector& angles) = 0;
    virtual double absorption_cross_section() const = 0;
    virtual double scattering_cross_section() const = 0;
    virtual stdvector scattering_matrix_element(size_t row,
						size_t col) const = 0;

    double radius() const {
      return radius_;
    }
  };  
}

#endif
