#ifndef flick_material_phytoplankton
#define flick_material_phytoplankton

#include "../material.hpp"
#include "../../numeric/flist.hpp"

namespace flick {
namespace material {
  class phytoplankton : public base {
    double asymmetry_factor_ = 0.98;
    double chl_concentration_ = 1e-6; // [kg/m^3]    
    //flick::henyey_greenstein hg_{asymmetry_factor_};
    flick::fournier_forand ff_{asymmetry_factor_};
    pl_function A_;
    pl_function B_;
  public:
    phytoplankton(double chl_concentration = 1e-6)
      : chl_concentration_{chl_concentration} {
      const std::string path_{"/material/water"};
      pl_flist AB = read<pl_flist>(path_+"/phytoplankton.txt");
      A_ = AB(0);
      B_ = AB(1);
      A_.add_extrapolation_points();
      B_.add_extrapolation_points();
    }
    void chl_concentration(double c) {
      chl_concentration_ = c;
    } 
    double absorption_coefficient() const {
      return a_star() * chl_mg();
    }
    double scattering_coefficient() const {
      return b_star() * chl_mg() * b_wavelength_factor(); 
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      return m.add(0,0,ff_.value(scattering_direction.mu()));
    }
  private:
    double a_star() const
    // Chlorophyll-specific absorption coefficient [m^2/mg]
    // Bricaud, A., Babin, M., Morel, A. and Claustre, H.,
    // 1995. Variability in the chlorophyll‐specific absorption
    // coefficients of natural phytoplankton: Analysis and
    // parameterization. Journal of Geophysical Research: Oceans,
    // 100(C7), pp.13321-13332.
    {
      double wl_nm = wavelength()*1e9;
      return A_.value(wl_nm) * pow(chl_mg(), -B_.value(wl_nm));
    }
    double b_star() const
    // Chlorophyll-specific scattering coefficient [m^2/mg]     
    // Morel, A., 1987. Chlorophyll-specific scattering coefficient of
    // phytoplankton. A simplified theoretical approach. Deep Sea
    // Research Part A. Oceanographic Research Papers, 34(7),
    // pp.1093-1105.
    {
      return 0.3 * pow(chl_mg(), -0.38);
    }
    double b_wavelength_factor() const
    // Babin, M., Morel, A., Fournier-Sicre, V., Fell, F. and Stramski,
    // D., 2003. Light scattering properties of marine particles in
    // coastal and open ocean waters as related to the particle mass
    // concentration. Limnology and oceanography, 48(2), pp.843-859.
    {
      double gamma = 0.5;
      return pow(wavelength()/550e-9, -gamma);
    }    
    double chl_mg() const {
      return chl_concentration_ * 1e6;
    }
  };
}
}

#endif
