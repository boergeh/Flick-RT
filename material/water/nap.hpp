#ifndef flick_material_nap
#define flick_material_nap

#include "../material.hpp"

namespace flick {
namespace material {
  class nap : public base
  // Nonalgal particulate matter 
  {
    double mass_concentration_; // [kg/m^3]
    double asymmetry_factor_ = 0.8;
    flick::fournier_forand ff_{asymmetry_factor_};
  public:
    nap(double mass_concentration = 1e-3)
      : mass_concentration_{mass_concentration} {
    }
    void mass_concentration(double c)
    // Dry mass concentration [kg/m^3]
    {  
      mass_concentration_ = c;
    } 
    double absorption_coefficient() const {     
      return a_m_star_443() * mass_concentration_ * abs_wavelength_factor();
    }
    double scattering_coefficient() const {
      return b_m_star_555() * mass_concentration_ * scat_wavelength_factor(); 
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      mueller m;
      return m.add(0,0,ff_.value(scattering_direction.mu()));
    }
  private:
    double a_m_star_443() const
    // Mass-specific absorption coefficient at 443 nm for non-algal
    // particles suspendend in water [m^2/kg]
    // See Babin et al. 2003.  
    {
      return 0.035 * 1e3; 
    }   
    double b_m_star_555() const
    // Mass-specific scattering coefficient at 555 nm for non-algal
    // particles suspendend in water [m^2/kg]
    //
    // Babin, M., Morel, A., Fournier-Sicre, V., Fell, F. and
    // Stramski, D., 2003. Light scattering properties of marine
    // particles in coastal and open ocean waters as related to the
    // particle mass concentration. Limnology and oceanography, 48(2),
    // pp.843-859.
    {
      return 0.5 * 1e3; // [m^2/kg]
    }
    double abs_wavelength_factor() const  
    // Babin, M., Stramski, D., Ferrari, G.M., Claustre, H., Bricaud,
    // A., Obolensky, G. and Hoepffner, N., 2003. Variations in the
    // light absorption coefficients of phytoplankton, nonalgal
    // particles, and dissolved organic matter in coastal waters
    // around Europe. Journal of Geophysical Research: Oceans,
    // 108(C7).   
    {
      double slope_per_nm = 0.0123;
      double wl_nm = wavelength()*1e9;
      double wl0_nm = 443; 
      return 0.75 * exp(-slope_per_nm * (wl_nm - wl0_nm));
    }
    double scat_wavelength_factor() const
    // Babin, M., Morel, A., Fournier-Sicre, V., Fell, F. and Stramski,
    // D., 2003. Light scattering properties of marine particles in
    // coastal and open ocean waters as related to the particle mass
    // concentration. Limnology and oceanography, 48(2), pp.843-859.
    {
      double gamma = 0.5;
      return pow(wavelength()/555e-9, -gamma);
    }    
  };
}
}

#endif
