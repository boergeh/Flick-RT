#ifndef flick_material_pure_water
#define flick_material_pure_water

#include "../material.hpp"
#include "../../polarization/rayleigh_mueller.hpp"

namespace flick {
namespace material {
  class pure_water : public base {
    pp_function absorption_coefficient_;
    pp_function segelstein_real_refractive_index_;
    pl_function temperature_correction_;
    pl_function salinity_correction_;
    pl_function salinity_psu_{0};
    pl_function temperature_{constants::T_ntp};
    const double pope_fry_temperature_{295};
    const std::string path_{"/material/water"};
  public:
    pure_water() {
     absorption_coefficient_ = read<pp_function>
       (path_+"/absorption_coefficient.txt"); 
     segelstein_real_refractive_index_ = read<pp_function>
       (path_+"/refractive_index.txt"); 
     temperature_correction_ = read<pl_function>
       (path_+"/temperature_correction.txt"); 
     salinity_correction_ = read<pl_function>
       (path_+"/salinity_correction.txt"); 
     absorption_coefficient_.add_extrapolation_points(1);
     segelstein_real_refractive_index_.add_extrapolation_points(1);
     temperature_correction_.add_extrapolation_points(0);
     salinity_correction_.add_extrapolation_points(0);
    }
    pure_water(double salinity_psu, double temperature)
      : pure_water()  {
      salinity_psu_ = salinity_psu;
      temperature_ = temperature;
    }  
    void salinity(const pl_function& s) {
      salinity_psu_ = s;
    }
    void temperature(const pl_function& temperature) {
      temperature_ = temperature;
    }
    double absorption_coefficient() const {
      double T = temperature_.value(pose().position().z());
      double S = salinity_psu_.value(pose().position().z());
      double delta_T = T - pope_fry_temperature_;
      double delta_S = S;
      double da_dT = temperature_correction_.value(wavelength());
      double da_dS = salinity_correction_.value(wavelength());
      double a0 = absorption_coefficient_.value(wavelength());
      return a0 + da_dT * delta_T + da_dS * delta_S;
    }
    double scattering_coefficient() const
    // Consider update with temp. corr according to: Pure water
    // spectral absorption, scattering, and real part of refractive
    // index model. Algorithm Technical Basis Document
    {
      return pow(129.0 / (wavelength() * 1e9), 4.32);
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const
    // Depolarization ratio from R. S. Farinato and R. L. Rowell, “New
    // values of the light scattering depolarization and anisotropy of
    // water,” J. Chem. Phys. 65, 593–595 (1976).
    {
      return rayleigh_mueller(angle(scattering_direction), 0.039);
    }
    double real_refractive_index() const
    // Partly following https://spot.colorado.edu/
    // ~braup/glims/Rottgers-liquid_water-
    // absorption-spectral-spectroscopy-ATBD_waterradiance_watermodel_v2.pdf
    {
      double wl_1 = 280e-9;
      double wl_2 = 1600e-9;
      double wl = wavelength();
      if (wl < wl_1)
	return segelstein_real_refractive_index_.value(wl) + sal_temp_shift(wl_1);
      else if (wl > wl_2)
	return segelstein_real_refractive_index_.value(wl) + sal_temp_shift(wl_2);
      else
	return quan_fry_water_refractive_index(wl);
    }
  private:
    double sal_temp_shift(double wl) const {
      return quan_fry_water_refractive_index(wl) -
	segelstein_real_refractive_index_.value(wl);
    }
      
    double quan_fry_water_refractive_index(double wavelength) const
    // Quan, X. and Fry, E.S., 1995. Empirical equation for the
    // index of refraction of seawater. Applied optics, 34(18),
    // pp.3477-3480.
    {
      std::vector<double> n = {1.31405, 1.31405e-4, -1.05e-6, 1.6e-8, -2.02e-6, 15.868,
	0.01155, -0.00423, -4382, 1.1455e6};
      double T = temperature_.value(pose().position().z())-273.15; // [Celsius]
      double S = salinity_psu_.value(pose().position().z());
      double lambda = wavelength*1e9; // [nm]
      double n_sw = n[0]+(n[1]+n[2]*T+n[3]*pow(T,2))*S+n[4]*pow(T,2)+
	(n[5]+n[6]*S+n[7]*T)/lambda+n[8]/pow(lambda,2)+n[9]/pow(lambda,3);
      return n_sw;// * air_refractive_index();
    }
    /*
    double air_refractive_index() const
    // Zhang, X., and L. Hu. 2009. Estimating scattering of pure water
    // from density fluctuation of the refractive
    // index. Opt. Expr. 17: 1671- 1678.
    {
      double nu = 1/(wavelength()*1e6); // [per microns]
      std::vector<double> k = {238.0185, 5792105, 57.362, 167917}; // [microns^-2]
      return 1 + (k[1]/(k[0]-pow(nu,2))+k[3]/(k[2]-pow(nu,2)))/1e8;
    }
    */
  };
}
}

#endif
