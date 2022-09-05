#ifndef flick_material_pure_water
#define flick_material_pure_water

#include "../material.hpp"
#include "../../polarization/rayleigh_mueller.hpp"

namespace flick {
namespace material {
  class pure_water : public base {
    pp_function absorption_coefficient_;
    pp_function refractive_index_;
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
     refractive_index_ = read<pp_function>
       (path_+"/refractive_index.txt"); 
     temperature_correction_ = read<pl_function>
       (path_+"/temperature_correction.txt"); 
     salinity_correction_ = read<pl_function>
       (path_+"/salinity_correction.txt"); 
     absorption_coefficient_.add_extrapolation_points(1);
     refractive_index_.add_extrapolation_points(1);
     temperature_correction_.add_extrapolation_points(0);
     salinity_correction_.add_extrapolation_points(0);
    }
    void salinity(const pl_function& s) {
      salinity_psu_ = s;
    }
    void temperature(const pl_function& temperature) {
      temperature_ = temperature;
    }
    double absorption_coefficient() {
      double T = temperature_.value(pose().position().z());
      double S = salinity_psu_.value(pose().position().z());
      double delta_T = T - pope_fry_temperature_;
      double delta_S = S;
      double da_dT = temperature_correction_.value(wavelength());
      double da_dS = salinity_correction_.value(wavelength());
      double a0 = absorption_coefficient_.value(wavelength());
      return a0 + da_dT * delta_T + da_dS * delta_S;
    }
    double scattering_coefficient()
    // Consider update with temp. corr according to: Pure water
    // spectral absorption, scattering, and real part of refractive
    // index model. Algorithm Technical Basis Document
    {
      return pow(129.0 / (wavelength() * 1e9), 4.32);
    }
    mueller mueller_matrix(const unit_vector& scattering_direction)
    // Depolarization ratio from R. S. Farinato and R. L. Rowell, “New
    // values of the light scattering depolarization and anisotropy of
    // water,” J. Chem. Phys. 65, 593–595 (1976).
    {
      return rayleigh_mueller(angle(scattering_direction), 0.039);
    }
    double refractive_index()
    // Consider update to include temp. and sal. corr.
    {
      return refractive_index_.value(wavelength());
    }
  };
}
}

#endif
