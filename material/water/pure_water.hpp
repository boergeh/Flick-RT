#ifndef flick_material_pure_water
#define flick_material_pure_water

#include "../material.hpp"
#include "../../polarization/rayleigh_mueller.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
  class zhang_hu_scattering
  // Zhang, X. and Hu, L., 2009. Estimating scattering of pure water
  // from density fluctuation of the refractive index. Optics Express,
  // 17(3), pp.1671-1678.
  {
    double S_;
    double T_;
    double wavelength_;
  public:
    zhang_hu_scattering(double S, double T, double wavelength)
      : S_{S}, T_{T}, wavelength_{wavelength} {
    }
    static double depolarization_ratio() {
      // R. S. Farinato and R. L. Rowell, “New values of the light
      // scattering depolarization and anisotropy of water,”
      // J. Chem. Phys. 65, 593–595 (1976).
      return 0.039;
    }
    double real_refractive_index() const
    // Quan, X. and Fry, E.S., 1995. Empirical equation for the
    // index of refraction of seawater. Applied optics, 34(18),
    // pp.3477-3480.
    {
      std::vector<double> n = {1.31405, 1.31405e-4, -1.05e-6, 1.6e-8, -2.02e-6, 15.868,
	0.01155, -0.00423, -4382, 1.1455e6};
      double Tc = constants::to_celsius(T_);
      double lambda = wavelength_*1e9; // [nm]
      double n_sw = n[0]+(n[1]+n[2]*Tc+n[3]*pow(Tc,2))*S_+n[4]*pow(Tc,2)+
	(n[5]+n[6]*S_+n[7]*Tc)/lambda+n[8]/pow(lambda,2)+n[9]/pow(lambda,3);
      return n_sw * air_refractive_index();
    }
    double volume_scattering_90() const {
      using namespace constants;
      double lambda = wavelength_;
      double beta_T = isothermal_compressibility();
      double f = cabannes_factor();
      return pow(pi,2)/(2*pow(lambda,4))*pow(PMH_factor(),2)*k_B*T_*beta_T*f;
    }
    double scattering_coefficient() const {
      double beta_90 = volume_scattering_90();
      double d = depolarization_ratio();
      double k = (1-d)/(1+d);
      return (1+k/3)*beta_90*4*constants::pi;
    }
  private:
    double cabannes_factor() const {
      double d = depolarization_ratio();
      return (6+6*d)/(6-7*d);
    }
    double PMH_factor() const {
      // Proutiere, A., Megnassan, E. and Hucteau, H.,
      // 1992. Refractive index and density variations in pure
      // liquids: A new theoretical relation. The Journal of Physical
      // Chemistry, 96(8), pp.3485-3489.
      double n = real_refractive_index();
      double n2 = pow(n,2);
      return (n2-1)*(1+2./3*(n2+2)*pow((n2-1)/(3*n),2));
    }
    double air_refractive_index() const
    // Ciddor, P.E., 1996. Refractive index of air: new equations for
    // the visible and near infrared. Applied optics, 35(9),
    // pp.1566-1573.
    {
      double nu = 1/(wavelength_*1e6); // [per microns]
      std::vector<double> k = {238.0185, 5792105, 57.362, 167917}; // [microns^-2]
      return 1 + (k[1]/(k[0]-pow(nu,2))+k[3]/(k[2]-pow(nu,2)))/1e8;
    }
    /*
    double scattering_slope() const {
      // Morel, A., 1974. Optical properties of pure water and pure
      // seawater. Optical aspects of oceanography.
      return -4;//.32;
    }
    */
    double isothermal_compressibility() const {
      // Kell, G.S., 1970. Isothermal compressibility of liquid water
      // at 1 atm. Journal of Chemical and Engineering Data, 15(1),
      // pp.119-122.
      double Tc = constants::to_celsius(T_);
      std::vector<double> a = {50.88630, 0.7171582, 0.7819867e-3, 31.62214e-6,
	-0.1323594e-6, 0.6345750e-9};
      std::vector<double> b = {1, 21.65928e-3};
      return sum_terms(a,Tc)/sum_terms(b,Tc)*1e-6/constants::P_stp;
    }
    double sum_terms(const std::vector<double>& coefs, double Tc) const {
      double v = 0;
      for (size_t i = 0; i < coefs.size(); i++) {
	v += coefs[i]*pow(Tc,i);
      }
      return v;
    }
  };
  
  class pure_water : public base {
    double volume_fraction_ = 1;
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
    pure_water(double salinity_psu, double temperature, double volume_fraction = 1)
      : pure_water() {
      salinity_psu_ = salinity_psu;
      temperature_ = temperature;
      volume_fraction_ = volume_fraction;
    }   
    void salinity(const pl_function& s) {
      salinity_psu_ = s;
    }
    void temperature(const pl_function& temperature) {
      temperature_ = temperature;
    }
    double absorption_coefficient() const {
      double delta_T = temperature() - pope_fry_temperature_;
      double delta_S = salinity();
      double da_dT = temperature_correction_.value(wavelength());
      double da_dS = salinity_correction_.value(wavelength());
      double a0 = absorption_coefficient_.value(wavelength());
      return (a0 + da_dT * delta_T + da_dS * delta_S) * volume_fraction_;
    }
    double scattering_coefficient() const {
      zhang_hu_scattering s(salinity(),temperature(),wavelength());
      return s.scattering_coefficient() * volume_fraction_;
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      return rayleigh_mueller(angle(scattering_direction),
			      zhang_hu_scattering::depolarization_ratio());
    }
    double real_refractive_index() const {
      double wl_1 = 280e-9;
      double wl_2 = 1600e-9;
      double wl = wavelength();
      if (wl < wl_1)
	return segelstein_real_refractive_index_.value(wl) + sal_temp_shift(wl_1);
      else if (wl > wl_2)
	return segelstein_real_refractive_index_.value(wl) + sal_temp_shift(wl_2);
      else {
	zhang_hu_scattering s(salinity(),temperature(),wavelength());
	return s.real_refractive_index();
      }
    }
  private:
    double salinity() const {
      return salinity_psu_.value(pose().position().z());
    }
    double temperature() const {
      return temperature_.value(pose().position().z());
    }
    double sal_temp_shift(double wl) const {
      zhang_hu_scattering s(salinity(),temperature(),wl);
      return s.real_refractive_index()-segelstein_real_refractive_index_.value(wl);
    } 
  };
}
}

#endif
