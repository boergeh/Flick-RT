#ifndef flick_material_pure_water
#define flick_material_pure_water

#include "../material.hpp"
#include "../../polarization/rayleigh_mueller.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
namespace material {
namespace water {
  using vec = std::vector<double>;
  class base {
  protected:
    double S_; // [PSU]
    double T_; // [Kelvin]
    double Tc_; // [Celsius]
  public:
    base(double salinity, double temperature)
      : S_{salinity}, T_{temperature} {
      Tc_ = constants::to_celsius(T_);
    }
  protected:
    double power_sum(const vec& coefficients, double variable) const {
      double v = 0;
      for (size_t i = 0; i < coefficients.size(); i++) {
	v += coefficients[i]*pow(variable,i);
      }
      return v;
    }
  };

  class density : public base
  // [Kg/m^3], from UNESCO, 38, 1981. Adapted from Matlab script by
  // Zhang, X. and Hu, L., 2009.
  {
  public:
    using base::base;
    double value() const {
      return fresh_water() + salt_water_correction(); 
    }
  private:
    double fresh_water() const {
      vec a = {999.842594, 6.793952e-2, -9.09529e-3, 1.001685e-4,
	-1.120083e-6, 6.536332e-9};
      return power_sum(a,Tc_);
    }
    double salt_water_correction() const {
      vec a = {8.24493e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9};
      vec b = { -5.72466e-3, 1.0227e-4, -1.6546e-6};
      vec c = {4.8314e-4};
      return power_sum(a,Tc_)*S_+power_sum(b,Tc_)*pow(S_,1.5)+
	power_sum(c,Tc_)*pow(S_,2);
    }    
  };

  class activity : public base
  // Partial derivative of natural logarithm w.r.t.salinity
  // Millero and Leung (1976), American Journal of Science, 276,
  // 1035-1077. Adapted from Matlab script by Zhang, X. and Hu, L.,
  // 2009.
  {
    vec a = {-5.58651e-4, 2.40452e-7, -3.12165e-9, 2.40808e-11};
    vec b = {1.79613e-5, -9.9422e-8, 2.08919e-9, -1.39872e-11};
    vec c = {-2.31065e-6, -1.37674e-9, -1.93316e-11};
  public:
    using base::base;
    double value() const {
      return power_sum(a,Tc_)+1.5*power_sum(b,Tc_)*pow(S_,0.5)+
	2*power_sum(c,Tc_)*S_;
    }
  };

  class compressibility : public base {
  public:
    using base::base;
    double value() {
      return 1/(fresh_water_coefficient_millero()+salinity_correction_coefficient())/
	constants::P_stp; 
    }
  private:
    double fresh_water_coefficient_kell() const
    // Kell, G.S., 1970. Isothermal compressibility of liquid water at 1
    // atm. Journal of Chemical and Engineering Data, 15(1),
    // pp.119-122. Adapted from Matlab script by Zhang, X. and Hu, L., 
    // 2009.
    {
      vec a = {50.88630, 0.7171582, 0.7819867e-3, 31.62214e-6,
	-0.1323594e-6, 0.6345750e-9};
      vec b = {1, 21.65928e-3};
      double c = power_sum(a,Tc_)/power_sum(b,Tc_)*1e-6;
      return 1/c;
    }
    double fresh_water_coefficient_millero() const
    // Millero (1980), Deep-sea Research
    {
      vec a = {19652.21, 148.4206, -2.327105, 1.360477e-2, -5.155288e-5};
      return power_sum(a,Tc_);
    }
    double salinity_correction_coefficient() const {
      vec a = {54.6746, -0.603459, 1.09987e-2, -6.167e-5};
      vec b = {7.944e-2, 1.6483e-2, -5.3009e-4};
      return (power_sum(a,Tc_)*S_ + power_sum(b,Tc_)*pow(S_,1.5));
    }
  };
  
  class refractive_index : public base
  // Quan, X. and Fry, E.S., 1995. Empirical equation for the
  // index of refraction of seawater. Applied optics, 34(18),
  // pp.3477-3480.
  {
    vec a = {1.31405, 1.779e-4, -1.05e-6, 1.6e-8, -2.02e-6,
      15.868, 0.01155, -0.00423, -4382, 1.1455e6};
  public:
    using base::base;
    double at(double wavelength) const {
      double lambda = wavelength*1e9; // [nm]
      return (a[0]+(a[1]+a[2]*Tc_+a[3]*pow(Tc_,2))*S_+a[4]*pow(Tc_,2)+
	      (a[5]+a[6]*S_+a[7]*Tc_)/lambda+a[8]/pow(lambda,2)+a[9]/pow(lambda,3))*
	air_refractive_index(wavelength);
    }
    double dn_dS(double wavelength) {
      double wl = wavelength*1e9; // [nm]
      return (a[1]+a[2]*Tc_+a[3]*pow(Tc_,2)+a[6]/wl)*
	air_refractive_index(wavelength);
    }
    double density_variation(double wavelength) const {
      // Proutiere, A., Megnassan, E. and Hucteau, H.,
      // 1992. Refractive index and density variations in pure
      // liquids: A new theoretical relation. The Journal of Physical
      // Chemistry, 96(8), pp.3485-3489.
      double n = at(wavelength);
      double n2 = pow(n,2);
      return (n2-1)*(1+2./3*(n2+2)*pow(n/3-1/(3*n),2));
    }  
  private:
    double air_refractive_index(double wavelength) const
    // Ciddor, P.E., 1996. Refractive index of air: new equations for
    // the visible and near infrared. Applied optics, 35(9),
    // pp.1566-1573.
    {
      vec a = {238.0185, 5792105, 57.362, 167917}; // [microns^-2]
      double nu = 1/(wavelength*1e6); // [microns^-1]
      return 1+(a[1]/(a[0]-pow(nu,2))+a[3]/(a[2]-pow(nu,2)))/1e8;
    }      
  };
 
  class scattering : public base
  // Zhang, X. and Hu, L., 2009. Estimating scattering of pure water
  // from density fluctuation of the refractive index. Optics Express,
  // 17(3), pp.1671-1678.
  {
  public:
    using base::base;
    static double depolarization_ratio() {
      // R. S. Farinato and R. L. Rowell, “New values of the light
      // scattering depolarization and anisotropy of water,”
      // J. Chem. Phys. 65, 593–595 (1976).
      return 0.039;
    }
    double vsf90_density(double wavelength) const {
      refractive_index n(S_,T_);
      double wl = wavelength;
      double beta_T = compressibility(S_,T_).value();
      return vsf90_factor(wl) / 2 * pow(n.density_variation(wl),2) *
	constants::k_B * T_ * beta_T;
    }
    double vsf90_salinity(double wavelength) const {
      refractive_index n(S_,T_);
      double wl = wavelength;
      return vsf90_factor(wl) * 2 * pow(n.at(wl),2) * S_ *
	water_molecular_weight() * pow(n.dn_dS(wl),2) /
      	(density(S_,T_).value() * (-activity(S_,T_).value())*constants::N_A);
    }
    double coefficient(double wavelength) const {
      double wl = wavelength;
      double beta90 = vsf90_density(wl) + vsf90_salinity(wl);
      double d = depolarization_ratio();
      double k = (1-d)/(1+d);
      return (1+k/3) *beta90 * 4 * constants::pi;
    }
  private:
    double vsf90_factor(double wavelength) const {
      return pow(constants::pi,2)/pow(wavelength,4)*cabannes_factor();
    }
    double water_molecular_weight() const
    // [kg/mol]
    {
      return 18e-3;
    }
    double cabannes_factor() const {
      double d = depolarization_ratio();
      return (6+6*d)/(6-7*d);
    }
  };
  }
    
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
      absorption_coefficient_.add_constant_extrapolation();
      segelstein_real_refractive_index_.add_constant_extrapolation();
      temperature_correction_.add_zero_extrapolation();
      salinity_correction_.add_zero_extrapolation();
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
      double da_dT = temperature_correction();
      double da_dS = salinity_correction_.value(wavelength());
      double a0 = absorption_coefficient_.value(wavelength());
      return (a0 + da_dT * delta_T + da_dS * delta_S) * volume_fraction_;
    }
    double scattering_coefficient() const {
      return water::scattering(salinity(),temperature()).coefficient(wavelength())*
	volume_fraction_;
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const {
      return rayleigh_mueller(angle(scattering_direction),
			      water::scattering::depolarization_ratio());
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
	water::scattering s(salinity(),temperature());
	return water::refractive_index(salinity(),temperature()).at(wavelength());
      }
    }
  private:
    double temperature_correction() const {
      double noise_cutoff = 480e-9;
      double wl = wavelength();
      if (wl > noise_cutoff)
	return temperature_correction_.value(wl);
      return 0;
    }
    double salinity() const {
      return salinity_psu_.value(pose().position().z());
    }
    double temperature() const {
      return temperature_.value(pose().position().z());
    }
    double sal_temp_shift(double wl) const {
      return water::refractive_index(salinity(),temperature()).at(wavelength())-
	segelstein_real_refractive_index_.value(wl);
    } 
  };
}
}

#endif
