#ifndef flick_material_absorption_smoother
#define flick_material_absorption_smoother

#include "air.hpp"
#include "../../numeric/value_collection.hpp"
#include "../../numeric/distribution.hpp"

namespace flick {
  template<class Source_spectrum>
  class absorption_smoother
  // Smooths atmospheric absorption spectra, assuming no scattering
  {
    atmospheric_state atm_;
    const double atmosphere_thickness_ = 100e3;
    std::string gas_;
    Source_spectrum source_spectrum_;
    std::shared_ptr<material::hitran_air> air_all_;
    std::shared_ptr<material::hitran_air> air_rest_;
    bool print_progress_ = false;
  public:
    absorption_smoother(Source_spectrum s, atmospheric_state a,
			const std::string& gas)
      : source_spectrum_{s}, atm_{a}, gas_{gas} {
      air_all_ = std::make_shared<material::hitran_air>(atm_);
      atm_.remove_gas(gas_);
      air_rest_ = std::make_shared<material::hitran_air>(atm_);
      unit_vector d{0,0};
      vector p{0,0,0};
      air_all_->set_direction(d);
      air_all_->set_position(p);
      air_rest_->set_direction(d);
      air_rest_->set_position(p);
    }
    void print_progress(bool tf) {
      print_progress_ = tf;
    }
    double cross_section(const distribution::basic_distribution& d,
			 double accuracy=0.01) {
      value_collection collection(accuracy);
      collection.set_log_values(false);
      collection.noise_floor(1e-33);
      collection.initial_set(2);
      size_t n = 128;
      double max_duration = 60;
      double duration = 0;
      if (print_progress_) {
	std::cout << "Center wavelenth: "<< std::setprecision(5)
		  << vec::mean(d.quantiles(2))*1e9 << " nm"<< std::endl;
	std::cout << std::setprecision(6);
      }
      while(not collection.last_accurate() and duration < max_duration) {
	auto time0 = std::chrono::system_clock::now();
	stdvector wls = d.quantiles(n);
	double T_rest = transmittance(wls,d,*air_rest_);
	double T_all = transmittance(wls,d,*air_all_);
	double optical_depth = -log(T_all/T_rest);
	double cross_sec = optical_depth / atm_.per_area(gas_);
	collection.add(cross_sec,n);
	n *= 2;
	auto time1 = std::chrono::system_clock::now();
	duration = std::chrono::duration<double>(time1-time0).count();
	if (print_progress_) {
	  std::cout << collection << std::endl;
	  if (duration >= max_duration and not collection.last_accurate())
	    std::cout << "timed out after "<<max_duration<<" s" << std::endl;
	}
      }
      return collection.last_value();
    }
  private:
    double transmittance(const stdvector& wls, const distribution::basic_distribution& d, material::hitran_air& a) {
      pl_function f;
      for (const auto& wl:wls) {
	double R = d.pdf(wl);
	double F = source_spectrum_.value(wl);
	a.set_wavelength(wl);
	double tau = a.absorption_optical_depth(atmosphere_thickness_);
	f.append({wl,R*F*exp(-tau)});
      }
      return f.integral();
    }
  };
  
  template<class Source_spectrum>
  stdvector gaussian_smooth_cross_section(absorption_smoother<Source_spectrum>& as, const stdvector& center_wls, double band_width_sigma, double accuracy) {
    stdvector c;
    for (size_t i=0; i<center_wls.size(); i++) {
      double mu = log(center_wls[i]);
      double sig = band_width_sigma;
      c.push_back(as.cross_section(distribution::log_normal(mu,sig),accuracy));
    }
    return c;
  }
  
  template<class Source_spectrum>
  stdvector triangular_smooth_cross_section(absorption_smoother<Source_spectrum>& as, const stdvector& center_wls, double band_width_sigma, double accuracy) {
    stdvector c;
    for (size_t i=0; i<center_wls.size(); i++) {
      double mu = center_wls[i];
      double sig = band_width_sigma;
      double x_low = mu*exp(-sig);
      double x_high = mu*exp(sig);
      double x_mode = mu;
      c.push_back(as.cross_section(distribution::triangular(x_low,x_high,x_mode),accuracy));
    }
    return c;
  }
}

#endif
