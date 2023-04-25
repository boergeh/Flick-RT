#ifndef flick_filter
#define flick_filter

#include "../environment/input_output.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"
#include "../environment/exception.hpp"
#include "radiator.hpp"

namespace flick {
  namespace filter {     
    class filter {
      //protected:
      //pp_function spectrum_;
    public:
      virtual double transmittance(double wavelength) const = 0;
      //{
      //return spectrum_.value(wavelength);
      //}
      /*
      pp_function spectrum(size_t n_points) const {
	return importance_sampled(spectrum_,n_points);
      }
      */
    };
    
    class gaussian : public filter {
      double wl0_;
      double sigma_; 
    public:
      gaussian(double center_wavelength, double fwhm) {
	wl0_ = center_wavelength;
	sigma_ = fwhm/(2*sqrt(2*log(2)));
	//double dwl = 4 * sigma;
	//std::vector<double> wl = range(wl0-dwl,wl0+dwl,n_points).logspace();
	//std::vector<double> tr(n_points);
	//using namespace constants;
	//for (size_t i=0; i < n_points; ++i) {
	//  tr[i] = 1/(sigma*sqrt(2*pi))*exp(-0.5*pow((wl[i]-wl0)/sigma,2));
	//}
	//spectrum_ = pp_function{wl,tr};

      }
      double transmittance(double wavelength) const {
	using namespace constants;
	return 1/(sigma_*sqrt(2*pi))*exp(-0.5*pow((wavelength-wl0_)/sigma_,2));
      }
    };
    
    class cut_ends : public filter {
      double l_;
      double u_;
    public:
      cut_ends(double lower_edge, double upper_edge)
	: l_{lower_edge},u_{upper_edge} {}
      double transmittance(double wavelength) const {
	if (wavelength < l_ or wavelength > u_)
	  return 0;	
	return 1;
      }
    };

    class erythema : public filter
    //  CIE-standard McKinlay–Diffey erythemal action
    //  spectrum. McKinlay, A. F. & Diffey, B. L. (1987). "A reference
    //  action spectrum for ultraviolet induced erythema in human
    //  skin". CIE Journal. 6 (1): 17–22.
    {
    public:
      double transmittance(double wavelength) const {
	if (wavelength < 298e-9)
	  return 1;
	else if (wavelength < 328e-9)
	  return pow(10, 0.094*(298-wavelength*1e9));
	else
	  return pow(10, 0.015*(139-wavelength*1e9));
      }
    };
    
    class tabulated : public filter {
      pl_function spectrum_;
    public:
      tabulated(const std::string& file_name) {
	std::ifstream ifs(file_name);
	ifs >> spectrum_;
      }
      double transmittance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
    };

    /*
    class stack : public filter {
      std::vector<std::shared_ptr<filter>> filters_;
    public:
      void add(const filter& f) {
      }
      double transmittance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
    };
    */
  }
  
  pl_function transmit(const pl_function& radiation_spectrum, const filter::filter &f) {
    const std::vector<double>& wl = radiation_spectrum.x();
    const std::vector<double>& rs = radiation_spectrum.y();
    pl_function new_rs;
    for (size_t i=0; i<rs.size(); ++i)
      new_rs.append({wl[i], rs[i]*f.transmittance(wl[i])});
    return new_rs;
  }

  double uv_index(const pl_function& radiation_spectrum) {
    return 40*transmit(radiation_spectrum, filter::erythema()).integral();
  }
  double uva_index(const pl_function& radiation_spectrum) {
    return 40*transmit(radiation_spectrum, filter::erythema()).integral(315e-9,400e-9);
  }
  double uvb_index(const pl_function& radiation_spectrum) {
    return 40*transmit(radiation_spectrum, filter::erythema()).integral(280e-9,315e-9);
  }

  /*
  pp_function transmit(const radiator::radiator& r, const filter::filter &f,
		       size_t n_points) {
    return transmit(r.spectrum(),f,n_points);
  }
  */
}

#endif

