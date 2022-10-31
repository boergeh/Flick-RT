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
    protected:
      pp_function spectrum_;
    public:
      double transmittance(double wavelength) const {
	return spectrum_.value(wavelength);
      }
      pp_function spectrum(size_t n_points) const {
	return importance_sampled(spectrum_,n_points);
      }
    };
    
    class gaussian : public filter {
    public:
      gaussian(double center_wavelength, double fwhm, size_t n_points) {
	double wl0 = center_wavelength;
	double sigma = fwhm/(2*sqrt(2*log(2)));
	double dwl = 4 * sigma;
	std::vector<double> wl = range(wl0-dwl,wl0+dwl,n_points).logspace();
	std::vector<double> tr(n_points);
	using namespace constants;
	for (size_t i=0; i < n_points; ++i) {
	  tr[i] = 1/(sigma*sqrt(2*pi))*exp(-0.5*pow((wl[i]-wl0)/sigma,2));
	}
	spectrum_ = pp_function{wl,tr};
      }
    };
    
    class tabulated : public filter {
    public:
      tabulated(const std::string& file_name) {
	std::ifstream ifs(file_name);
	ifs >> spectrum_;
      }
    };
    
    class stack : public filter {
    public:
      void add(const filter& f, size_t n_points) {
	if (spectrum_.size() == 0)
	  spectrum_ = f.spectrum(n_points);
	else {
	  pp_function s = f.spectrum(n_points);
	  spectrum_ = multiply(spectrum_, s, s.x());
	}
      }      
    };
  }
  
  pp_function transmit(const radiator::radiator& r, const filter::filter &f,
		       size_t n_points) {
    pp_function fs = f.spectrum(n_points);
    return multiply(r.spectrum(),fs,fs.x());
  }
}

#endif

