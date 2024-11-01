#ifndef flick_filter
#define flick_filter

#include "../environment/input_output.hpp"
#include "../numeric/function.hpp"
#include "../numeric/range.hpp"
#include "../numeric/distribution.hpp"
#include "../environment/exception.hpp"
#include "radiator.hpp"

namespace flick {
  namespace filter {     
    class filter {
    public:
      virtual double transmittance(double wavelength) const = 0;
    };
    
    class gaussian : public filter {
      double wl0_;
      double sigma_; 
    public:
      gaussian(double center_wavelength, double fwhm) {
	wl0_ = center_wavelength;
	sigma_ = fwhm/(2*sqrt(2*log(2)));
      }
      double transmittance(double wavelength) const {
	using namespace constants;
	return 1/(sigma_*sqrt(2*pi))*exp(-0.5*pow((wavelength-wl0_)/sigma_,2));
      }
    };

    class triangular : public filter {
      distribution::triangular t_;
    public:
      triangular(double center_wavelength, double fwhm)
	: t_(center_wavelength-fwhm,center_wavelength+fwhm,center_wavelength) {
 
      }
      double transmittance(double wavelength) const {
	return t_.pdf(wavelength);
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
    //  CIE/ISO-standard erythemal action spectrum. Schmalwieser,
    //  A.W., Wallisch, S. and Diffey, B., 2002. A library of action
    //  spectra for erythema and pigmentation. Photochemical &
    //  Photobiological Sciences, 1, pp.251-268.
    {
    public:
      double transmittance(double wavelength) const {
	if (wavelength < 298e-9)
	  return 1;
	else if (wavelength < 328e-9)
	  return pow(10, 0.094*(298-wavelength*1e9));
	else
	  return pow(10, 0.015*(140-wavelength*1e9));
      }
    };
    
    class tabulated : public filter {
      pl_function t_;
      double wavelength_shift_ = 0;
    public:
      tabulated(const pl_function& filter_transmittance)
	: t_{filter_transmittance} {
      }
      double transmittance(double wavelength) const {
	return t_.value(wavelength-wavelength_shift_);
      }
      void shift(double wavelength) {
	wavelength_shift_ = wavelength;
      }
    };

    class sentinel3 : public filter
    // Use closest normalized sentinel3 OLIC spectral response
    // function around a given user wavelength
    {
      pl_function centers_;
      std::shared_ptr<tabulated> srf_;
      double user_center_wavelength_;
    public:
      sentinel3(double user_center_wavelength)
	: user_center_wavelength_{user_center_wavelength} {
	std::string p = path()+"/radiator/filter_data/sentinel3/srf";
	centers_ = read<pl_function>(p+"/center_wavelength.txt");
	centers_.scale_x(1e-9);
	size_t n = closest_srf();
	pl_function f = read<pl_function>(p+"/band_"+std::to_string(n)+".txt");
	f.scale_x(1e-9);
	f.scale_y(1/f.integral());
	f.add_extrapolation_points(0);
	srf_ = std::make_shared<tabulated>(tabulated(f));
	srf_->shift(user_center_wavelength_-centers_.x()[n]);
      }
      double transmittance(double wavelength) const {
	return srf_->transmittance(wavelength);
      }
      size_t closest_srf() {
	double n = centers_.value(user_center_wavelength_);
	if (n > centers_.size()-1)
	  return centers_.size()-1;
	if (n < 0)
	  return 0;
	return std::round(n);
      }
    };
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
  double gaussian_mean(const pl_function& radiation_spectrum, double wl0, double fwhm) {
    return transmit(radiation_spectrum, filter::gaussian(wl0,fwhm)).integral();
  }
  double triangular(const pl_function& radiation_spectrum, double wl0, double fwhm) {
    return transmit(radiation_spectrum, filter::triangular(wl0,fwhm)).integral();
  }
  double square_mean(const pl_function& radiation_spectrum, double wl0, double full_width) {
    return transmit(radiation_spectrum, filter::cut_ends(wl0-full_width/2,wl0+full_width/2)).integral()/full_width;
  }
  double weighted_integral(const pl_function& radiation_spectrum,
			  const pl_function& filter_transmittance) {
    return transmit(radiation_spectrum, filter::tabulated(filter_transmittance)).integral();
  }
  double sentinel3(const pl_function& radiation_spectrum, double wl0) {
    return transmit(radiation_spectrum, filter::sentinel3(wl0)).integral();
  }  
}

#endif

