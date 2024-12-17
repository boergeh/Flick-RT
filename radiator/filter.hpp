#ifndef flick_filter
#define flick_filter

#include "../environment/input_output.hpp"
#include "../numeric/linalg/matrix.hpp"
#include "../numeric/function.hpp"
#include "../numeric/flist.hpp"
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

    class photons : public filter
    // Converts to number of photons per wavelength
    {
    public:
      double transmittance(double wavelength) const {
	using namespace constants;
	return wavelength / (h * c);
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

    template<int Lms_no>
    class cone_lms : public filter
    // Normalized human cone cell spectral sensitivity
    {
      std::string p = path()+"/radiator/filter_data/cie_photometry";
      pl_flist flist = read<pl_flist>(p+"/cie_lms_cf_2deg_1nm.txt");
      pl_function f = flist(Lms_no);
    public:
      cone_lms() {
	f.scale_x(1e-9);
      }
      const std::vector<double>& wavelength_grid() const {
	return f.x();
      }
      double transmittance(double wavelength) const {
	return f.value(wavelength);
      }
      friend std::ostream& operator<<(std::ostream &os, const cone_lms<Lms_no>& c) {
	os << c.f;
      return os;
      }
    };
    
    template<int Xyz_no>
    class xyz_bar : public filter
    // Color matching functions
    {
      std::string p = path()+"/radiator/filter_data/cie_photometry";
      pl_flist flist = read<pl_flist>(p+"/cie_xyz_1931_2deg.txt");
      pl_function f = flist(Xyz_no);
    public:
      xyz_bar() {
	f.scale_x(1e-9);
      }
      void use_cone_fundamentals() {
	cone_lms<0> L;
	cone_lms<1> M;
	cone_lms<2> S;
	std::vector<std::vector<double>> m =
	  {{1.94735469, -1.41445123, 0.36476327},
	   {0.68990272, 0.34832189, 0},
	   {0, 0, 1.93485343}};
	const std::vector<double>& x = L.wavelength_grid();
	f.clear();
	for (size_t i = 0; i < x.size(); i++) {
	  double y = L.transmittance(x[i]) * m[Xyz_no][0] +
	    M.transmittance(x[i]) * m[Xyz_no][1] +
	    S.transmittance(x[i]) * m[Xyz_no][2];
	  f.append(point(x[i]*1e-9,y));
	}
      }
      double transmittance(double wavelength) const {
	return f.value(wavelength);
      }
      friend std::ostream& operator<<(std::ostream &os, const xyz_bar<Xyz_no>& c) {
	os << c.f;
	return os;
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
	f.add_zero_extrapolation();
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

  template<typename Function>
  pl_function transmit(const Function& radiation_spectrum, const filter::filter &f,
		    const std::vector<double>& wavelengths={}) {
    std::vector<double> wl = wavelengths;
    if (wl.empty())
      wl = radiation_spectrum.x();
    pl_function new_rs;
    for (size_t i=0; i < wl.size(); ++i)
      new_rs.append({wl[i], radiation_spectrum.value(wl[i])*f.transmittance(wl[i])});
    return new_rs;
  }

  template<typename Function>
  double n_photons(const Function& radiation_spectrum, double wl_low, double wl_high) {
    return transmit(radiation_spectrum, filter::photons()).integral(wl_low, wl_high);
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
  template<typename Function>
  std::vector<double> chromaticity(const Function& radiation_spectrum) {
    double wl1 = 380e-9;
    double wl2 = 780e-9;
    const Function& s = radiation_spectrum;
    std::vector<double> xyz(3);
    xyz[0] = transmit(s, filter::xyz_bar<0>()).integral(wl1,wl2);
    xyz[1] = transmit(s, filter::xyz_bar<1>()).integral(wl1,wl2);
    xyz[2] = transmit(s, filter::xyz_bar<2>()).integral(wl1,wl2);
    double sum = 0;
    for (int i = 0; i < xyz.size(); i++)
      sum += xyz[i];
    //double sum = std::reduce(xyz.begin(), xyz.end()); 
    for (int i = 0; i < xyz.size(); i++) {
      xyz[i] = xyz[i]/sum;
    }
    return xyz;
  }
  template<typename Function>
  std::vector<double> rgb(const Function& radiation_spectrum) {
    using namespace linalg;
    linalg::matrix sRGB_D65 =  // White for CIE D65 spectrum 
      {{3.2404542,-1.5371385,-0.4985314},
       {-0.9692669,1.8760108,0.0415560},
       {0.0556434,-0.2040259,1.0572252}};
    linalg::matrix xyz = {chromaticity(radiation_spectrum)};
    std::vector<double> rgb = linalg::t(sRGB_D65*linalg::t(xyz))[0];
    double max = *std::max_element(rgb.begin(), rgb.end());
    double gamma = 1/2.2;
    for (int i = 0; i < 3; i++) {
      rgb[i] = rgb[i]/max;
      rgb[i] = pow(rgb[i],gamma);
    }
    return rgb;
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

