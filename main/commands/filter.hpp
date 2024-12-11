#ifndef flick_command_filter
#define flick_command_filter

#include "basic_command.hpp"
#include "../../radiator/filter.hpp"

namespace flick {
  namespace command {
    struct filter : public basic_command {
      filter():basic_command("filter"){};
      void run() {
	if (a(1)=="cone_L") {
	  std::cout << flick::filter::cone_lms<0>();
	  return;
	}
	if (a(1)=="cone_M") {
	  std::cout << flick::filter::cone_lms<1>();
	  return;
	}
	if (a(1)=="cone_S") {
	  std::cout << flick::filter::cone_lms<2>();
	  return;
	}
	if (a(1)=="x_bar") {
	  std::cout << flick::filter::xyz_bar<0>();
	  return;
	}
	if (a(1)=="y_bar") {
	  std::cout << flick::filter::xyz_bar<1>();
	  return;
	}
	if (a(1)=="z_bar") {
	  std::cout << flick::filter::xyz_bar<2>();
	  return;
	}

	std::string fname = a(1);
	auto f = flick::read<flick::pl_function>(fname);
	if (a(2)=="chromaticity") {
	  std::cout << flick::chromaticity(f);
	}
	else if (a(2)=="rgb") {
	  std::cout << flick::rgb(f);
	}
	else if (a(2)=="uv_index") {
	  std::cout << flick::uv_index(f);
	}
	else if (a(2)=="uva_index") {
	  std::cout << flick::uva_index(f);
	}
	else if (a(2)=="uvb_index") {
	  std::cout << flick::uvb_index(f);
	}
	else if (a(2)=="n_photons") {
	  std::cout << flick::n_photons(f,std::stod(a(3)),std::stod(a(4)));
	}
	else if (a(2)=="curvature_sampled") {
	  int n = std::stoi(a(3));
	  pl_function pdf = absolute(derivative(derivative(f))).normalize();
	  std::vector<double> fx = pdf.x();
	  std::vector<double> fy = pdf.y();
	  double delta_x = fx.back()-fx.front();
	  for (size_t i=0; i<fx.size(); i++) {
	    fy[i] = 1-exp(-fy[i]*delta_x);
	  }
	  pdf = pl_function{fx,fy};
	  pdf.normalize();
	  pl_function cdf = accumulate(pdf);
	  pl_function g = remove_non_increasing_values(cdf);
	  pl_function quantile = invert(g);
	  std::vector<double> x(n);
	  std::vector<double> y(n);
	  std::vector<double> p = range(0,1,n).linspace();
	  for (size_t i=0; i<y.size(); i++) {
	    double wl = quantile.value(p[i]);
	    x[i] = wl;
	    y[i] = f.value(wl);
	  }
	  std::cout << pl_function{x,y};	  
	}
	else if (a(2)=="gaussian_mean") {
	  double wl0 = std::stod(a(3));
	  double fwhm = std::stod(a(4));
	  std::cout << flick::gaussian_mean(f,wl0,fwhm);
	}
	else if (a(2)=="triangular") {
	  double wl0 = std::stod(a(3));
	  double fwhm = std::stod(a(4));
	  std::cout << flick::triangular(f,wl0,fwhm);
	}
	else if (a(2)=="weighted_integral") {
	  std::string fname = a(3);
	  auto f2 = flick::read<flick::pl_function>(fname);
	  std::cout << flick::weighted_integral(f,f2);
	}
	else if (a(2)=="sentinel3") {
	  double wl0 = std::stod(a(3));
	  std::cout << flick::sentinel3(f,wl0);
	}
	else
	  error();
      }
    };    
  }
}

#endif
