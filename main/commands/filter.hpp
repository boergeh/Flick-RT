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
	  std::cout << flick::filter::cone_L();
	  return;
	}
	if (a(1)=="cone_M") {
	  std::cout << flick::filter::cone_M();
	  return;
	}
	if (a(1)=="cone_S") {
	  std::cout << flick::filter::cone_S();
	  return;
	}
	if (a(1)=="tristimulus_x") {
	  std::cout << flick::filter::tristimulus_x();
	  return;
	}
	if (a(1)=="tristimulus_y") {
	  std::cout << flick::filter::tristimulus_y();
	  return;
	}
	if (a(1)=="tristimulus_z") {
	  std::cout << flick::filter::tristimulus_z();
	  return;
	}

	std::string fname = a(1);
	auto f = flick::read<flick::pl_function>(fname);
	if (a(2)=="tristimulus") {
	  std::cout << flick::tristimulus(f);
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
