#ifndef flick_command_iop
#define flick_command_iop

#include "../basic_command.hpp"
#include "../../material/water/pure_water.hpp"

namespace flick {
  namespace command {
    struct iop : public basic_command {
      iop():basic_command("iop"){};
      void run() {
	if (a(1)=="pure_water") {
	  double from_wl = std::stod(a(3));
	  double to_wl = std::stod(a(4));
	  double n_points = std::stod(a(5));	
	  double T = std::stod(a(6));	
	  double S = std::stod(a(7));	
	  auto wls = range(from_wl, to_wl, n_points).logspace();
	  material::pure_water pw;
	  pw.temperature(T);
	  pw.salinity(S);
	  double value = 0;
	  for (auto wl:wls) {
	    pw.set(wavelength(wl));
	    if (a(2)=="absorption_length")
	      value = 1/pw.absorption_coefficient();
	    else if (a(2)=="scattering_length")
	      value = 1/pw.scattering_coefficient();
	    else if (a(2)=="refractive_index")
	      value = pw.real_refractive_index();
	    else
	      error();
	    std::cout << std::setprecision(4) << wl << " " << value << '\n';
	  }
	} else
	  error();
      }
    };    
  }
}

#endif
