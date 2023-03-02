#ifndef flick_command_delta_fit
#define flick_command_delta_fit

#include "basic_command.hpp"
#include "../../numeric/legendre/delta_fit.hpp"

namespace flick {
  namespace command {
    struct delta_fit : public basic_command {
      delta_fit():basic_command("delta_fit"){};
      void run() {
	size_t n_terms = std::stoi(a(1));
	std::string fname = a(3);
	auto f = flick::read<flick::pl_function>("./"+fname);
	if (a(2)=="coefficients") {
	   std::cout << flick::delta_fit(f,n_terms).coefficients();
	}
	else if (a(2)=="function_value") {
	  std::cout << flick::delta_fit(f,n_terms).function_values({-1, -0.766,0.999,0.9999,0.99999758 ,1});
	}
	else
	  error();
      }
    };    
  }
}

#endif
