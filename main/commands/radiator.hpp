#include "basic_command.hpp"
#include "../../radiator/planck.hpp"
#include "../../radiator/toa_solar.hpp"

namespace flick {
  namespace command {
    struct radiator : public basic_command {
      radiator():basic_command("radiator"){};
      void run() {
	if (a(1)=="planck") {
	  double T = std::stod(a(2));
	  size_t n_points = std::stoi(a(3));
	  std::cout << flick::radiator::planck(T).spectrum(n_points);
	}
	else if (a(1)=="toa-solar") {
	  std::cout << flick::radiator::toa_solar().spectrum();
	}
	else {
	  error();
	}
      }
    };
  }
}
