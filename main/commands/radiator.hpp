#include "../basic_command.hpp"
#include "../../radiator/planck.hpp"

namespace flick {
  namespace command {
    struct radiator : public basic_command {
      radiator():basic_command("radiator"){};
      void run() {
	if (a(1)=="planck") {
	  double T = std::stod(a(2));
	  size_t n_points = std::stoi(a(3));
	  std::cout << planck(T).irradiance_spectrum(n_points);
	} else {
	  error();
	}
      }
    };
  }
}
