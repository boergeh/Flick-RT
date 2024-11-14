#include "basic_command.hpp"
#include "../../radiator/planck.hpp"
#include "../../radiator/toa_solar.hpp"
#include "../../radiator/cie_d65.hpp"
#include "../../radiator/cie_a.hpp"

namespace flick {
  namespace command {
    struct radiator : public basic_command {
      radiator():basic_command("radiator"){};
      void run() {
	if (a(1)=="planck") {
	  double T = std::stod(a(2));
	  size_t n_points = std::stoi(a(3));
	  std::cout << flick::radiator::planck(T).importance_spectrum(n_points);
	}
	else if (a(1)=="toa_solar") {
	  auto s = flick::radiator::toa_solar().spectrum();
	  std::cout << s.header("");
	}
	else if (a(1)=="cie_d65") {
	  auto s = flick::radiator::cie_d65().spectrum();
	  std::cout << s.header("");
	}
	else if (a(1)=="cie_a") {
	  auto s = flick::radiator::cie_a().spectrum();
	  std::cout << s.header("");
	}
	else {
	  error();
	}
      }
    };
  }
}
