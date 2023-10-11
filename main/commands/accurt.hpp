#ifndef flick_command_accurt
#define flick_command_accurt

#include "basic_command.hpp"
#include "../../material/atmosphere_ocean.hpp"
#include "../../environment/input_output.hpp"
#include "../../environment/configuration.hpp"
#include "../../accurt_api/accurt.hpp"

namespace flick {
  namespace command {
    class accurt : public basic_command {
    public:
      accurt() : basic_command("accurt") {};
      void run() {
	if (size() == 4 && a(1) == "-g" ) {
	  if (a(2) == "toa_reflectance") {
	    configuration_template::toa_reflectance().write(a(3));
	  } else if (a(2) == "boa_transmittance") {
	    configuration_template::boa_transmittance().write(a(3));
	  } else if (a(2) == "rs_reflectance") {
	    configuration_template::rs_reflectance().write(a(3));
	  } else if (a(2) == "ramses_reflectance") {
	    configuration_template::ramses_reflectance().write(a(3));
	  }
	} else if (size() == 2) {
	  auto c = read<configuration_template::basic_accurt>("./"+a(1));
	  auto m = std::make_shared<material::atmosphere_ocean>(c);
	  std::cout << flick::accurt(c,m).relative_radiation();
	}
	else {
	  error();
	}
      }
    };    
  }
}

#endif
