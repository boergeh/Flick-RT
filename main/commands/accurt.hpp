#ifndef flick_command_accurt
#define flick_command_accurt

#include "basic_command.hpp"
#include "../../material/atmosphere_ocean.hpp"
#include "../../environment/input_output.hpp"
#include "../../environment/configuration.hpp"
#include "../../accurt_api/accurt.hpp"
#include <filesystem>

namespace flick {
  namespace command {
    class accurt : public basic_command {
    public:
      accurt() : basic_command("accurt") {};
      void run() {
	if (size() == 3 && a(1) == "toa_reflectance") {
	  if (not std::filesystem::exists(a(2))) {
	    configuration_template::toa_reflectance().write(a(2));
	  }
	  auto c = read<configuration_template::toa_reflectance>("./"+a(2));
	  auto m = std::make_shared<material::atmosphere_ocean>(c);
	  std::cout << flick::accurt(c,m).relative_radiation();
	}
	else if (size() == 3 && a(1) == "boa_transmittance") {
	  if (not std::filesystem::exists(a(2))) {
	    configuration_template::boa_transmittance().write(a(2));
	  }
	  auto c = read<configuration_template::toa_reflectance>("./"+a(2));
	  auto m = std::make_shared<material::atmosphere_ocean>(c);
	  std::cout << flick::accurt(c,m).relative_radiation();
	}
      }
    };    
  }
}

#endif
