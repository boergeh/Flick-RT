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
	  }
	  else if (a(2) == "ocean_radiance") {
	    configuration_template::ocean_radiance().write(a(3));
	  }
	  else if (a(2) == "boa_transmittance") {
	    configuration_template::boa_transmittance().write(a(3));
	  }
	  else if (a(2) == "rs_reflectance") {
	    configuration_template::rs_reflectance().write(a(3));
	  }
	  else if (a(2) == "ramses_reflectance") {
	    configuration_template::ramses_reflectance().write(a(3));
	  }
	}
	else if (size() == 2) {
	  auto c0 = read<configuration_template::basic_accurt>(a(1));
	  std::string name = c0.get<std::string>("material_name");
	  std::shared_ptr<material::base> m;
	  if (name == "atmosphere") {
	    auto c = read<configuration_template::atmosphere>(a(1));
	    m = std::make_shared<material::atmosphere>(c);
	  }
	  else if(name == "atmosphere_ocean") {
	    auto c = read<configuration_template::atmosphere_ocean>(a(1));
	    m = std::make_shared<material::atmosphere_ocean>(c);
	  }
	  else {
	    throw std::runtime_error("accurt command: material '"+name+"' not found");
	  }
	  flick::accurt a(c0,m);
	  if (a.radiance_distribution_override()) {
	    auto r = a.relative_distributed_radiance(); 
	    a.print_radiance_distribution(r);
	  }
	  else
	    std::cout << a.relative_radiation();
	}
	else {
	  error();
	}
      }
    };    
  }
}

#endif
