#ifndef flick_accurt_input
#define flick_accurt_input

#include "../environment/input_output.hpp"
#include "../environment/configuration.hpp"
#include "../numeric/function.hpp"
#include "../numeric/std_operators.hpp"

namespace flick {
namespace accurt {
  struct basic_api : public configuration {
    basic_api() {
      add<double>("SOURCE_ZENITH_ANGLE",0);
    }
  };
  
  struct transmittance_api : public basic_api {
  };

  class basic_config : public configuration {
  public:
    basic_config() {
      add<std::string>("SOURCE_SCALING_FACTOR","constant_one");
      add<std::string>("BOTTOM_BOUNDARY_SURFACE","loamy_sand");      
    }
  };

  struct transmittance_config : public basic_config {
    transmittance_config() {
      add_configuration(read<transmittance_api>("./transmittance_api.txt"));
	  
    }    
  };

  using reflectance_config = transmittance_config;

  struct irradiance {
    std::vector<double> wavelengths;
    std::vector<double> depths;
    std::vector<std::vector<double>> values;
  };
  
  irradiance read_irradiance(const std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs)
      throw std::invalid_argument(file_name+" not found");
    size_t n_runs, n_streams, n_depths, n_wavelengths;
    ifs >> n_runs >> n_streams >> n_depths >> n_wavelengths;
    std::cout << n_wavelengths << std::endl;
    std::vector<double> depths(n_depths);
    std::vector<double> wavelengths(n_wavelengths);
    std::vector<std::vector<double>> values(n_depths,std::vector<double>(n_wavelengths));
    for (size_t i=0; i < n_depths; i++)
      ifs >> depths[i];
    for (size_t i=0; i < n_wavelengths; i++)
      ifs >> wavelengths[i];
    for (size_t i=0; i < n_depths; i++) {
      for (size_t j=0; j < n_wavelengths; j++) {
	ifs >> values[i][j];
      }
    }
    ifs.close();
    irradiance irr;
    irr.wavelengths = wavelengths;
    irr.depths = depths;
    irr.values = values;
    return irr;
  }

  pe_function transmittance() {
    std::string path = "./tmpOutput";
    std::string f_down = path+"/cosine_irradiance_total_downward.txt";
    irradiance irr_down = read_irradiance(f_down);
    size_t n = 1;
    if (irr_down.depths.at(n) > 100e3)
      n++;
    return pe_function{irr_down.wavelengths,irr_down.values.at(n+1)/irr_down.values.at(n)};
  }
  pe_function reflectance() {
    std::string path = "./tmpOutput";
    std::string f_down = path+"/cosine_irradiance_total_downward.txt";
    std::string f_up = path+"/cosine_irradiance_total_upwnward.txt";
    irradiance irr_down = read_irradiance(f_down);
    irradiance irr_up = read_irradiance(f_up);
    size_t n = 1;
    if (irr_down.depths.at(n) > 100e3)
      n++;
    return pe_function{irr_down.wavelengths,irr_up.values.at(n)/irr_up.values.at(n)};
  }

  pe_function run(configuration c) {
    c.set_text_qualifiers("#","#");
    //write(c,"./tmp.txt");
    //system("mkdir ./tmpMaterials");
    //system("mkdir ./tmpOutput");
    //write(material, "./tmpMaterials/"+material);

    //system("DYLD_LIBRARY_PATH=$ACCURT_PATH/lib AccuRT tmp");
    //std::string file_name = "./tmpOutput/cosine_irradiance_total_downward.txt";

    //std::cout << reflectance
    
    return transmittance();
  }
}
}

#endif
