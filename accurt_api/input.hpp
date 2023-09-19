#ifndef flick_accurt_input
#define flick_accurt_input

#include "../environment/input_output.hpp"
#include "../environment/configuration.hpp"
#include "../numeric/function.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/table.hpp"
#include "../material/material.hpp"
#include "../material/layered_iops.hpp"

namespace flick {
  struct accurt_configuration : public configuration {
    accurt_configuration() {
      add<double>("DETECTOR_HEIGHT",0);
      add<std::string>("DETECTOR_ORIENTATION","up");
      add<double>("REFERENC_DETECTOR_HEIGHT",100e3);
      add<std::string>("REFERENCE_DETECTOR_ORIENTATION","up");
      add<double>("SOURCE_ZENITH_ANGLE",0);
      add<std::string>("BOTTOM_BOUNDARY_SURFACE","loamy_sand");
      add<double>("LAYER_DEPTHS_UPPER_SLAB",{30.0e3, 50.0e3, 60.0e3, 70.0e3, 
					     76.0e3, 80.0e3, 84.0e3, 88.0e3, 
					     90.0e3, 92.0e3, 94.0e3, 96.0e3, 
					     98.0e3, 100.0e3});
      add<double>("DETECTOR_WAVELENGTHS",{400e-9,500e-9});
      
    }
  };

  class accurt {
    configuration c_;
    std::shared_ptr<material::base> material_;
    const std::string path_ = "./tmpOutput";
    const size_t n_detector_ = 1;
    const size_t n_reference_ = 2;
    stdvector wavelengths_;
  public:
    accurt(const configuration& c, std::shared_ptr<material::base> material)
      : c_{c}, material_{material} {
      c_.add<std::string>("SOURCE_SCALING_FACTOR","constant_one");
      c_.set_text_qualifiers("#","#");
      wavelengths_ = c_.get_vector<double>("WAVELENGTHS");
    }
    /*
    void set_vertical_radiance() {
      c_.add<std::string>("SAVE_RADIANCE","true");
      c_.add<std::string>("DETECTOR_AZIMUTH_ANGLES","nan");
      c_.add<double>("DETECTOR_POLAR_ANGLES",{0, 180});  
    }
    */
    //void set_radiance(size_t n_polar, size_t n_azimuth) {
    //c_.add<std::string>("SAVE_RADIANCE","true");
      //c_.add<std::string>("DETECTOR_AZIMUTH_ANGLES","nan");
      //c_.add<double>("DETECTOR_POLAR_ANGLES",{0, 180});  
    // }
    void make_material_files() {
      stdvector d = c_.get_vector<double>("LAYER_DEPTHS_UPPER_SLAB");
      stdvector layer_bottoms = {0, 10e3};//100e3+(-1.0)*std::reverse(d.begin(),d.end());
      size_t n_terms = c_.get<size_t>("STREAM_UPPER_SLAB_SIZE") + 1;
      layered_iops layered_atmosphere(*material_,layer_bottoms,n_terms);
      write(accurt_user_specified(layered_atmosphere,wavelengths_),
	    "./tmpMaterial/user_specified_atmosphere");
      
      //bottoms = {-1, -10, -100};
      //layered_iops layered_ocean(m,bottoms,n_terms); // more terms needed
      //write(accurt_user_specified(layered_ocean,wavelengths_),
      //    "./tmpMaterial/user_specified_ocean");
    }
    pp_function relative_irradiance() {
      run();
      pe_table irr_down =
	read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      pe_table irr_up =
	read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      stdvector irr = irr_down.row(n_detector_).y();
      stdvector irr_ref = irr_down.row(n_reference_).y();   
      if (c_.get<std::string>("DETECTOR_ORIENTATION")=="up") {
	irr = irr_down.row(n_detector_).y();
      }
      if (c_.get<std::string>("REFERENCE_DETECTOR_ORIENTATION")=="up") {
	irr_ref = irr_down.row(n_reference_).y();
      }
      return pp_function{wavelengths_, irr/irr_ref};	  
    }
  private:
    void run() {
      write(c_,"./tmp");
      system("mkdir ./tmpMaterials");
      system("mkdir ./tmpOutput");
      make_material_files();
      system("DYLD_LIBRARY_PATH=$ACCURT_PATH/lib AccuRT tmp");
    }
    pe_table read_irradiance(const std::string& file_name) {
      std::ifstream ifs(file_name);
      if (!ifs)
	throw std::invalid_argument(file_name+" not found");
      size_t n_runs, n_streams;
      ifs >> n_runs >> n_streams;
      pe_table t;
      ifs >> t;
      std::cout << t;
      return t;
    }
  };
  
  /*
  pe_table read_irradiance(const std::string& file_name) {
    std::ifstream ifs(file_name);
    if (!ifs)
      throw std::invalid_argument(file_name+" not found");
    size_t n_runs, n_streams;
    ifs >> n_runs >> n_streams;
    pe_table t;
    ifs >> t;
    std::cout << t;
    return t;
  */
    /*
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
    */
  //}

  /*
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
  */
      /*
  void run(configuration c) {
    c.set_text_qualifiers("#","#");
    std::string path = "./tmpOutput";
    std::string f_down = path+"/cosine_irradiance_total_downward.txt";
    pe_table irr = read_irradiance(f_down);
      */
    //write(c,"./tmp.txt");
    //system("mkdir ./tmpMaterials");
    //system("mkdir ./tmpOutput");
    //write(material, "./tmpMaterials/"+material);

    //system("DYLD_LIBRARY_PATH=$ACCURT_PATH/lib AccuRT tmp");
    //std::string file_name = "./tmpOutput/cosine_irradiance_total_downward.txt";

    //std::cout << reflectance
    
    //return transmittance();
  //}
}

#endif
