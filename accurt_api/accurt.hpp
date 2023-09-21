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
      add<std::string>("SOURCE_TYPE","constant_one");
      add<double>("SOURCE_SCALING_FACTOR",1);
      add<double>("SOURCE_ZENITH_ANGLE",0);
      add<std::string>("BOTTOM_BOUNDARY_SURFACE","loamy_sand");
      add<double>("BOTTOM_BOUNDARY_SURFACE_SCALING_FACTOR",1);
      add<size_t>("STREAM_UPPER_SLAB_SIZE",16);
      add<double>("STREAM_LOWER_SLAB_PARAMETERS",{1,2});
      add<double>("LAYER_DEPTHS_UPPER_SLAB",{30.0e3, 50.0e3, 60.0e3, 70.0e3, 
					     76.0e3, 80.0e3, 84.0e3, 88.0e3, 
					     90.0e3, 92.0e3, 94.0e3, 96.0e3, 
					     98.0e3, 100.0e3});
      add<double>("LAYER_DEPTHS_LOWER_SLAB",{1});
      add<std::string>("MATERIALS_INCLUDED_UPPER_SLAB","user_specified_atmosphere");
      add<std::string>("MATERIALS_INCLUDED_LOWER_SLAB","vacuum");
      add<double>("DETECTOR_DEPTHS_UPPER_SLAB",{0,50e3});
      add<double>("DETECTOR_DEPTHS_LOWER_SLAB",{0.5});
      add<std::string>("DETECTOR_AZIMUTH_ANGLES","0:20:180");
      add<std::string>("DETECTOR_POLAR_ANGLES","0:2:180");
      add<double>("DETECTOR_WAVELENGTHS",{400,500});
      add<double>("DETECTOR_WAVELENGTH_BAND_WIDTHS",{270,0.01,4000,0.01});
      add<std::string>("SAVE_COSINE_IRRADIANCE","true");
      add<std::string>("SAVE_SINE_IRRADIANCE","false");
      add<std::string>("SAVE_SCALAR_IRRADIANCE","false");
      add<std::string>("SAVE_RADIANCE","false");
      add<std::string>("SAVE_IOPS","false");
      add<std::string>("SAVE_BOTTOM_BOUNDARY_SURFACE","false");
      add<std::string>("SAVE_MATERIAL_PROFILE","false");
      add<double>("PROFILE_OUTPUT_WAVELENGTH",500);
      add<std::string>("PRINT_PROGRESS_TO_SCREEN","true");
      add<size_t>("REPEATED_RUN_SIZE",1);
      add<std::string>("USE_POLARIZATION","false");
      add<std::string>("DO_OCEAN_PHASEMATRIX","false");
      add<double>("ACCURACY",0.0);
      add<double>("MIE_CALCULATOR",2);
      add<std::string>("DO_DELTA_FIT","false");
      add<double>("DELTA_FIT_TRUNCATION",3.0);
      add<double>("RESPONSE_FUNCTION_TYPE",1);
      add<std::string>("DO_SPHERICAL_CORRECTION","true");
      add<std::string>("DO_2D_ROUGH_SEA_SURFACE","false");
      add<double>("SURFACE_WIND_SPEED",6);
      add<double>("RELATIVE_WIND_DIRECTION",0);
      add<std::string>("USRANG","true");
      add<double>("LPICK",0);
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
      c_.set_text_qualifiers("#","#");
      wavelengths_ = c_.get_vector<double>("DETECTOR_WAVELENGTHS")*1e-9;
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
      stdvector layer_boundaries = {0, 10e3, 100e3};//100e3+(-1.0)*std::reverse(d.begin(),d.end());
      size_t n_terms = c_.get<size_t>("STREAM_UPPER_SLAB_SIZE") + 1;
      layered_iops layered_atmosphere(*material_,layer_boundaries,n_terms);
      write(accurt_user_specified(layered_atmosphere,wavelengths_),
      	    "./tmpMaterials/user_specified_atmosphere");
      
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
      system("mkdir -p ./tmpMaterials");
      system("mkdir -p ./tmpOutput");
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
