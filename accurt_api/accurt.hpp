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
    class accurt_user_specified {
    layered_iops& iops_;
    stdvector wls_;
  public:
    accurt_user_specified(layered_iops& iops, const stdvector& wavelengths)
      : iops_{iops}, wls_{wavelengths} {}
  private:
    friend std::ostream& operator<<(std::ostream &os, const accurt_user_specified& u) {
      os << "# AccuRT configuration file for the user_specified material #\n"
	 << "PROFILE_LABEL = layer_numbering #\n"
	 << "MATERIAL_PROFILE = 1 #\n"
	 << "TURN_OFF_DELTA_FIT = true #\n\n"
	 << "WAVELENGTHS = ";
      for (auto& wl:u.wls_)
	os << wl*1e9 << " ";
      os << "#\n\n";
      std::vector<stdvector> n;
      std::vector<stdvector> a;
      std::vector<stdvector> s;
      std::vector<std::vector<stdvector>> p;
      for (auto& wl:u.wls_) {
	u.iops_.set_wavelength(wl);
	for (size_t j=0; j<u.iops_.n_layers(); j++) {
	  n.push_back(u.iops_.refractive_index());
	  a.push_back(u.iops_.absorption_coefficient());
	  s.push_back(u.iops_.scattering_coefficient());
	  p.push_back(u.iops_.alpha_terms(0));
	}
      }
      os << "REFRACTIVE_INDEX = ";
      for (size_t i=0; i<u.wls_.size(); i++) {
	os << *std::max_element(std::begin(n[i]), std::end(n[i])) << " ";
      }
      os << "#\n\n";
      size_t l = u.iops_.n_layers();
      for (int i=l-1; i>=0; i--) {
	for (size_t j=0; j<u.wls_.size(); j++) {
	  os << "A_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << a[i][j] << " #\n";
	  os << "S_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << s[i][j] << " #\n";
	  os << "P_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = ";
	  for (size_t k=0; k<p[i][j].size(); k++) {
	    os << p[i][j][k]/p[i][j][0]/(2*k+1) << " ";
	  }
	  os << " #\n";
	  os << "\n";
	}
      }  
      return os;
    }
  };
  
  class accurt {
  public:
    struct configuration : public flick::configuration {
      configuration() {
	add<double>("DETECTOR_HEIGHT",1,"[m]");
	add<std::string>("DETECTOR_ORIENTATION","up","Vertically <up> or <down>");
	add<std::string>("DETECTOR_TYPE","plane_irradiance","<plane_irradiance>, <scalar_irradiance>, or <radiance>");
	add<double>("DETECTOR_WAVELENGTHS",{400e-9,500e-9},"[m]");
	add<double>("REFERENCE_DETECTOR_HEIGHT",100e3,"[m] Calculated detector signal is divided by the calculated reference detector plane irradiance signal at a give height.");
	add<std::string>("REFERENCE_DETECTOR_ORIENTATION","up","<up> or <down>");
	add<double>("SOURCE_ZENITH_ANGLE",0,"[degrees] Zero gives overhead source.");
      }
    };
  private:
    configuration c_;
    std::shared_ptr<material::base> material_;
    const std::string path_ = "./tmpOutput";
    size_t n_detector_;
    size_t n_reference_;
    const double toa_ = 100e3;
    const double bottom_depth_ = 200;
    stdvector wavelengths_;
  public:
    accurt(const configuration& c, std::shared_ptr<material::base> material)
      : c_{c}, material_{material} {
      c_.add<std::string>("SOURCE_TYPE","constant_one"); 
      c_.add<double>("SOURCE_SCALING_FACTOR",1);
      c_.add<std::string>("BOTTOM_BOUNDARY_SURFACE","loamy_sand");
      c_.add<double>("BOTTOM_BOUNDARY_SURFACE_SCALING_FACTOR",1);
      c_.add<size_t>("STREAM_UPPER_SLAB_SIZE",10);
      c_.add<double>("STREAM_LOWER_SLAB_PARAMETERS",{1,2});
      c_.add<double>("LAYER_DEPTHS_UPPER_SLAB",{30.0e3, 50.0e3, 60.0e3, 70.0e3, 
					     76.0e3, 80.0e3, 84.0e3, 88.0e3, 
					     90.0e3, 92.0e3, 94.0e3, 96.0e3, 
					     98.0e3, 100.0e3});
      c_.add<double>("LAYER_DEPTHS_LOWER_SLAB",{bottom_depth_});
      c_.add<std::string>("MATERIALS_INCLUDED_UPPER_SLAB","user_specified_atmosphere");
      c_.add<std::string>("MATERIALS_INCLUDED_LOWER_SLAB","vacuum");
      add_detector_depths();
      c_.add<std::string>("DETECTOR_AZIMUTH_ANGLES","0:20:180");
      c_.add<std::string>("DETECTOR_POLAR_ANGLES","0:2:180");   
      c_.add<double>("DETECTOR_WAVELENGTH_BAND_WIDTHS",{270,0.01,4000,0.01});
      c_.add<std::string>("SAVE_COSINE_IRRADIANCE","true");
      c_.add<std::string>("SAVE_SINE_IRRADIANCE","false");
      c_.add<std::string>("SAVE_SCALAR_IRRADIANCE","false");
      c_.add<std::string>("SAVE_RADIANCE","false");
      c_.add<std::string>("SAVE_IOPS","false");
      c_.add<std::string>("SAVE_BOTTOM_BOUNDARY_SURFACE","false");
      c_.add<std::string>("SAVE_MATERIAL_PROFILE","false");
      c_.add<double>("PROFILE_OUTPUT_WAVELENGTH",500);
      c_.add<std::string>("PRINT_PROGRESS_TO_SCREEN","true");
      c_.add<size_t>("REPEATED_RUN_SIZE",1);
      c_.add<std::string>("USE_POLARIZATION","false");
      c_.add<std::string>("DO_OCEAN_PHASEMATRIX","false");
      c_.add<double>("ACCURACY",0.0);
      c_.add<double>("MIE_CALCULATOR",2);
      c_.add<std::string>("DO_DELTA_FIT","false");
      c_.add<double>("DELTA_FIT_TRUNCATION",3.0);
      c_.add<double>("RESPONSE_FUNCTION_TYPE",1);
      c_.add<std::string>("DO_SPHERICAL_CORRECTION","true");
      c_.add<std::string>("DO_2D_ROUGH_SEA_SURFACE","false");
      c_.add<double>("SURFACE_WIND_SPEED",6);
      c_.add<double>("RELATIVE_WIND_DIRECTION",0);
      c_.add<std::string>("USRANG","true");
      c_.add<double>("LPICK",0);
     
      wavelengths_ = c_.get_vector<double>("DETECTOR_WAVELENGTHS");
      c_.set<double>("DETECTOR_WAVELENGTHS",wavelengths_*1e9);
      c_.set_text_qualifiers("#","#");
    }
    pp_function relative_radiation() {
      if (c_.get<std::string>("DETECTOR_TYPE")=="radiance")
	return relative_radiance();
      else if (c_.get<std::string>("DETECTOR_TYPE")=="scalar_irradiance")
	return relative_scalar_irradiance();
      else
	return relative_plane_irradiance();
    }
  private:
    void set_vertical_radiance() {
      c_.set<std::string>("SAVE_RADIANCE","true");
      c_.set<std::string>("DETECTOR_AZIMUTH_ANGLES","nan");
      c_.set<std::string>("DETECTOR_POLAR_ANGLES","0 180");
      c_.set<size_t>("STREAM_UPPER_SLAB_SIZE",30);
    }
    void make_material_files() {
      stdvector d = c_.get_vector<double>("LAYER_DEPTHS_UPPER_SLAB");
      std::reverse(d.begin(),d.end());
      stdvector layer_boundaries = toa_-d;
      layer_boundaries.push_back(toa_);
      
      size_t n_terms = c_.get<size_t>("STREAM_UPPER_SLAB_SIZE") + 1;
      layered_iops layered_atmosphere(*material_,layer_boundaries,n_terms);
      write(accurt_user_specified(layered_atmosphere,wavelengths_),
      	    "./tmpMaterials/user_specified_atmosphere");
           
      //layered_iops layered_ocean(m,bottoms,n_terms); // more terms needed
      //write(accurt_user_specified(layered_ocean,wavelengths_),
      //    "./tmpMaterial/user_specified_ocean");
    }
    pp_function relative_plane_irradiance() {
      run();
      return pp_function{wavelengths_,
	detector_plane_irradiance()/reference_detector_irradiance()};	  
    }
    pp_function relative_scalar_irradiance() {
      run();
      return pp_function{wavelengths_,
	detector_scalar_irradiance()/reference_detector_irradiance()};	  
    }
    pp_function relative_radiance() {
      set_vertical_radiance();
      run();
      return pp_function{wavelengths_,
	detector_radiance()/reference_detector_irradiance()};	  
    }
    stdvector detector_plane_irradiance() {
      pe_table t;
      if (c_.get<std::string>("DETECTOR_ORIENTATION")=="up") {
	t = read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      } else {
	t = read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      }
      return t.row(n_detector_).y();
    }
    stdvector detector_scalar_irradiance() {
      c_.set<std::string>("SAVE_SCALAR_IRRADIANCE","true");      
      pe_table t;
      if (c_.get<std::string>("DETECTOR_ORIENTATION")=="up") {
	t = read_irradiance(path_+"/scalar_irradiance_total_downward.txt");
      } else {
	t = read_irradiance(path_+"/scalar_irradiance_total_upward.txt");
      }
      return t.row(n_detector_).y();
    }
    stdvector reference_detector_irradiance() {
      pe_table t;
      if (c_.get<std::string>("REFERENCE_DETECTOR_ORIENTATION")=="up") {
	t = read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      } else {
	t = read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      }
      return t.row(n_reference_).y();
    }
    stdvector detector_radiance() {
      //pe_table t = read_radiance(path_+"/radiance.txt");
      if (c_.get<std::string>("REFERENCE_DETECTOR_ORIENTATION")=="up") {
      } else {
      }
      return stdvector(wavelengths_.size());
    }    
    void add_detector_depths() {
      double h1 = c_.get<double>("DETECTOR_HEIGHT");    
      double h2 = c_.get<double>("REFERENCE_DETECTOR_HEIGHT");
      if (h1 > h2) {
	n_detector_ = 0;
	n_reference_ = 1;
      } else {
	n_detector_ = 1;
	n_reference_ = 0;
      }
      stdvector h = {h1,h2};
      std::sort(h.begin(),h.end());
      stdvector h_upper, h_lower;
      if (h[0] < 0) {
	h_upper = {0};
	h_lower = h;
	n_detector_++;
	n_reference_++;
      }
      else if (h[1] > 0) {
	h_upper = h;
	h_lower = {-bottom_depth_/2};
      }
      else if (h[0] > 0 and h[1] < 0) {
	h_upper = {h[0]};
	h_lower = {h[1]};
      }
      std::reverse(h_upper.begin(), h_upper.end());
      c_.add<double>("DETECTOR_DEPTHS_UPPER_SLAB",toa_-h_upper);
      c_.add<double>("DETECTOR_DEPTHS_LOWER_SLAB",0-h_lower);
    }
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
}

#endif
