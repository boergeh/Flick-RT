#ifndef flick_accurt_input
#define flick_accurt_input

#include "../environment/input_output.hpp"
#include "../environment/configuration.hpp"
#include "../numeric/function.hpp"
#include "../numeric/std_operators.hpp"
#include "../numeric/table.hpp"
#include "../numeric/grid.hpp"
#include "../material/material.hpp"
#include "../material/water/pure_water.hpp"
#include "../material/layered_iops.hpp"
#include "../material/z_profile.hpp"

namespace flick {
  class accurt_user_specified {
    std::shared_ptr<layered_iops> iops_;
    stdvector wls_;
  public:
    accurt_user_specified(std::shared_ptr<layered_iops> iops, const stdvector& wavelengths)
      : iops_{iops}, wls_{wavelengths} {}
  private:
    friend std::ostream& operator<<(std::ostream &os, const accurt_user_specified& u) {
	os << "# AccuRT configuration file for the user_specified material #\n"
	 << "PROFILE_LABEL = layer_numbering #\n"
	 << "MATERIAL_PROFILE = ";
      for (size_t i = 0; i < u.iops_->n_layers(); i++) {
	os << i+1 << " 1  ";
      }
      os << " #\n"
	 << "TURN_OFF_DELTA_FIT = true #\n\n"
	 << "WAVELENGTHS = ";
      for (auto& wl:u.wls_)
	os << wl*1e9 << " ";
      os << "#\n\n";
      std::vector<stdvector> n(u.wls_.size());
      std::vector<stdvector> a(u.wls_.size());
      std::vector<stdvector> s(u.wls_.size());
      std::vector<std::vector<stdvector>> p(u.wls_.size());
      for (size_t i=0; i<u.wls_.size(); i++) {
	u.iops_->set_wavelength(u.wls_[i]);
	n[i] = u.iops_->refractive_index();
	a[i] = u.iops_->absorption_coefficient();
	s[i] = u.iops_->scattering_coefficient();
	p[i] = u.iops_->alpha_terms(0);
      }
      os << "REFRACTIVE_INDICES = ";
      for (size_t i=0; i<u.wls_.size(); i++) {
	os << *std::max_element(std::begin(n[i]), std::end(n[i])) << " ";
      }
      os << "#\n\n";
      size_t l = u.iops_->n_layers();
      for (int i=l-1; i>=0; i--) {
	for (size_t j=0; j<u.wls_.size(); j++) {
	  os << "A_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << a[j][i] << " #\n";
	  os << "S_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << s[j][i] << " #\n";
	  os << "P_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = ";
	  for (size_t k=0; k<p[j][i].size(); k++) {
	    os << p[j][i][k]/p[j][i][0]/(2*k+1) << " ";
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
    struct configuration : public basic_configuration {
      configuration() {
	add<std::string>("material_name","atmosphere_ocean","Name of material to be used with AccuRT.");
	add<std::string>("subtract_specular_radiance","false","Possible subtraction the nadir radiance specularly reflected on the water surface. Should be set 'true' when calculating remote sensing reflectance.");
	add<double>("detector_height",100e3,"Radiometer position relative to sea level [m]");
	add<std::string>("detector_orientation","up","Vertical orientation, <up> or <down>");
	add<std::string>("detector_type","irradiance","<plane_irradiance>, <scalar_irradiance>, or <radiance>");
	add<double>("DETECTOR_WAVELENGTHS",{400e-9,500e-9},"Wavelengths to be detected [m]");
	add<double>("reference_detector_height",100e3,"Calculated detector signal is divided by the calculated \nplane irradiance reference signal at a give height [m].");
	add<std::string>("reference_detector_orientation","up","<up> or <down>");
	add<double>("SOURCE_ZENITH_ANGLE",0,"Zero gives overhead source [degrees]");
	add<double>("BOTTOM_BOUNDARY_SURFACE_SCALING_FACTOR",1,"Zero gives black surface");
      }
    };
  private:
    basic_configuration c_;
    std::shared_ptr<material::base> material_;
    const std::string path_ = "./tmpOutput";
    size_t n_detector_;
    size_t n_reference_;
    double max_height_ = 1;
    double bottom_depth_ = 1;
    stdvector wavelengths_;
    const size_t precision_ = 12;
    const double minimum_height_ = 1e-6;
    const double minimum_depth_ = 1e-6;

  public:
    accurt(const basic_configuration& c, std::shared_ptr<material::base> material)
      : c_{c}, material_{material} {
      c_.add<std::string>("SOURCE_TYPE","constant_one"); 
      c_.add<double>("SOURCE_SCALING_FACTOR",1);
      c_.add<std::string>("BOTTOM_BOUNDARY_SURFACE","loamy_sand");
      c_.add<size_t>("STREAM_UPPER_SLAB_SIZE",10);
      c_.add<double>("STREAM_LOWER_SLAB_PARAMETERS",{1,2});
      add_layer_depths();
      c_.add<std::string>("MATERIALS_INCLUDED_UPPER_SLAB","user_specified_upper_slab");
      c_.add<std::string>("MATERIALS_INCLUDED_LOWER_SLAB","user_specified_lower_slab");
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
      c_.add<std::string>("PRINT_PROGRESS_TO_SCREEN","false");
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
      c_.add<std::string>("USRANG","false");
      c_.add<double>("LPICK",0);
     
      wavelengths_ = c_.get_vector<double>("DETECTOR_WAVELENGTHS");
      c_.set<double>("DETECTOR_WAVELENGTHS",wavelengths_*1e9);
      c_.set_text_qualifiers("#","##");
    }
    pp_function relative_radiation() {
      if (c_.get<std::string>("detector_type")=="radiance")
	return relative_radiance();
      else if (c_.get<std::string>("detector_type")=="scalar_irradiance")
	return relative_scalar_irradiance();
      else
	return relative_plane_irradiance();
    }
  private:
    void set_vertical_radiance() {
      c_.set<std::string>("SAVE_RADIANCE","true");
      c_.set<std::string>("DETECTOR_AZIMUTH_ANGLES","nan");
      c_.set<std::string>("DETECTOR_POLAR_ANGLES","0 180");
      c_.set<size_t>("STREAM_UPPER_SLAB_SIZE",20);
    }
    void make_material_files() {        
      size_t n_terms = c_.get<size_t>("STREAM_UPPER_SLAB_SIZE") + 1;

      stdvector b = depths_to_boundaries(c_.get_vector<double>("LAYER_DEPTHS_UPPER_SLAB"));
      auto layered_upper_slab = std::make_shared<layered_iops>(material_,b,n_terms);
      write(accurt_user_specified(layered_upper_slab, wavelengths_),
       	    "./tmpMaterials/user_specified_upper_slab", precision_);

      b = 0 - c_.get_vector<double>("LAYER_DEPTHS_LOWER_SLAB");
      std::reverse(b.begin(),b.end());
      b.push_back(0);
      auto layered_lower_slab = std::make_shared<layered_iops>(material_,b,n_terms);
      write(accurt_user_specified(layered_lower_slab, wavelengths_),
      	    "./tmpMaterials/user_specified_lower_slab", precision_);
    }
    stdvector depths_to_boundaries(stdvector depths) {
      std::reverse(depths.begin(),depths.end());
      stdvector layer_boundaries = max_height_-depths;
      layer_boundaries.push_back(max_height_);
      return layer_boundaries;
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
    stdvector detector_radiance() {
      grid_4d g = read_radiance(path_+"/radiance.txt");
      stdvector wls = g.x[1]*1e-9;
      stdvector Lu(wls.size());
      stdvector Ld(wls.size());
      for (size_t i=0; i<wls.size(); i++) {
	Lu[i] = g.f[0][i][1][0];
	Ld[i] = g.f[0][i][0][0];
      }
      if (c_.get<std::string>("subtract_specular_radiance")=="true") {
	Lu - Ld * nadir_fresnel_coefficient(wls);
      }
      if (c_.get<std::string>("detector_orientation")=="up") {
	return Ld;
      } else {
	return Lu;
      }
    }
    stdvector detector_plane_irradiance() {
      pe_table t;
      if (c_.get<std::string>("detector_orientation")=="up") {
	t = read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      } else {
	t = read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      }
      return t.row(n_detector_).y();
    }
    stdvector detector_scalar_irradiance() {
      c_.set<std::string>("SAVE_SCALAR_IRRADIANCE","true");      
      pe_table t;
      if (c_.get<std::string>("detector_orientation")=="up") {
	t = read_irradiance(path_+"/scalar_irradiance_total_downward.txt");
      } else {
	t = read_irradiance(path_+"/scalar_irradiance_total_upward.txt");
      }
      return t.row(n_detector_).y();
    }
    stdvector reference_detector_irradiance() {
      pe_table t;
      if (c_.get<std::string>("reference_detector_orientation")=="up") {
	t = read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      } else {
	t = read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      }
      return t.row(n_reference_).y();
    }
    stdvector nadir_fresnel_coefficient(const stdvector& wls) {
      material::pure_water pw;
      stdvector c(wls.size());
      for (size_t i=0; i<wls.size(); i++) {
	double n1 = 1;
	pw.set_wavelength(wls[i]);
	double n2 = pw.real_refractive_index(); 
	c[i] = pow((n1 - n2) / (n1 + n2),2); 
      }
      return c;
    }
    void add_layer_depths() {
      stdvector h = {-bottom_depth_, 0,max_height_};
      material::z_profile* zp = dynamic_cast<material::z_profile*>(&*material_);
      if (zp != NULL) {
	h = zp->height_grid();
      }      
      max_height_ = h.back();
      if (h.front()>=0)
      	h.insert(h.begin(),-bottom_depth_);
      bottom_depth_ = -h.front();
      stdvector depths = max_height_ - h;
      depths.pop_back();
      std::reverse(depths.begin(),depths.end());
      stdvector depths_upper, depths_lower;
      for (auto d:depths) {
	if (d > max_height_)
	  depths_lower.push_back(d-max_height_);
	else
	  depths_upper.push_back(d);
      }
      c_.add<double>("LAYER_DEPTHS_UPPER_SLAB",depths_upper);
      c_.add<double>("LAYER_DEPTHS_LOWER_SLAB",depths_lower);
    }
    void add_detector_depths() {
      double h1 = c_.get<double>("detector_height");    
      double h2 = c_.get<double>("reference_detector_height");
      bool in_atmosphere = (h1 >= 0 and h2 >= 0);
      bool in_ocean = (h1 < 0 and h2 < 0);
      if (in_atmosphere and fabs(h1-h2) < minimum_height_) {
	h2 -= minimum_height_/2;
	//throw std::runtime_error("Atmospheric detector and reference must be at least 0.01 m appart.");
      }
      else if (in_ocean and fabs(h1-h2) < minimum_depth_) {
	//throw std::runtime_error("Ocean detector and reference must be at least 0.001 m appart.");
	h2 += minimum_depth_/2;
      }
      if (h1 >= 0 and h1 < minimum_height_) {
	h1 += minimum_height_;
      }
      if (h2 >= 0 and h2 < minimum_height_) {
	h2 += minimum_height_;
	//throw std::runtime_error("Atmospheric detector and reference must be at least 0.01 m above surface.");
      }
      if (h1 < 0 and h1 > -minimum_depth_) {
	h1 -= minimum_depth_;
      }
      if (h2 < 0 and h2 > -minimum_depth_) {
	h2 -= minimum_depth_;
	//throw std::runtime_error("Ocean detector and reference must be at least 0.001 m below surface.");
      }
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
      if (h[1] < 0) {
	h_upper = {1};
	h_lower = h;
	n_detector_++;
	n_reference_++;
      }
      else if (h[0] > 0) {
	h_upper = h;
	h_lower = {-bottom_depth_/2};
      }
      else if (h[0] < 0 and h[1] > 0) {
	h_upper = {h[1]};
	h_lower = {h[0]};
      }
      std::reverse(h_upper.begin(), h_upper.end());
      std::reverse(h_lower.begin(), h_lower.end());
      c_.add<double>("DETECTOR_DEPTHS_UPPER_SLAB",max_height_-h_upper);
      c_.add<double>("DETECTOR_DEPTHS_LOWER_SLAB",0-h_lower);
    }
    void run() {
      write(c_,"./tmp",precision_);
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
      ifs.close();
      return t;
    }
    grid_4d read_radiance(const std::string& file_name) {
      std::ifstream ifs(file_name);
      if (!ifs)
	throw std::invalid_argument(file_name+" not found");
      size_t n_runs, n_streams;
      ifs >> n_runs >> n_streams;
      grid_4d g;
      ifs >> g;
      ifs.close();
      return g;
    }
  };

  namespace configuration_template {
    class basic_accurt : public basic_configuration {
      const size_t precision_ = 12;
    protected:
      const double toa_height_ = 120e3;
    public:
      basic_accurt() {
	add_configuration(accurt::configuration());
	set<double>("DETECTOR_WAVELENGTHS",{350e-9, 400e-9, 450e-9, 500e-9,
					    550e-9, 600e-9, 650e-9, 700e-9,
					    750e-9});
	set<double>("reference_detector_height",toa_height_);
	set<double>("SOURCE_ZENITH_ANGLE",45);
	set_text_qualifiers("#","##");
      }
      void write(const std::string& fname) {
	flick::write(*this, "./"+fname, precision_);
      }
    };

    struct atmosphere : public basic_accurt {
      atmosphere() : basic_accurt() {
	add_configuration(material::atmosphere::configuration());
      }
    };

    struct atmosphere_ocean : public basic_accurt {
      atmosphere_ocean() : basic_accurt() {
	add_configuration(material::atmosphere_ocean::configuration());
      }
    };

    struct toa_reflectance : public atmosphere_ocean {
      toa_reflectance() : atmosphere_ocean() {
	set<double>("detector_height",toa_height_-0.2);
	set<std::string>("detector_orientation","down");
	set<std::string>("detector_type","radiance");
      }
    };
    
    struct boa_transmittance : public atmosphere {
      boa_transmittance() : atmosphere() {
	set<std::string>("material_name","atmosphere");
	set<double>("detector_height",0.2);
	set<std::string>("detector_orientation","up");
	set<std::string>("detector_type","irradiance");
      }
    };
      
    struct rs_reflectance : public atmosphere_ocean {
      rs_reflectance() : atmosphere_ocean() {
	set<double>("reference_detector_height",0.4);
	set<double>("detector_height",0.2);
	set<std::string>("detector_orientation","down");
	set<std::string>("detector_type","radiance");
	set<std::string>("subtract_specular_radiance","true");	
      }
    };
    
    struct ramses_reflectance : public atmosphere_ocean {
      ramses_reflectance() : atmosphere_ocean() {
	double irradiance_top = -0.3;
	double detector_distance = 0.3;
	set<double>("reference_detector_height",irradiance_top);
	set<double>("detector_height",irradiance_top-detector_distance);
	set<std::string>("detector_orientation","down");
	set<std::string>("detector_type","radiance");
      }
    };
  }
}

#endif
