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
      : iops_{iops}, wls_{wavelengths} {
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const accurt_user_specified& u) {
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
	  double delta_fit_scaling = p[j][i][0]*4*constants::pi;
	  os << "A_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << a[j][i] << " #\n";
	  os << "S_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << s[j][i]*delta_fit_scaling << " #\n";
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
	add<std::string>("material_name","atmosphere_ocean", R"(Name of the material that gives inherent optical properties used in
the radiative transfer calculations, which could be 'atmosphere' or
'atmosphere_ocean')");

	add<std::string>("subtract_specular_radiance","false", R"(Possible subtraction the nadir radiance specularly reflected off the
water surface, which could be 'true' or 'false'. Note that it should
be set to 'true' when calculating remote sensing reflectance)");
	
	add<double>("detector_height", 120e3, R"(Detector height relative to sea level [m], where a positive value
gives height in atmosphere and a negative value gives depth in the
ocean)");

	add<std::string>("detector_orientation", "up", R"(Detector vertical orientation, looking 'up' or 'down')");
	
	add<double>("detector_orientation_override", {0,0}, R"(Polar and azimuth viewing angles override for slant radiance view
[degrees]. Note that '0 0' will deactivate the override, and a
radiance detector will still be looking vertically up when
detector_orientation is set to 'up' and be looking vertically down
when detector_orientations is set to 'down'. For other override
angles, a polar angle less than 90 degrees will make the detector look
towards the upper heimisphere, and zero degrees azimuth will make it
look towards the azimuth position of the source.)");
	
	add<size_t>("detector_radiance_distribution_override", {0,0}, R"(Number of polar and azimuth viewing angles override for detection of
angular distributed radiance. Note that '0 0' will deactivate the
override, and a radiance detector will still be looking vertically up
when detector_orientation is set to 'up' and be looking vertically
down when detector_orientations is set to 'down'. For other override
number of angles, the first value gives number of polar viewing angles
linearly spaced from 0 to 180 degrees, and the second value gives
number of azimuth viewing angles linearly spaced from -180 to 180
degrees. If detector_radiance_distribution_override is activated,
Flick will return a 3D array with radiance values for each wavelength,
polar angle, and azimuth angel, respectively.)");
	
	add<std::string>("detector_type", "irradiance", R"(Type of radiation to be detected, which could be 'plane_irradiance',
'scalar_irradiance', or 'radiance')");
	
	add<double>("detector_wavelengths", {400e-9,500e-9}, R"(Space-separated list of radiation wavelengths to be calculated and
saved [m]. Note that a wavelength of e.g., 400 nm may be written as
400e-9 m for clarity)");
	
	add<double>("reference_detector_height", 120e3, R"(Radiation is calculated relative to a reference plane
irradiance at a given height [m]. For example, a reference
detector height of 100e3 gives calculated radiation relative to the
top-of-atmosphere irradiance)");
	
	add<std::string>("reference_detector_orientation","up", R"(Reference detector vertical orientation, looking 'up' or 'down')");
			 
	add<double>("source_zenith_angle", 0, R"(Source zenith angle [degrees], where zero gives vertically
downward-directed incident irradiance. The source is usually the solar
beam.)");
	
	add<double>("bottom_boundary_surface_scaling_factor", 1, R"(Darkness scaling of the bottom boundary, where '0' gives no bottom
reflection and '1' gives loamy sand reflection)");
		    
	add<size_t>("stream_upper_slab_size", 34, R"(Number of streams used when solving the radiative transfer equation)");
      }
      size_t to_streams(size_t n_angles) {
	size_t n_streams = pow(n_angles,1/1.6); 
	if (n_streams % 2 == 0) 
	  return n_streams;
	return ++n_streams;
      }
    };
  private:
    basic_configuration c_;
    std::shared_ptr<material::base> material_;
    const std::string path_ = "./.tmp_accurtOutput";
    size_t n_detector_;
    size_t n_reference_;
    double max_height_ = 1;
    double bottom_depth_ = 1;
    stdvector wavelengths_;
    const size_t precision_ = 12;
    const double dh_ = 1e-4;
  public:
    accurt(const basic_configuration& c, std::shared_ptr<material::base> material)
      : c_{c}, material_{material} {
      c_.add<std::string>("source_type","constant_one"); 
      c_.add<double>("source_scaling_factor",1);
      c_.add<std::string>("bottom_boundary_surface","loamy_sand");
      c_.add<double>("stream_lower_slab_parameters",{1,2});
      add_layer_depths();
      c_.add<std::string>("materials_included_upper_slab","user_specified_upper_slab");
      c_.add<std::string>("materials_included_lower_slab","user_specified_lower_slab");
      add_detector_depths();
      c_.add<std::string>("detector_azimuth_angles","0:20:180");
      c_.add<std::string>("detector_polar_angles","0:2:180");   
      c_.add<double>("detector_wavelength_band_widths",{270,0.01,4000,0.01});
      c_.add<std::string>("save_cosine_irradiance","true");
      c_.add<std::string>("save_sine_irradiance","false");
      c_.add<std::string>("save_scalar_irradiance","false");
      c_.add<std::string>("save_radiance","false");
      c_.add<std::string>("save_iops","false");
      c_.add<std::string>("save_bottom_boundary_surface","false");
      c_.add<std::string>("save_material_profile","false");
      c_.add<double>("profile_output_wavelength",500);
      c_.add<std::string>("print_progress_to_screen","false");
      c_.add<size_t>("repeated_run_size",1);
      c_.add<std::string>("use_polarization","false");
      c_.add<std::string>("do_ocean_phasematrix","false");
      c_.add<double>("accuracy",0.0);
      c_.add<double>("mie_calculator",2);
      c_.add<std::string>("do_delta_fit","false");
      c_.add<double>("delta_fit_truncation",3.0);
      c_.add<double>("response_function_type",1);
      c_.add<std::string>("do_spherical_correction","true");
      c_.add<std::string>("do_2d_rough_sea_surface","false");
      c_.add<double>("surface_wind_speed",6);
      c_.add<double>("relative_wind_direction",0);
      c_.add<std::string>("usrang","false");
      c_.add<double>("lpick",0);
      wavelengths_ = c_.get_vector<double>("detector_wavelengths");
      c_.set<double>("detector_wavelengths",wavelengths_*1e9);
    }
    pp_function relative_radiation() {
      std::string t = c_.get<std::string>("detector_type"); 
      if (t=="radiance")
	return relative_radiance();
      else if (t=="scalar_irradiance")
	return relative_scalar_irradiance();
      else if (t=="plane_irradiance")
	return relative_plane_irradiance();
      else
	throw std::runtime_error("detector_type");
    }
    bool radiance_distribution_override() {
      return (c_.get_vector<size_t>("detector_radiance_distribution_override").at(0) > 0);
    }
    bool detector_orientation_override() {
      return (c_.get_vector<double>("detector_orientation_override").at(0) > 1e-5);
    }
    void print_radiance_distribution(const std::vector<std::vector<std::vector<double>>>& r) {
      for (size_t i=0; i<r.size(); ++i) {
	for (size_t j=0; j<r[i].size(); ++j) {
	  for (size_t k=0; k<r[i][j].size(); ++k) {
	    std::cout << r[i][j][k] << " ";
	  }
	  std::cout << '\n';
	}
      }      
    }
    std::vector<std::vector<std::vector<double>>> relative_distributed_radiance() {
      auto n_d = c_.get_vector<size_t>("detector_radiance_distribution_override");
      size_t n_polar = n_d.at(0);
      size_t n_azimuth = n_d.at(1);
      set_distributed_radiance(n_polar,n_azimuth);
      run();
      grid_4d g = read_radiance(path_+"/radiance.txt");
      stdvector wls = g.x[1];
      std::vector<std::vector<std::vector<double>>> L;
      
      L.resize(wls.size());
      for (size_t i=0; i<wls.size(); i++) {
	L[i].resize(n_polar);
	for (size_t j=0; j<n_polar; j++) {
	  L[i][j].resize(n_azimuth);
	}
      }
      for (size_t i=0; i<wls.size(); i++) {
	for (size_t j=0; j<n_polar; j++) {
	  for (size_t k=0; k<n_azimuth; k++) {
	    L[i][j][k] = g.f[n_detector_][i][j][k];
	  }
	}
      }    
      return L;
    }
  private:
    void set_vertical_radiance() {
      c_.set<std::string>("save_radiance","true");
      c_.set<std::string>("detector_azimuth_angles","nan");
      c_.set<std::string>("detector_polar_angles","0 180");
    }
    void set_slant_radiance(double polar_view, double azimuth_view) {
      const double polar_max = 179.1;
      if (polar_view > polar_max)
	polar_view = polar_max;
      c_.set<std::string>("save_radiance","true");
      c_.set<std::string>("detector_polar_angles","0 "+std::to_string(180-polar_view));
      c_.set<std::string>("detector_azimuth_angles",std::to_string(azimuth_view));
    }
    void set_distributed_radiance(size_t n_polar, size_t n_azimuth) {
      c_.set<std::string>("save_radiance","true");
      std::string polar = vector_to_string(range(0,180,n_polar).linspace());
      std::string azimuth = vector_to_string(range(-180,180,n_azimuth).linspace());
      c_.set<std::string>("detector_polar_angles",polar);
      c_.set<std::string>("detector_azimuth_angles",azimuth);
    }
    std::string vector_to_string(const stdvector& v) {
      std::string s;
      for (size_t i=0; i<v.size(); ++i) {
	s += std::to_string(v[i]) + " ";
      }
      return s;
    }
    void make_material_files() {        
      size_t n_terms = c_.get<size_t>("stream_upper_slab_size") + 1;

      stdvector b = depths_to_boundaries(c_.get_vector<double>("layer_depths_upper_slab"));
      auto layered_upper_slab = std::make_shared<layered_iops>(material_,b,n_terms);
      write(accurt_user_specified(layered_upper_slab, wavelengths_),
       	    "./.tmp_accurtMaterials/user_specified_upper_slab", precision_);

      b = 0 - c_.get_vector<double>("layer_depths_lower_slab");
      std::reverse(b.begin(),b.end());
      b.push_back(0);
      auto layered_lower_slab = std::make_shared<layered_iops>(material_,b,n_terms);
      write(accurt_user_specified(layered_lower_slab, wavelengths_),
      	    "./.tmp_accurtMaterials/user_specified_lower_slab", precision_);
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
      if (detector_orientation_override()) {
	auto view = c_.get_vector<double>("detector_orientation_override");     
	double polar_view = view.at(0);
	double azimuth_view = view.at(1);
	set_slant_radiance(polar_view, azimuth_view);
      } else {
	set_vertical_radiance();
      }
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
	Lu[i] = g.f[n_detector_][i][1][0];
	Ld[i] = g.f[n_detector_][i][0][0];
      }
      std::string subtract = c_.get<std::string>("subtract_specular_radiance"); 
      if (subtract=="true") {
	Lu =  Lu - Ld * nadir_fresnel_coefficient(wls);
      } else if (subtract != "false") {
	throw std::runtime_error("subtract_specular_radiance");
      }
      if (detector_orientation_override()) {
	return Ld;
      }
      std::string ori = c_.get<std::string>("detector_orientation");
      if (ori=="up") {
	return Ld;
      } else if (ori=="down") {
	return Lu;
      } else {
	throw std::runtime_error("detector_orientation");
      }
    }
    stdvector detector_plane_irradiance() {
      pe_table t;
      if (c_.get<std::string>("detector_orientation")=="up") {
	t = read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      } else if (c_.get<std::string>("detector_orientation")=="down") {
	t = read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      } else
	throw std::runtime_error("detector_orientation");
      return t.row(n_detector_).y();
    }
    stdvector detector_scalar_irradiance() {
      c_.set<std::string>("SAVE_SCALAR_IRRADIANCE","true");      
      pe_table t;
      if (c_.get<std::string>("detector_orientation")=="up") {
	t = read_irradiance(path_+"/scalar_irradiance_total_downward.txt");
      } else if (c_.get<std::string>("detector_orientation")=="down") {
	t = read_irradiance(path_+"/scalar_irradiance_total_upward.txt");
      } else
	throw std::runtime_error("detector_orientation");
      return t.row(n_detector_).y();
    }
    stdvector reference_detector_irradiance() {
      pe_table t;
      if (c_.get<std::string>("reference_detector_orientation")=="up") {
	t = read_irradiance(path_+"/cosine_irradiance_total_downward.txt");
      } else if (c_.get<std::string>("reference_detector_orientation")=="down") {
	t = read_irradiance(path_+"/cosine_irradiance_total_upward.txt");
      } else
	throw std::runtime_error("reference_detector_orientation");
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
      stdvector h = {-bottom_depth_, 0, max_height_};
      material::z_profile<pe_function>* zp = dynamic_cast<material::z_profile<pe_function>*>(&*material_);
      if (zp != NULL) {
	h = zp->height_grid();
      }      
      max_height_ = h.back();
      if (h.front()>=0)
      	h.insert(h.begin(), -bottom_depth_);
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
      depths_upper.back() = max_height_; 
      c_.add<double>("layer_depths_upper_slab", depths_upper);
      c_.add<double>("layer_depths_lower_slab", depths_lower);
    }
    double ensure_distance_to_surface(double h) {
      if (h >= 0 and h < dh_)
	h = dh_;
      if (h < 0 and h > -dh_)
	h = -dh_;
      return h;
    }
    double ensure_distance_to_max_height(double h) {
      if (h >=  max_height_)
	h = max_height_-dh_;
      return h;
    }
    double ensure_distance_to_bottom(double h) {
      if (h <=  -bottom_depth_)
	h = -bottom_depth_+dh_;
      return h;
    }
    void add_detector_depths() {
      double h_d = c_.get<double>("detector_height");    
      double h_r = c_.get<double>("reference_detector_height");
      h_d = ensure_distance_to_surface(h_d);
      h_r = ensure_distance_to_surface(h_r);
      h_d = ensure_distance_to_max_height(h_d);
      h_r = ensure_distance_to_max_height(h_r);
      h_d = ensure_distance_to_bottom(h_d);
      h_r = ensure_distance_to_bottom(h_r);
      
      bool both_in_atmosphere = (h_d > 0 and h_r > 0);
      bool both_in_ocean = (h_d < 0 and h_r < 0);
      bool same_height = (fabs(h_d-h_r) < dh_);

      double mh = max_height_;
      stdvector depths = {0, mh-dh_/2, mh+dh_/2, mh+bottom_depth_-dh_/2};      
      if (h_d > h_r) {
	n_detector_ = 0;
	n_reference_ = 1;
      } else {
	n_detector_ = 1;
	n_reference_ = 0;
      }
      if (not (both_in_atmosphere or both_in_ocean)) {
	n_detector_ += 1;
	n_reference_ += 1;
      } else if (both_in_ocean) {
	n_detector_ += 2;
	n_reference_ += 2;
      }
      if (same_height && both_in_atmosphere) {
	n_reference_ = 0;
	n_detector_ = 0;
	depths[1] = mh-dh_/2;
      } else if (same_height && both_in_ocean) {
	n_reference_ = 3;
	n_detector_ = 3;
	depths[2] = mh+dh_/2;
      }
      depths[n_detector_] = mh - h_d;
      depths[n_reference_] = mh - h_r;
      c_.add<double>("detector_depths_upper_slab",{depths[0],depths[1]});
      c_.add<double>("detector_depths_lower_slab",{depths[2]-mh,depths[3]-mh});
    }
    void run() {
      c_.set_text_qualifiers("#","##");
      c_.set_uppercase(true);
      write(c_,"./.tmp_accurt",precision_);
      system("mkdir -p ./.tmp_accurtMaterials");
      system("mkdir -p ./.tmp_accurtOutput");
      make_material_files();
      int s=system("DYLD_LIBRARY_PATH=$ACCURT_PATH/lib AccuRT .tmp_accurt");
      if (s!=0)
	throw std::runtime_error("accurt_api");
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

    struct ocean_radiance : public atmosphere_ocean {
      ocean_radiance() : atmosphere_ocean() {
	set<double>("detector_height",-0.5);
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
