#ifndef flick_lines
#define flick_lines

#include <fstream>
#include <numbers>
#include "../../environment/input_output.hpp"
#include "../../numeric/constants.hpp"
#include "../../numeric/range.hpp"
#include "../../numeric/sorted_vector.hpp"


namespace flick {  
  class lines
  // Implemented from https://hitran.org
  // Note: background spectrum for water vapor not included 
  {    
  public:
    lines(const std::string& gas_file_name) {
      read_lines_input(path()+"/material/gas/lines_input/"+gas_file_name);
    }
    void total_pressure(double p) {
      total_pressure_ = p;
    }
    void partial_pressure(double p) {
      partial_pressure_ = p;
    }
    void wing_cutoff(double wavenumber_per_cm)
    {
      wing_cutoff_ = wavenumber_per_cm;
    }
    void temperature(double t) {
      temperature_ = t;
    }
    static double molecules_per_volume(double pressure, double temperature)
    // Ideal gas law
    {
      return pressure / constants::k_B / temperature;
    }    
    double absorption_coefficient(double wavelength) const {
      double readout_wavenumber = 1/(wavelength*1e2);
      int above = 1;
      int below = -1;
      size_t close_line = wavenumbers_.find(readout_wavenumber);
      double abs = absorption_from(close_line, readout_wavenumber);
      abs += add_line_wings(above, close_line, readout_wavenumber);
      abs += add_line_wings(below, close_line, readout_wavenumber);
      return abs;
    }
    std::vector<double> absorption_coefficient(const std::vector<double>&
					       wavelengths) const
    // Returns absorption coefficient spectrum [1/m] 
    {
      std::vector<double> v(wavelengths.size());
      for (size_t i = 0; i < wavelengths.size(); ++i) {
	v[i] = absorption_coefficient(wavelengths[i]);
      }
      return v;
    }
  private:
    // All wavenumbers are internally stored in units of 1/cm
    struct line_parameters {
      double line_intensity;
      float air_half_width;
      float self_half_width;
      float temperature_exponent;
      float line_shift;
      float lower_state_energy;
    };
    const double reference_temperature_{296};
    double total_pressure_{constants::P_ntp};
    double partial_pressure_{constants::P_ntp};
    double wing_cutoff_{250};
    double temperature_{300};   
    sorted_vector wavenumbers_;
    std::vector<line_parameters> lines_;

    double add_line_wings(int step, size_t close_line,
			  double readout_wavenumber) const {
      double total_abs = 0;
      double wing_abs = 1;
      while(wing_abs > 0) {
	close_line += step;
	wing_abs = absorption_from(close_line,readout_wavenumber);
	total_abs += wing_abs;
      }
      return total_abs;
    }
    double absorption_from(int line_number,
			   double readout_wavenumber) const {
      if (line_number < 0)
	return 0;
      if (line_number >= lines_.size())
	return 0;
      double nu0 = wavenumbers_[line_number];
      double nu = readout_wavenumber;
      if (fabs(nu-nu0) > wing_cutoff_)
      	return 0;      
      double p = total_pressure_ / constants::P_ntp;
      double p_self = partial_pressure_ / constants::P_ntp;
      double T = temperature_;
      double T_ref = reference_temperature_;
      double S = lines_.at(line_number).line_intensity;
      double gam_air = lines_.at(line_number).air_half_width;
      double gam_self = lines_.at(line_number).self_half_width;
      double E_low = lines_.at(line_number).lower_state_energy;
      double n_air = lines_.at(line_number).temperature_exponent;
      double delta_air = lines_.at(line_number).line_shift;
      double gam = pow(T_ref/T,n_air)*(gam_air*(p-p_self)+gam_self*p_self);
      double nu_star = nu0+delta_air*p;
      double f_L = gam/std::numbers::pi/(pow(gam,2)+pow(nu-nu_star,2));
      if (nu < 200) { // Van Vleck-Weisskopf lineshape
	f_L = gam/std::numbers::pi*nu/nu0*
	  (1/(pow(gam,2)+pow(nu-nu_star,2))+1/(pow(gam,2)+pow(nu+nu_star,2)));
      }
      double cm2_to_m2 = 1e-4;
      double N = molecules_per_volume(partial_pressure_, T);
      return S * f_L * cm2_to_m2 * N;
    }   
    void read_lines_input(const std::string& file) {
      std::ifstream ifs(file);
      if (!ifs)
	throw std::invalid_argument(file+" not found");
      size_t rows;
      ifs >> rows;
      lines_.reserve(rows);
      std::vector<double> waven_vec;
      waven_vec.reserve(rows);
      line_parameters lp;
      double wn;
      double epsilon = 1e-6;
      for (size_t i=0; i<rows; ++i) {
	ifs >> wn;
	if (i > 0 && wn < waven_vec[i-1]+epsilon) {
	  wn = waven_vec[i-1]+epsilon;
	}
	waven_vec.push_back(wn);	
	ifs >> lp.line_intensity;
	ifs >> lp.air_half_width;
	ifs >> lp.self_half_width;
	ifs >> lp.temperature_exponent;
	ifs >> lp.line_shift;
	ifs >> lp.lower_state_energy;
	if (!(wn > 0) || !(lp.line_intensity > 0)
	    || !(lp.air_half_width > 0) || !(lp.self_half_width > 0))
	  throw std::invalid_argument(file+" in line "+std::to_string(i+2));
	lines_.push_back(lp);
      }
      ifs.close();
      wavenumbers_ = sorted_vector(waven_vec);
    }    
  };
}

#endif
  
