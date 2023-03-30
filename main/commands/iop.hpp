#ifndef flick_command_iop
#define flick_command_iop

#include "basic_command.hpp"
#include "../../material/water/pure_water.hpp"
#include "../../material/ice/pure_ice.hpp"
#include "../../material/spheres.hpp"
#include "../../material/henyey_greenstein.hpp"
#include "../../material/tabulated.hpp"
#include "../../numeric/legendre/delta_fit.hpp"
#include "../../environment/input_output.hpp"

namespace flick {
  namespace command {
    class iop : public basic_command {
      std::vector<double> wls_;
    public:
      iop():basic_command("iop"){};
      void run() {
	double from_wl = std::stod(a(2));
	double to_wl = std::stod(a(3));
	double n_points = std::stod(a(4));	
	wls_ = range(from_wl, to_wl, n_points).logspace();
	
	if (a(5)=="pure_water") {
	  double T = std::stod(a(6));
	  double S = std::stod(a(7));	
	  material::pure_water m;
	  m.temperature(T);
	  m.salinity(S);
	  stream_iops(m, a(1));
	}
	else if (a(5)=="pure_ice") {
	  material::pure_ice m;
	  stream_iops(m, a(1));
	}
	else if (a(5)=="water_cloud") {
	  double mu = log(std::stod(a(8)));
	  double sigma = std::stod(a(9));
	  double volfrac = std::stod(a(7));
	  if (a(6)=="full_mie") {
	    material::water_cloud<monodispersed_mie> m(volfrac,mu,sigma);
	   stream_iops(m, a(1));
	  } else {
	    material::water_cloud<parameterized_monodispersed_mie>
	      m(volfrac,mu,sigma);
	    stream_iops(m, a(1));
	  }
	  /*
	  log_normal_distribution sd(log(std::stod(a(8))),std::stod(a(9))); 
	  material::spheres<log_normal_distribution,
			    material::vacuum,
			    material::pure_water,
			    parameterized_monodispersed_mie>
	    m(std::stod(a(7)), sd, material::vacuum(), material::pure_water());
	  */
	  //stream_iops(m, a(1));
	}
	else if (a(5)=="bubbles_in_ice") {
	  double mu = log(std::stod(a(8)));
	  double sigma = std::stod(a(9));
	  double volfrac = std::stod(a(7));
	  if (a(6)=="full_mie") {
	    material::bubbles_in_ice<monodispersed_mie> m(volfrac,mu,sigma);
	   stream_iops(m, a(1));
	  } else {
	    material::bubbles_in_ice<parameterized_monodispersed_mie>
	      m(volfrac,mu,sigma);
	    stream_iops(m, a(1));
	  }
	  
	  //material::bubbles_in_ice<parameterized_monodispersed_mie> m(volfrac,r_mean,width);
	  /*
	  log_normal_distribution sd(log(std::stod(a(8))),std::stod(a(9)));
	  material::spheres<log_normal_distribution,
			    material::pure_ice,
			    material::vacuum,
			    monodispersed_mie>
	    m(std::stod(a(7)), sd, material::pure_ice(), material::vacuum());
	  */
	  //stream_iops(m, a(1));
	}
	else if (a(5)=="brines_in_ice") {
	  double mu = log(std::stod(a(8)));
	  double sigma = std::stod(a(9));
	  double volfrac = std::stod(a(7));
	  double salinity = std::stod(a(10));
	  if (a(6)=="full_mie") {
	    material::brines_in_ice<monodispersed_mie> m(volfrac,mu,
							 sigma,salinity);
	    stream_iops(m, a(1));
	  } else {
	    material::brines_in_ice<parameterized_monodispersed_mie>
	      m(volfrac,mu, sigma,salinity);
	    stream_iops(m, a(1));
	  }
	  
	  /*
	  log_normal_distribution sd(log(std::stod(a(8))),std::stod(a(9)));
	  double T = 273.15;
	  double S = std::stod(a(10));	
	  material::pure_water w;
	  w.temperature(T);
	  w.salinity(S);
	  material::spheres<log_normal_distribution,
			    material::pure_ice,
			    material::pure_water,
			    parameterized_monodispersed_mie>
	    m(std::stod(a(7)), sd, material::pure_ice(), w);
	  */
	  
	  //stream_iops(m, a(1));
	}
	else if (a(5)=="henyey_greenstein") {
	  absorption_coefficient abs{std::stod(a(6))};
	  scattering_coefficient sca{std::stod(a(7))};
	  asymmetry_factor g{std::stod(a(8))};
	  double real_refractive_index{std::stod(a(9))};
	  material::henyey_greenstein m(abs,sca,g,real_refractive_index);
	  stream_iops(m, a(1));
	}
	else if (a(5)=="tabulated") {
	  absorption_coefficient abs{std::stod(a(6))};
	  scattering_coefficient sca{std::stod(a(7))};
	  std::string  file_name{a(8)};
	  double real_refractive_index{std::stod(a(9))};
	  tabulated_phase_function p = read<pe_function>("./"+file_name);
	  material::tabulated m(abs,sca,p,real_refractive_index);
	  stream_iops(m, a(1));
	}
	else
	  error();
      }
    private:
      template<class Material>
      void stream_iops(Material& m, const std::string& property) const {
	const std::string& p = property;
	if (p.substr(0,6)=="AccuRT") {
	  size_t i1 = p.find("_");
	  size_t i2 = p.rfind("_");
	  size_t layer_n = std::stoi(p.substr(i1+1,i2-i1-1));
	  size_t n_terms = std::stoi(p.substr(i2+1,p.size()));
	  stream_AccuRT(m,layer_n, n_terms);
	} else {
	  double value = 0;
	  for (auto wl:wls_) {
	    m.set(wavelength(wl));
	    if (p=="absorption_length")
	      value = 1/m.absorption_coefficient();
	    else if (p=="scattering_length")
	      value = 1/m.scattering_coefficient();
	    else if (p=="refractive_index")
	      value = m.real_refractive_index();
	    else {
	      error();
	      break;
	    }
	    std::cout << std::setprecision(5) << wl << " " << value << '\n';
	  }
	}
      }
      template<class Material>
      void stream_AccuRT(Material& m, size_t layer_n, size_t n_terms) const {
	std::cout
	  << "# AccuRT configuration file for a one-layer user-specified material #\n"
	  << "PROFILE_LABEL = layer_numbering #\n"
	  << "MATERIAL_PROFILE = " << layer_n << " 1 #\n"
	  << "WAVELENGTHS = ";
	for (auto wl:wls_)
	  std::cout << wl*1e9 << " ";
	std::cout << "#\n";
	delta_fit df(material::phase_function(m),n_terms);
	auto normalized_scaled_coef = df.coefficients()/df.scaling_factor()
	  * 4 * constants::pi;
	for (size_t i=0; i<wls_.size(); ++i) {
	  m.set(wavelength(wls_[i]));
	  std::cout << "A_"<<layer_n<<"_" << std::to_string(i+1) << " = "
		    << m.absorption_coefficient() << " #\n"
		    << "S_"<<layer_n<<"_" << std::to_string(i+1) << " = "
		    << m.scattering_coefficient()*df.scaling_factor() << " #\n"
		    << "P_"<<layer_n<<"_" <<  std::to_string(i+1) << " = ";
	  for (size_t j=0; j < normalized_scaled_coef.size(); ++j) {
	    std::cout << normalized_scaled_coef[j]/(2*j+1) << " "; 
	  }
	  std::cout << " #\n";
	}
	std::cout << "REFRACTIVE_INDICES = ";
	for (size_t i=0; i<wls_.size(); ++i) {
	  m.set(wavelength(wls_[i]));
	  std::cout << m.real_refractive_index() << " ";
	}
	std::cout << "#\n";
	std::cout << "TURN_OFF_DELTA_FIT = true #\n";
      }
    };    
  }
}

#endif
