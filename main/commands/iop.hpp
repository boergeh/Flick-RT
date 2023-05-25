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
	}
	else if (a(5)=="ice_cloud" or a(5)=="snow") {
	  double mu = log(std::stod(a(8)));
	  double sigma = std::stod(a(9));
	  double volfrac = std::stod(a(7));
	  if (a(6)=="full_mie") {
	    material::ice_cloud<monodispersed_mie> m(volfrac,mu,sigma);
	   stream_iops(m, a(1));
	  } else {
	    material::ice_cloud<parameterized_monodispersed_mie>
	      m(volfrac,mu,sigma);
	    stream_iops(m, a(1));
	  }
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
      std::vector<size_t> sub_script_numbers(const std::string& s) const {
	std::vector<size_t> numbers;
	size_t i1 = s.find("_",0);
	size_t i2 = s.find("_",i1+1);
	while(i2 != std::string::npos) {
	  numbers.push_back(std::stoi(s.substr(i1+1,i2-i1-1)));
	  i1 = s.find("_",i1+1);
	  i2 = s.find("_",i2+1);
	}
	numbers.push_back(std::stoi(s.substr(i1+1,s.size()-1)));
	return numbers;
      }
      template<class Material>
      void stream_iops(Material& m, const std::string& property) const {
	const std::string& p = property;
	if (p.substr(0,6)=="AccuRT") {
	  auto n = sub_script_numbers(p);
	  stream_AccuRT(m, n.at(0), n.at(1));
	} else if (p.substr(0,13)=="scattering_ab") {
	  auto n = sub_script_numbers(p.substr(13));	  
	  size_t n_angles = n.at(0);
	  size_t n_terms = 0;
	  auto [a,b,x] = m.mueller_ab_functions(n_angles);
	  if (n.size()==2) {
	    n_terms = n.at(1);
	    std::tie(a,b,x) = m.fitted_mueller_ab_functions(n_angles,n_terms);
	  }
	  std::cout << std::setprecision(7);
	  for (size_t i = 0; i < wls_.size(); ++i) {
	    m.set(wavelength(wls_[i]));
	    double k = m.scattering_coefficient();
	    for (size_t j = 0; j < x.size(); ++j) {
	      std::cout << x[j] << " " << a[0][j]*k << " " << a[1][j]*k
			<< " " << a[2][j]*k
			<< " " << a[3][j]*k
			<< " " << b[0][j]*k
			<< " " << b[1][j]*k << "\n";
	    }
	  }
	} else if (p.substr(0,17)=="wigner_alpha_beta") {
	  auto n = sub_script_numbers(p.substr(17));
	  size_t n_terms = n.at(0);
	  auto [alpha,beta] = m.fitted_mueller_alpha_beta(n_terms);
	  std::cout << std::setprecision(7);
	  for (size_t i = 0; i < wls_.size(); ++i) {
	    m.set(wavelength(wls_[i]));
	    double k = m.scattering_coefficient();
	    for (size_t j = 0; j < n_terms; ++j) {
	      std::cout << alpha[0][j]*k << " " << alpha[1][j]*k
			<< " " << alpha[2][j]*k
			<< " " << alpha[3][j]*k
			<< " " << beta[0][j]*k
			<< " " << beta[1][j]*k << "\n";
	    }
	  }
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
	    std::cout << std::setprecision(6) << wl << " " << value << '\n';
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
