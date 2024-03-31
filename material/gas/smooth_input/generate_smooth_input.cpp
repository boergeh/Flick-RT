#include "../absorption_smoother.hpp"
#include "../../../numeric/units.hpp"
#include "../../../numeric/distribution.hpp"

int main(int argc, char* argv[]) {
  const std::vector<std::string> arg(argv + 1, argv + argc);
  int gas_number = std::stoi(arg.at(0));
  using namespace flick;  
  using namespace units;

  struct meta_data {
    std::string short_name;
    std::vector<std::string> gas_names;
    double surface_T0;
    double surface_P0;
    double relative_TP_step;
    double from_wl;
    double to_wl;
    double band_width_sigma;
    double accuracy;
    size_t number_of_wls() const {
      return log(to_wl / from_wl) / band_width_sigma;
    }
    stdvector center_wls() const {
      return range(from_wl,to_wl,number_of_wls()).logspace();
    }
    std::string full_name(size_t gas_number) const {
      return gas_names.at(gas_number)+"_"+short_name;
    }
  };

  std::vector<meta_data> meta_data_vector;
  meta_data m;

  m.short_name = "uv";
  m.gas_names = {atmospheric_state().gas_names().at(gas_number)};
  m.surface_T0 = constants::T_ntp;
  m.surface_P0 = constants::P_ntp;
  m.relative_TP_step = -0.1;
  m.from_wl = 280e-9;
  m.to_wl = 400e-9;
  m.band_width_sigma = log((300+0.1)/(300-0.1));
  m.accuracy = 0.005;
  //meta_data_vector.push_back(m);

  m.short_name = "uv_vis";
  m.to_wl = 950e-9;
  //meta_data_vector.push_back(m);

  m.short_name = "uv_vis_nir";
  m.band_width_sigma = log((300+1)/(300-1));
  m.to_wl = 1100e-9;
  meta_data_vector.push_back(m);

  m.short_name = "solar";
  m.to_wl = 2500e-9;
  //meta_data_vector.push_back(m);

  m.short_name = "terrestrial";
  m.from_wl = 4e-6;
  m.to_wl = 100e-6;
  //m.band_width_sigma = log((4+0.1)/(4-0.1));
  //meta_data_vector.push_back(m);

  struct source {
    double value(double wl) {
      return 1;
    }
  };

  const auto& v = meta_data_vector;
  stdvector scale_T = {1,1+m.relative_TP_step,1};
  stdvector scale_P = {1,1,1+m.relative_TP_step};
  std::vector<std::string> ext = {"T0_P0","T1_P0","T0_P1"};
  std::cout <<"Total number of wavelengths: "<<m.number_of_wls() << std::endl;
  for (size_t k=0; k < ext.size(); k++) {
    for (size_t i=0; i < v.size(); i++) {
      for (size_t j=0; j < v[i].gas_names.size(); j++) {
	std::cout << v[i].short_name << ":" << std::setprecision(4)<<std::endl;
	double T = v[i].surface_T0*scale_T[k];
	double P = v[i].surface_P0*scale_P[k];
	std::cout << "T: " << T << " K"<<std::endl;
	std::cout << "P: " << P/100 << " hPa" << std::endl;
	
	atmospheric_state state(T, P, 8);
	
	std::cout << v[i].gas_names[j] << ":" << std::endl;
	auto as = absorption_smoother(source(), state, v[i].gas_names[j]);
	as.print_progress(true);
	stdvector c = gaussian_smooth_cross_section(as,
						    v[i].center_wls(),
						    v[i].band_width_sigma,
						    v[i].accuracy);
	std::string name = v[i].full_name(j)+"_"+ext[k];
	write(pl_function(v[i].center_wls(),c),"./"+name+".txt");
      }
    }
  }
  return 0;
}
