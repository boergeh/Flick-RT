#include "../absorption_smoother.hpp"
#include "../../../numeric/units.hpp"
#include "../../../numeric/distribution.hpp"

int main(int argc, char* argv[]) {
  const std::vector<std::string> arg(argv + 1, argv + argc);
  std::string wl_region_name = arg.at(0);
  std::string gas_name = arg.at(1);
  using namespace flick;  
  using namespace units;

  struct meta_data {
    std::string short_name;
    std::vector<std::string> gas_names;
    const double relative_TP_step = -0.1; // must be same as in smooth_air class
    double surface_T0;
    double surface_P0;
    double from_wl;
    double to_wl;
    double band_width_sigma;
    double accuracy;
    double length_factor;
    size_t number_of_wls() const {
      return 2*log(to_wl / from_wl) / band_width_sigma;
    }
    stdvector center_wls() const {
      return range(from_wl,to_wl,number_of_wls()).logspace();
    }
    std::string full_name(const std::string& gas_name) const {
      return gas_name+"_"+short_name;
    }
    double get_band_width_sigma(double wl0, double delta_wl) {
      return log((wl0+delta_wl/2)/(wl0-delta_wl/2));
    } 
  };

  std::vector<meta_data> meta_data_vector;
  meta_data m;

  m.short_name = wl_region_name;
  m.gas_names = {gas_name};
  m.surface_T0 = constants::T_ntp;
  m.surface_P0 = constants::P_ntp;
  m.accuracy = 0.01;
  m.length_factor = 1.5;

  if (wl_region_name=="uv") {
    m.from_wl = 280e-9;
    m.to_wl = 400e-9;
    m.length_factor = 1.5;
    m.band_width_sigma = m.get_band_width_sigma(300,0.2);
  }
  else if (wl_region_name=="uv_vis") { 
    m.from_wl = 280e-9;
    m.to_wl = 1030e-9;
    m.band_width_sigma = m.get_band_width_sigma(300,0.5);
  }
  else if (wl_region_name=="uv_vis_toa") {
    m.from_wl = 280e-9;
    m.to_wl = 1030e-9;
    m.band_width_sigma = m.get_band_width_sigma(300,0.5);
    m.length_factor = 2.5;
  }
  else if (wl_region_name=="solar") {
    m.from_wl = 280e-9;
    m.to_wl = 2500e-9;
    m.band_width_sigma = m.get_band_width_sigma(300,5);
  }
  else if (wl_region_name=="terrestrial") {
    m.from_wl = 4e-6;
    m.to_wl = 100e-6;
    m.band_width_sigma = m.get_band_width_sigma(4000,10);
  } else {
    throw std::runtime_error("unknown spectral region name");
  }
  struct source {
    pp_function F = read<pp_function>("./toa_solar_spectrum.txt");
    double value(double wl) {
      return F.value(wl);
    }
  };

  const auto& v = meta_data_vector;
  stdvector scale_T = {1,1+m.relative_TP_step,1};
  stdvector scale_P = {1,1,1+m.relative_TP_step};
  std::vector<std::string> ext = {"T0_P0","T1_P0","T0_P1"};
  std::cout <<"Total number of wavelengths: "<<m.number_of_wls() << std::endl;
  for (size_t k=0; k < ext.size(); k++) {
    for (size_t j=0; j < m.gas_names.size(); j++) {
      std::cout << m.short_name << ":" << std::setprecision(4)<<std::endl;
      double T = m.surface_T0*scale_T[k];
      double P = m.surface_P0*scale_P[k];
      std::cout << "T: " << T << " K"<<std::endl;
      std::cout << "P: " << P/100 << " hPa" << std::endl;
      
      atmospheric_state state(T, P, 8);
      
      std::cout << m.gas_names[j] << ":" << std::endl;
      auto as = absorption_smoother(source(), state, m.gas_names[j]);
      as.set_length_factor(m.length_factor);
      as.print_progress(true);
      stdvector c = triangular_smooth_cross_section(as,
						    m.center_wls(),
						    m.band_width_sigma,
						    m.accuracy);
      std::string name = m.full_name(gas_name)+"_"+ext[k];
      write(pl_function(m.center_wls(),c),"./"+name+".txt");
    }
  }
  return 0;
}
