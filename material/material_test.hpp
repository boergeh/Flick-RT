#include "water/pure_water.hpp"
#include "henyey_greenstein.hpp"
#include "tabulated.hpp"
#include "z_profile.hpp"
#include "../numeric/range.hpp"
#include "aerosols/aerosols.hpp"

namespace flick {
  begin_test_case(material_test) {
    using namespace constants;
    using namespace flick;
    material::pure_water pw;
    material::henyey_greenstein hg(absorption_coefficient{1},
				   scattering_coefficient{1},
				   asymmetry_factor{0.5});
    double pi = constants::pi;
    double g = 0.9;
    phase_function p = hg_phase_function(g,100);
    check_close(2*pi*p.integral(),1,0.03,"a");
    check_close(p.asymmetry_factor(),g,0.03,"b");
    
    material::tabulated tab(absorption_coefficient{1},
    			    scattering_coefficient{1},
    			    p);

    const std::string path_{"/material"};
    phase_function p2 = read<pe_function>(path_+"/tabulated.txt"); 
    check_close(2*pi*p2.integral(),1,0.3,"c");

    material::rural_aerosols ra;
    material::aggregate_z_profile ag({0,1.5e3,100e3},{0,1,3.14});
    ag.set(pose{{0,0,0},unit_vector{0,0,1}});
    ag.add(ra);
    check_close(ag.scattering_optical_depth(100e3),
		ra.scattering_optical_depth(100e3),1e-9,"a");
    unit_vector u{0,0};
    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0),1e-9,"b");
    ag.add(ra);
    check_close(ag.absorption_optical_depth(100e3),
    		2*ra.absorption_optical_depth(100e3),1e-9,"c");
    ag.add(ag);
    check_close(ag.absorption_optical_depth(100e3),
    		4*ra.absorption_optical_depth(100e3),1e-9,"d");

    u = unit_vector{3.14,0};
    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0),1e-9,"e");
    material::urban_aerosols ua;
    ag.add(ua);
    check_close(ag.absorption_optical_depth(100e3),
    		4*ra.absorption_optical_depth(100e3)+
		ua.absorption_optical_depth(100e3),1e-9,"f");

    check_close(ag.mueller_matrix(u).value(0,0),
		ra.mueller_matrix(u).value(0,0)*0.999,0.07,"g");
   
    //auto m = pw.mueller_matrix();
    //std::cout << m.element(0,0);
    /*
    water w;
    w.wavelength(500e-9);
    check_close(w.abs_coef(),0.02,5);
    w.salinity(pl_function{35});
    w.temperature({273});
    check_close(w.abs_coef(),0.021,2);
    w.wavelength(835e-9);
    check_close(w.abs_coef(),2.99,1);
    w.salinity(pl_function{0});
    check_close(w.abs_coef(),3.02,1);
    w.temperature({273+30});
    check_close(w.abs_coef(),3.47,1);
    */
    //w.temperature({273});
    //check_close(w.abs_coef(),0.02,1);
    /*
    pp_function fm = read<pp_function>("/material/water/temperature_correction.txt");
    fm.scale_x(1e-9);
    std::cout << fm;
    write<pp_function>(fm,"/material/water/temperature_correction_new.txt",6);

    pp_function fm2 = read<pp_function>("/material/water/salinity_correction.txt");
    //fm2.scale_x(1e-9);
    //std::cout << fm2;
    write<pp_function>(fm2,"/material/water/salinity_correction_new.txt",6);
    pp_function fm3 = significant_digits(fm2, 10, 4);
    std::cout << fm3;
    */
    /*
    double epsilon = 1e-12;
    constant_profile cp{1};
    exponential_profile ep{{0,1},{5.9e3,0.5}};
    pose horizontal{{0,0,0},{1,0,0}};
    pose vertical{{0,0,0},{0,0,1}};
    check_close(cp.column_density(7,horizontal),
		ep.column_density(7,horizontal),epsilon);
    double scale_height = ep.column_density(100e3,vertical); 
    std::optional<double> column_length = ep.column_length(scale_height,vertical); 
    check_close(scale_height, 8.5e3, 1);
    check_close(*column_length, 100e3, 1e-10);
    */
  } end_test_case()
}
