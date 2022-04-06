#include "material.hpp"
namespace flick {
  begin_test_case(material_test) {
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
