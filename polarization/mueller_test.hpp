#include "algorithm.hpp"
#include "isotropic_mueller.hpp"
#include "rayleigh_mueller.hpp"
#include "mueller.hpp"
#include "../environment/input_output.hpp"

namespace flick {
  begin_test_case(mueller_test_A) {    
    using namespace constants;  
    pl_function f;
    size_t n = 150;
    auto ang = range(0,pi,n).linspace();
    for (size_t i=0; i<n; ++i) {
      mueller rm = rayleigh_mueller(ang[i],0.2);
      f.append({ang[i], rm(0).value * sin(ang[i])});
    }
    check_close(2*pi*f.integral(), 1, 0.01_pct);
    angular_mueller am(hg_phase_function<pe_function>(0.9,100));
    
    check_close(am.value(0,0,1),15.12,1.0_pct);
    check_small(am.value(1,0,0));

    angular_mueller am2(hg_phase_function<pe_function>(0.0,100));
    am2.add(3,3,pl_function{{0,3.14159},{1,1}});
   
    check_close(am2.value(3,3,0),1,1e-9_pct);
    am.add(am2,0.5);
    check_close(am.value(3,3,0),0.5,3e-3_pct);
     
  } end_test_case()
    
  begin_test_case(mueller_test_B) {   
    double pi = constants::pi;
    double g = 0.95;
    pe_function hg = hg_phase_function<pe_function>(g,60,0.9);
    tabulated_phase_function p{hg};
    check_close(p.asymmetry_factor(), g, 1_pct);
    check_close(p.integral_4pi(), 1, 1_pct);
        
    const std::string path_{"/polarization"};
    pe_function  coestner = read<pe_function>(path_+"/tabulated_phase_function.txt");
    coestner = to_cos_x(coestner);
    tabulated_phase_function p2{coestner};
    check_close(p2.integral_4pi(),1,0.4_pct);
    check_close(p2.asymmetry_factor(), 0.95, 0.1_pct);
  } end_test_case()
  
  begin_test_case(mueller_test_C) {   
    pe_function hg1 = hg_phase_function<pe_function>(0.9,10,0.7);
    tabulated_phase_function p1{hg1};
    check_close(p1.y().back(),hg1.y().back());
    angular_mueller am(p1);
    std::cout << am.value(0,0,1);
    check_close(am.value(0,0,1),hg1.y().back());
    angular_mueller am2 = am;
    am.add(am2, 0.3);  
    check_close(am.value(0,0,1),hg1.y().back());
    pe_function hg2 = hg_phase_function<pe_function>(0,10,0.7);
    tabulated_phase_function p2{hg2};
    am.add(angular_mueller(p2),0.5);
    check_close(am.value(0,0,1),0.5*(hg1.y().back()+hg2.y().back()));  
  } end_test_case()
}

