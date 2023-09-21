#include "atmosphere.hpp"
#include "../numeric/units.hpp"
#include "ab_functions.hpp"
#include "../numeric/legendre/delta_fit.hpp"
#include "aerosols/aerosols.hpp"

namespace flick {
  begin_test_case(atmosphere_test) {
    using namespace units;
    using namespace material;
   
    atmosphere::config c;
    c.set<size_t>("angles",100);
    c.set<size_t>("heights",8);
    auto atm = atmosphere(c);
    //auto wls = range(280e-9,950e-9,2).logspace();
    //std::cout << std::setprecision(5)
    //	      << optical_depth(*atm,100e3,wls).scattering() << std::endl;

    atm.set_position({0,0,0});
    double s1 = atm.scattering_coefficient();
    atm.set_position({0,0,5e3});
    double s2 = atm.scattering_coefficient();
    check(s1 < 0.1*s2);

    //auto [alpha,beta] = material::fitted_mueller_alpha_beta(atm,8);
    //std::cout << "Phase function terms:" << alpha[0]/(4*constants::pi) << std::endl;
    auto f = phase_function(atm);
    //std::cout << "Phase function terms: "<< std::setprecision(5) << delta_fit(f,32).coefficients()*(4*constants::pi);
  } end_test_case()
}
