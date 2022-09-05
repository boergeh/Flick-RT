#include "../../model/single_layer_slab.hpp"

int main() {
  using namespace flick;  
  absorption_coefficient a{0.01};
  scattering_coefficient b{0.01};
  asymmetry_factor g{0.5};

  model::single_layer_slab slab{thickness{1}};
  slab.fill<material::henyey_greenstein>(a,b,g);    
  slab.set(bottom_albedo{1});
  slab.set(incidence_angle{constants::pi/4});
  slab.set(number_of_packages{10000});
  
  std::cout << "slab heimispherical reflectance: "
	    << slab.hemispherical_reflectance() << std::endl;

  return 0;
}
