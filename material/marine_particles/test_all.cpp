#include "../../environment/unit_test.hpp"
#include "marine_particles_test.hpp"

int main() {
  using namespace flick;
  unit_test t("marine_particles");
  t.include<marine_particles_test>();
 
  t.run_test_cases();
  return 0;
}
