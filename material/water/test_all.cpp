#include "../../environment/unit_test.hpp"
#include "pure_water_test.hpp"
#include "phytoplankton_test.hpp"
#include "nap_test.hpp"

int main() {
  using namespace flick;
  unit_test t("water");
  t.include<pure_water_test>();
  t.include<phytoplankton_test>();
  t.include<nap_test>();
 
  t.run_test_cases();
  return 0;
}
