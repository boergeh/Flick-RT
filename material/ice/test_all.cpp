#include "../../environment/unit_test.hpp"
#include "pure_ice_test.hpp"

int main() {
  using namespace flick;
  unit_test t("ice");
  t.include<pure_ice_test>("pure_ice_test");
  t.run_test_cases();
  return 0;
}
