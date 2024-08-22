#include "../../environment/unit_test.hpp"
#include "spherical_harmonics_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric/spherical_harmonics");
  t.include<spherical_harmonics_test>();
  t.run_test_cases();
  return 0;
} 
