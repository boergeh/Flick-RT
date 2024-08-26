#include "../../environment/unit_test.hpp"
#include "spherical_harmonics_test.hpp"
#include "sh_expansion_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric/spherical_harmonics");
  t.include<spherical_harmonics_test>();
  t.include<sh_expansion_test>();
  t.run_test_cases();
  return 0;
} 
