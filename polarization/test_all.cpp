#include "../environment/unit_test.hpp"
#include "algorithm_test.hpp"
#include "fresnel_test.hpp"

int main() {
  using namespace flick;
  unit_test t("polarization");
  t.include<algorithm_test>("algorithm_test");
  t.include<fresnel_test_A>("fresnel_test_A");
  t.include<fresnel_test_B>("fresnel_test_B");
  t.run_test_cases();
  return 0;
}
