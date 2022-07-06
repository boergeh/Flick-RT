#include "../environment/unit_test.hpp"
#include "algorithm_test.hpp"

int main() {
  using namespace flick;
  unit_test t("polarization");
  t.include<algorithm_test>("algorithm_test");
  t.run_test_cases();
  return 0;
}
