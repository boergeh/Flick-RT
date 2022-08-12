#include "../environment/unit_test.hpp"
#include "boundary_test.hpp"
#include "volume_test.hpp"

int main() {
  using namespace flick;
  using namespace flick::geometry;
  unit_test t("geometry");
  t.include<boundary_test>("boundary_test");
  t.include<volume_test_A>("volume_test_A");
  t.include<volume_test_B>("volume_test_B");
  t.run_test_cases();
  return 0;
}
