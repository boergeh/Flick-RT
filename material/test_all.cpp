#include "../environment/unit_test.hpp"
#include "material_test.hpp"

int main() {
  using namespace flick;
  unit_test t("material");
  t.include<material_test>("material_test");
  t.run_test_cases();
  return 0;
}
