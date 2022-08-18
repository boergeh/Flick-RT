#include "../environment/unit_test.hpp"
#include "single_layer_slab_test.hpp"

int main() {
  using namespace flick;
  unit_test t("model");
  t.include<single_layer_slab_test>("single_layer_slab_test");
  t.run_test_cases();
  return 0;
}
