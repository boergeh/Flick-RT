#include "../environment/unit_test.hpp"
#include "single_layer_slab_test.hpp"

int main() {
  using namespace flick;
  unit_test t("model");
  t.include<single_layer_slab_test_A>("single_layer_slab_test_A");
  t.include<single_layer_slab_test_B>("single_layer_slab_test_B");
  t.include<single_layer_slab_test_C>("single_layer_slab_test_C");
  t.include<single_layer_slab_test_D>("single_layer_slab_test_D");
  t.include<single_layer_slab_test_E>("single_layer_slab_test_E");
  t.include<single_layer_slab_test_F>("single_layer_slab_test_F");
  t.run_test_cases();
  return 0;
}
