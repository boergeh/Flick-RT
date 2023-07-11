#include "../environment/unit_test.hpp"
#include "single_layer_slab_test.hpp"
#include "multilayer_test.hpp"

int main() {
  using namespace flick;
  unit_test t("model");
  t.include<multilayer_test>("multilayer");    
  t.include<single_layer_slab_test_A>();  
  t.include<single_layer_slab_test_B>();
  t.include<single_layer_slab_test_C>();
  t.include<single_layer_slab_test_D>();
  t.include<single_layer_slab_test_E>();
  t.include<single_layer_slab_test_F>();    
  t.include<single_layer_slab_test_G>(); 
  t.include<single_layer_slab_test_H>(); 
  t.include<single_layer_slab_test_I>(); 
  t.run_test_cases();
  return 0;
}
