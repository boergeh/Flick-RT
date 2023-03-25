#include "../environment/unit_test.hpp"
#include "mie_test.hpp"

int main() {
  using namespace flick;
  unit_test t("mie");
  t.include<bessel_test>("bessel_test");
  t.include<mie_test_A>("mie_test_A");
  t.include<mie_test_B>("mie_test_B");
  t.include<mie_test_C>("mie_test_C");
  t.include<mie_test_D>("mie_test_D");
  t.include<mie_test_E>("mie_test_E");
  t.run_test_cases();
  return 0;
} 
