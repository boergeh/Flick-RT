#include "../environment/unit_test.hpp"
#include "ordinary_mc_test.hpp"

int main() {
  using namespace flick;
  unit_test t("transporter");
  /*
  t.include<ordinary_mc_test_A>("ordinary_mc_test_A");
  t.include<ordinary_mc_test_B>("ordinary_mc_test_B");
  t.include<ordinary_mc_test_C>("ordinary_mc_test_C");
  t.include<ordinary_mc_test_D>("ordinary_mc_test_D");
  */
  t.include<ordinary_mc_test_E>("ordinary_mc_test_E");
  t.run_test_cases();
  return 0;
}
