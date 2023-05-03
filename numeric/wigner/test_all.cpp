#include "../../environment/unit_test.hpp"
#include "wigner_fit_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric/wigner");
  t.include<wigner_fit_test_A>("wigner_fit_test_A");
  t.include<wigner_fit_test_B>("wigner_fit_test_B");
  t.run_test_cases();
  return 0;
} 
