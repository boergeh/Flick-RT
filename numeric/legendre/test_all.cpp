#include "../../environment/unit_test.hpp"
#include "legendre_test.hpp"
#include "delta_fit_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric/legendre");
  t.include<legendre_test_A>();
  t.include<legendre_test_B>();
  t.include<delta_fit_test>();
  t.include<legendre_test_accumulated_integral>();
  t.run_test_cases();
  return 0;
} 
