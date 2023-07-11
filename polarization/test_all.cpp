#include "../environment/unit_test.hpp"
#include "algorithm_test.hpp"
#include "fresnel_test.hpp"
#include "stokes_test.hpp"
#include "mueller_test.hpp"

int main() {
  using namespace flick;
  unit_test t("polarization");
  t.include<stokes_test>();
  t.include<algorithm_test>();
  t.include<fresnel_test_A>();
  t.include<fresnel_test_B>();
  t.include<mueller_test>();
  t.run_test_cases();
  return 0;
}
