#include "../environment/unit_test.hpp"
#include "boundary_test.hpp"
#include "volume_test.hpp"

int main() {
  using namespace flick;
  using namespace flick::geometry;
   auto start = std::chrono::steady_clock::now();
  unit_test t("geometry");
  t.include<boundary_test>("boundary_test");
  t.include<volume_test>("volume_test");
  t.run_test_cases();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = end-start;
  std::cout << std::setprecision(2) << duration.count() << "s\n";
  return 0;
}
