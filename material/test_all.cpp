#include "../environment/unit_test.hpp"
#include "material_test.hpp"

int main() {
  using namespace flick;
  unit_test t("material");
  auto start = std::chrono::steady_clock::now();
  t.include<material_test>("material_test");
  t.run_test_cases();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = end-start;
  std::cout << std::setprecision(2) << duration.count() << "s\n";
  return 0;
}
