#include "../../environment/unit_test.hpp"
#include "lines_test.hpp"
#include "profile_test.hpp"

int main() {
  using namespace flick;
  auto start = std::chrono::steady_clock::now();
  unit_test t("gases");
  t.include<profile_test>("profile_test");
  t.include<lines_test>("lines_test");
  t.run_test_cases();
  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> duration = end-start;
  std::cout << std::setprecision(2) << duration.count() << "s\n";
  return 0;
}
