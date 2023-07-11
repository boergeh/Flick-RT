#include "../../environment/unit_test.hpp"
#include "matrix_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric/linalg");
  t.include<matrix_test>();
  t.run_test_cases();
  return 0;
} 
