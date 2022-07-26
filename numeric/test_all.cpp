#include "../environment/unit_test.hpp"
#include "direction_generator_test.hpp"
#include "vector_test.hpp"
#include "histogram_test.hpp"
#include "pose_test.hpp"
#include "range_test.hpp"
#include "sorted_vector_test.hpp"
#include "function_test.hpp"
#include "physics_function_test.hpp"

int main() {
  using namespace flick;
  unit_test t("numeric");
  t.include<sorted_vector_test>("sorted_vector_test");
  t.include<function_test>("function_test");
  t.include<direction_generator_test>("direction_generator_test");
  t.include<vector_test>("vector_test");
  t.include<histogram_test>("histogram_test");
  t.include<pose_test>("pose_test");
  t.include<range_test>("range_test");
  t.include<physics_function_test>("physics_function_test");
  t.run_test_cases();
  return 0;
} 
