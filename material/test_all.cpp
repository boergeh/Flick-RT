#include "../environment/unit_test.hpp"
#include "material_test.hpp"
#include "iop_profile_test.hpp"
#include "spheres_test.hpp"
#include "normalized_scattering_matrix_fit_test.hpp"
#include "ab_functions_test.hpp"
#include "layered_iops_test.hpp"
#include "z_profile_test.hpp"

int main() {
  using namespace flick;
  unit_test t("material");
  t.include<material_test_A>();
  t.include<iop_profile_test>();
  t.include<spheres_test_A>();
  t.include<spheres_test_B>();
  t.include<normalized_scattering_matrix_fit_test>();
  t.include<ab_functions_test_A>();
  t.include<ab_functions_test_B>();
  t.include<layered_iops_test>();
  t.include<z_profile_aggregate_test>();
  //t.include<accurt_user_specified_test>();
  t.run_test_cases();
  return 0;
}
