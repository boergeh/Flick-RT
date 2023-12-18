#include "search_and_replace.hpp"

namespace flick {
  begin_test_case(search_and_replace_test) {
    std::stringstream ss("# test ## \n list = 0 1 2 \n # ## \n cdom_440 = 1");
    parameter_text t;
    t.set_begin_qualifier("#");
    ss >> t;
    std::string l = t.get("list");
    check(l=="0 1 2");
    t.set("list","  0 1 2 3 4#");
    l = t.get("list");
    check(l=="0 1 2 3 4");

    l = t.get("cdom_440");
    check(l=="1");
    t.set("cdom_440","  0.2");
    l = t.get("cdom_440");
    check(l=="0.2");

  } end_test_case()
}
