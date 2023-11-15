#include "search_and_replace.hpp"

namespace flick {
  begin_test_case(search_and_replace_test) {
    std::stringstream ss("# test ## \n list = 0 1 2");
    parameter_text t;
    ss >> t;
    std::string l = t.get("list");
    check(l=="0 1 2");
    t.set("list","  0 1 2 3 4#");
    l = t.get("list");
    check(l=="0 1 2 3 4");
  } end_test_case()
}
