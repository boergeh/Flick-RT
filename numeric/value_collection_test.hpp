#include "value_collection.hpp"

namespace flick {
  begin_test_case(value_collection_test) {
    value_collection vc(0.01);
    vc.add(2,1);
    vc.add(2,1);
    vc.add(1,1);
    vc.add(1,1);
    check_close(vc.mean(),1.5);
    check_close(vc.std(),0.5);
    check(not vc.accurate());
  } end_test_case()
}
