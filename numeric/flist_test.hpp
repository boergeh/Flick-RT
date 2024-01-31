#include "flist.hpp"

namespace flick {
  begin_test_case(flist_test) {
    std::stringstream ss("/* Header text */\n\n 0 0.1 0.2 \n 1 1.1 1.2");
    pl_flist fl;
    ss >> fl;
    check_close(fl(0).value(0),0.1);
    check_close(fl(1).value(1),1.2);
  } end_test_case()
}
