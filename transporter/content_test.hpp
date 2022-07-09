#include "content.hpp"
//#include "../geometry/volume.hpp" // should not do this
namespace flick {
  begin_test_case(content_test) {
    //using sphere = geometry::sphere<content>;
    //sphere sp{1};
    //std::cout << sp();
    emitter em({0,0,0},stokes{1,0,0,0},1e3);
    //sp.content().insert(em);
  } end_test_case()
}
