#include "content.hpp"
#include "../geometry/volume.hpp"
namespace flick {
  begin_test_case(content_test) {
    using sphere = geometry::sphere<content>;
    sphere sp{1};
    //std::cout << sp();
    emitter em{1000,stokes{1,0,0,0}};
    sp.content().insert(em);
    
  } end_test_case()
}
