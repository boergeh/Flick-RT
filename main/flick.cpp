#include <string>
#include <vector>

#include "commands/basic_command.hpp"
#include "commands/help.hpp"
#include "commands/radiator.hpp"
#include "commands/time.hpp"
#include "commands/sun_position.hpp"
#include "commands/filter.hpp"
#include "commands/mie.hpp"
#include "commands/iop.hpp"
#include "commands/delta_fit.hpp"
#include "commands/text.hpp"
#include "commands/accurt.hpp"

int main(int argc, char* argv[]) {
  using namespace flick;
  try {
    const std::vector<std::string> a(argv + 1, argv + argc);  
    if (a.empty()) {run<command::help>({"help"}); return 0;}
    if (run<command::help>(a)) return 0;
    if (run<command::radiator>(a)) return 0;
    if (run<command::time>(a)) return 0;
    if (run<command::sun_position>(a)) return 0;
    if (run<command::filter>(a)) return 0;
    if (run<command::mie>(a)) return 0;
    if (run<command::iop>(a)) return 0;
    if (run<command::delta_fit>(a)) return 0;
    if (run<command::text>(a)) return 0;
    if (run<command::accurt>(a)) return 0;
    throw std::runtime_error("cannot recognize command, try\n\n flick help\n\n");
  }
  catch (const flick::exception& e) {
    std::cerr << "flick exception: " << e.what() << "\n";
    return -1;
  }
  catch(std::exception& e) {
    std::cerr << "flick std exception: " << e.what() << "\n";
    return -1;
  }
  catch(...) {
    std::cerr << "flick unknown exception \n";
    return -1;
  }
  return 0;
}

