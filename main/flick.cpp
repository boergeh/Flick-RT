#include <string>
#include <vector>

#include "commands/basic_command.hpp"
#include "commands/help.hpp"
#include "commands/radiator.hpp"
#include "commands/filter.hpp"
#include "commands/mie.hpp"
#include "commands/iop.hpp"
#include "commands/delta_fit.hpp"
#include "commands/text.hpp"

int main(int argc, char* argv[]) {
  using namespace flick;
  const std::vector<std::string> a(argv + 1, argv + argc);
  if (a.empty()) {run<command::help>({"help"}); return 0;}
  if (run<command::help>(a)) return 0;
  if (run<command::radiator>(a)) return 0;
  if (run<command::filter>(a)) return 0;
  if (run<command::mie>(a)) return 0;
  if (run<command::iop>(a)) return 0;
  if (run<command::delta_fit>(a)) return 0;
  if (run<command::text>(a)) return 0;
  std::cout << "\n Cannot recognize command, try 'flick help'.\n" << std::endl;
  return -1;
}

