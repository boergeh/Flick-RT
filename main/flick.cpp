#include <string>
#include <vector>

#include "commands/help.hpp"
#include "commands/radiator.hpp"
#include "commands/iop.hpp"
#include "commands/text.hpp"

namespace flick {
  template<typename Command>
  bool run(const std::vector<std::string>& args) {
    Command c;
    c.set_arguments(args);
    if (!args.empty() && c.has_name(args.at(0))) {
      try {
	c.run();
      } catch (const flick::exception& e) {
	std::cout << "\n Flick exception in " << e.what() << "\n\n";
      } catch (const std::exception& e) {
	std::cout << "\n c++ standard exception: " << e.what() << "\n\n";
      } catch (...) {
	std::cout << "\n Unknown exception" << "\n\n";
      }
      return true;
    }
    return false;
  }
}

int main(int argc, char* argv[]) {
  using namespace flick;
  const std::vector<std::string> a(argv + 1, argv + argc);
  if (a.empty()) {run<command::help>({"help"}); return 0;}
  if (run<command::help>(a)) return 0;
  if (run<command::radiator>(a)) return 0;
  if (run<command::iop>(a)) return 0;
  if (run<command::text>(a)) return 0;
  std::cout << "\n Cannot recognize command. Try flick help.\n" << std::endl;
  return 0;
}

