#ifndef flick_command_help
#define flick_command_help

#include "basic_command.hpp"

namespace flick {
  namespace command {
    struct help : public basic_command {
      help():basic_command("help"){};
      void run() {
	if (a(1).empty())
	  show("help.txt");
	else if (!a(1).empty()) {
	  try {
	    show(a(1)+".txt");
	  } catch (const std::exception& e) {
	    std::cout << "Cannot recognize '" + a(1) << "'." << std::endl;
	  }
	} else {
	  error();
	}
      }
    };
  }
}
       
#endif
