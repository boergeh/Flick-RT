#include "basic_command.hpp"
#include "../../environment/input_output.hpp"
#include "../../environment/search_and_replace.hpp"
#include "../../numeric/function.hpp"

namespace flick {
  namespace command {
    struct text : public basic_command {
      text():basic_command("text"){};
      void run() {
	auto fname = a(1);
	auto command = a(2);
	if (command=="set") {
	  auto parameter = a(3);
	  auto value = a(4);
	  auto t = read<parameter_text>("./"+fname);
	  t.set(parameter, value);
	  std::cout << t;
	}
	else if (command=="xy-precision")  {
	  size_t nx = std::stoi(a(3));
	  size_t ny = std::stoi(a(4));
	  auto f = read<pl_function>("./"+fname);
	  std::cout << significant_digits(f,nx,ny);
	}
	else if (command=="xy-scale")  {
	  double kx = std::stod(a(3));
	  double ky = std::stod(a(4));
	  auto f = read<pl_function>("./"+fname);
	  auto x = f.x();
	  auto y = f.y();
	  for (size_t i=0; i < x.size(); ++i)
	    std::cout << kx*x[i] << " " << ky*y[i] << '\n';
	}
	else {
	  error();
	}
      }
    };
  }
}
