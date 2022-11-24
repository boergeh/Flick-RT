#include "basic_command.hpp"
#include "../../environment/input_output.hpp"
#include "../../numeric/function.hpp"

namespace flick {
  namespace command {
    struct text : public basic_command {
      text():basic_command("text"){};
      void run() {
	if (a(1)=="xy") {
	  std::string fname = a(2);
	  auto f = read<pl_function>("./"+fname);
	  auto x = f.x();
	  auto y = f.y();
	  for (size_t i=0; i < x.size(); ++i)
	    std::cout << x[i] << " " << y[i] << '\n';
	} else if (a(1)=="xy-precision")  {
	  size_t nx = std::stoi(a(2));
	  size_t ny = std::stoi(a(3));
	  std::string fname = a(4);
	  auto f = read<pp_function>("./"+fname);
	  std::cout << significant_digits(f,nx,ny);
	} else if (a(1)=="xy-scale")  {
	  double kx = std::stod(a(2));
	  double ky = std::stod(a(3));
	  std::string fname = a(4);
	  auto f = read<pp_function>("./"+fname);
	  auto x = f.x();
	  auto y = f.y();
	  for (size_t i=0; i < x.size(); ++i)
	    std::cout << kx*x[i] << " " << ky*y[i] << '\n';
	} else {
	  error();
	}
      }
    };
  }
}
