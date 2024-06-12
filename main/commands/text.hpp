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
	  std::string parameter = a(3);
	  std::string value = a(4);
	  auto t = read<parameter_text>(fname);
	  t.set(parameter, value);
	  std::cout << t;
	}
	else if (command=="xy-precision")  {
	  size_t nx = std::stoi(a(3));
	  size_t ny = std::stoi(a(4));
	  auto f = read<pl_function>(fname);
	  std::cout << std::setprecision(std::max(nx,ny)) << std::showpoint <<
	    significant_digits(f,nx,ny);
	  std::cout << std::setprecision(6);
	}
	else if (command=="xy-scale")  {
	  if (size() > 5)
	    std::cout << std::showpoint << std::setprecision(std::stoi(a(5)));
	  double kx = std::stod(a(3));
	  double ky = std::stod(a(4));
	  auto f = read<pl_function>(fname);
	  auto x = f.x();
	  auto y = f.y();
	  for (size_t i=0; i < x.size(); ++i)
	    std::cout << kx*x[i] << " " << ky*y[i] << '\n';
	}
	else if (command=="xy")  {
	  if (size() > 3)
	    std::cout << std::showpoint << std::setprecision(std::stoi(a(3)));
	  auto f = read<pl_function>(fname);
	  auto x = f.x();
	  auto y = f.y();
	  for (size_t i=0; i < x.size(); ++i)
	    std::cout << x[i] << " " << y[i] << '\n';
	}
	else if (command=="matrix")  {
	  using namespace linalg;
	  auto f = read<pl_flist>(fname);
	  if (size() > 3) {
	    std::cout << f.header();
	    std::cout << std::showpoint << std::setprecision(std::stoi(a(3)));
	  }
	  std::cout << f.matrix();
	}
	else {
	  error();
	}
      }
    };
  }
}
