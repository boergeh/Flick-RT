#include <iostream>
#include <fstream>
#include <string>
#include "../numeric/function.hpp"
#include "../environment/input_output.hpp"
#include "../environment/input_output.hpp"
#include "../material/water/pure_water.hpp"
#include "../numeric/range.hpp"

using namespace flick;

void set_precision(const std::string& fname, size_t nx, size_t ny)
{
  auto f = read<pp_function>(fname);
  std::cout << significant_digits(f,nx,ny);
}
void list_xy(std::string fname)
{
  auto f = read<pl_function>("./"+fname);
  auto x = f.x();
  auto y = f.y();
  for (size_t i=0; i < x.size(); ++i)
    std::cout << x[i] << " " << y[i] << '\n';
}
/*
void help(const std::string& a) {
  if (a=="help")
    show("help.txt");
  else if (a=="list")
    show("list.txt");
}
*/
void list(const std::vector<std::string>& a) {
  if (a.at(0)=="--precision")
    set_precision(a.at(1),std::stoi(a.at(2)),std::stoi(a.at(3)));
  else if (a.at(0)=="--xy")
    list_xy(a.at(1));
}

class basic_run {
  std::vector<std::string> arguments_;
  std::string name_;
public:
  basic_run(std::string name) : name_{name}{}
  std::string a(size_t n) {
    if (n >= arguments_.size())
      return "";
    return arguments_.at(n);
  }
  void set_arguments(const std::vector<std::string>& args) {
    arguments_ = args;
  }
  bool has_name(const std::string& name) {
    return (name_ == name);
  }
  void error() {
    std::cout << "\nCannot recognize "+name_+
    " arguments. Try 'flick help "+name_+"'.\n" << std::endl;  
  }
  void show(const std::string& fname) {
    auto t = read<text>("main/help/"+fname);
    std::cout << t << std::endl;  
  }

};

struct iop : public basic_run {
  iop():basic_run("iop"){};
  void run() {
    if (a(1)=="pure_water") {
      auto wls = range(std::stod(a(3)),std::stod(a(4)),std::stoi(a(5))).logspace();
      pure_water pw;
      if (a(2)=="absorption_length") {
	for (auto wl:wls) {
	  pw.wavelength(wl);
	  std::cout << std::setprecision(4) << wl << " "
		    << 1/pw.absorption_coefficient() << '\n';
	}
      }
    } else {
      error();
    }
  }
};

struct help : public basic_run {
  help():basic_run("help"){};
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

template<typename Command>
bool run(const std::vector<std::string>& args) {
  Command c;
  c.set_arguments(args);
  if (c.has_name(args.at(0))) {
    c.run();
    return true;
  }
  return false;
}

int main(int argc, char* argv[]) {
  const std::vector<std::string> a(argv + 1, argv + argc);
  if (a.empty()) {run<help>({"help"}); return 0;}
  if (run<help>(a)) return 0;
  if (run<iop>(a)) return 0;
  std::cout << "\nCannot recognize command. Try 'flick help'.\n" << std::endl;
  /*
  if (a.size()==0 || (a.size()==1 && a.at(0)=="help"))
    show("help.txt");
  else if (a.at(0)=="help" && a.size()==2)
    help(a.at(1));
  else if (a.at(0)=="list")
    list({a.begin()+1, a.end()});
  else (a.at(0)=="list-xy")
    list_xy(a.at(1));
  */
  return 0;
}
