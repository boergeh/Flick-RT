#include <iostream>
#include <fstream>
#include <string>
#include "../numeric/function.hpp"
#include "../environment/input_output.hpp"

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
void show(const std::string& fname) {
  auto t = read<text>("main/help/"+fname);
  std::cout << t << std::endl;  
}

void help(const std::string& a) {
  if (a=="help")
    show("help.txt");
  else if (a=="list")
    show("list.txt");
}

void list(const std::vector<std::string>& a) {
  if (a.at(0)=="--precision")
    set_precision(a.at(1),std::stoi(a.at(2)),std::stoi(a.at(3)));
  else if (a.at(0)=="--xy")
    list_xy(a.at(1));
}

int main(int argc, char* argv[]) {
  const std::vector<std::string> a(argv + 1, argv + argc);
  if (a.size()==0 || (a.size()==1 && a.at(0)=="help"))
    show("help.txt");
  else if (a.at(0)=="help" && a.size()==2)
    help(a.at(1));
  else if (a.at(0)=="list")
    list({a.begin()+1, a.end()});
  else if (a.at(0)=="list-xy")
    list_xy(a.at(1));
  else
    std::cout << "\nCannot recognize command. Try flick help\n" << std::endl;  
  return 0;
}
