#include "input_output.hpp"

namespace flick {
  begin_test_case(input_output_test) {
    std::string file_name = "./tmp.txt";
    write("test_text",file_name);
    check(read<std::string>(file_name)=="test_text");
    check(std::filesystem::remove(file_name));
  } end_test_case()
}
