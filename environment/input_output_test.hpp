#include "input_output.hpp"

namespace flick {
  begin_test_case(input_output_test) {
    std::string file_name = "flick_tmp.txt";
    write("test_text",file_name);
    check(read<std::string>(file_name)=="test_text");
    check(std::filesystem::remove(file_name));

    file_name = path()+"/flick_tmp.txt";
    write("test_text",file_name);
    check(read<std::string>(file_name)=="test_text");
    check(std::filesystem::remove(file_name));

    file_name = "../flick_tmp.txt";
    write("test_text",file_name);
    check(read<std::string>("flick_tmp.txt")=="test_text");
    check(std::filesystem::remove(file_name));
  } end_test_case()
}
