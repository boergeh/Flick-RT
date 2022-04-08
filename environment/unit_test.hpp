#ifndef flick_unit_test
#define flick_unit_test

#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>
#include <vector>
#include <chrono>
#include <iomanip>

namespace flick {
  std::chrono::time_point<std::chrono::system_clock> time_1, time_2;
  void start_time() {
     time_1 = std::chrono::system_clock::now();
  }
  void show_time() {
    time_2 = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = time_2 - time_1;
    std::cout <<" - "<< std::setprecision(2) << duration.count() << " s\n";
  }
  
  class test_case {
    void write_begin() {
      std::cout << std::endl << " " << name_ << ": ";
    }
    void write_end() {
      std::cout << std::endl;
    }
  protected:
    std::string name_;
    bool do_printing_{true};
    int errors_{0};
    virtual ~test_case(){}
  public:
    test_case(const std::string& name) : name_{name}{}
    virtual void test() = 0;
    int errors(){return errors_;}
    void do_printing(bool b) {do_printing_ = b;}
    template<class T>
    void print(std::string s, T t) {
      if (do_printing_)
	std::cout << s << t << ",";
    }
    void print_progress() {
      if (do_printing_)
	std::cout << "[" << name_ << "]";
    }
    void check(bool b, const std::string& s="") {
      if (!b) {
	write_begin();
	std::cout << "\""<< s << "\" boolean test failed ";
	write_end();
	errors_++;
      }
    }
    void check_close(double value1, double value2, double percent, const std::string& s="") {
      double diff = fabs(value1/value2-1)*100;
      if (!std::isfinite(diff) | (diff > percent)) {
	write_begin();
	std::cout <<  "\""<< s << "\" " << " the difference between "
		  << value1 << " and " << value2 << " is " << diff
		  << " %, which is larger than the accepted "
		  << percent << " %";
	write_end();
	errors_++;
      }
    }
    void check_small(double value, double accepted_distance, const std::string& s="") {
      if (!std::isfinite(value) | (fabs(value) > accepted_distance)) {
	write_begin();
	std::cout << "\""<< s << "\" " << value
		  << " is further from zero than the accepted distance "
		  << accepted_distance;
	write_end();
	errors_++;
      }
    }
    void report_throw_failure() {
      write_begin();
      std::cout << "did not throw as expected";
      write_end();
    }
  };

  class unit_test {
    std::string name_;
    std::vector<std::shared_ptr<test_case>> test_cases_;
  public:
    unit_test(const std::string& name) : name_{name} {}
    template<class T>
    void include(const std::string& name) {
      test_cases_.push_back(std::make_shared<T>(name));
    }
    void run_test_cases() {
      start_time();
      std::string s;
      if (test_cases_.size() > 1)
	s = "s";
      std::cout << "Running " << test_cases_.size() << " test case"
		<< s << " in " << name_ << " unit test:";
      int errors = 0;
      for (auto& current_case: test_cases_) {
	current_case->print_progress();
	current_case->test();
	errors += current_case->errors();
      }
      std::cout << std::endl;
      if (errors <= 1) s = ""; else s = "s";
      std::stringstream ss;
      ss << "Found ";
       if (errors > 0)
	ss << errors;
      else
	ss << "no";
      ss << " error" << s << " in "
	 << name_ << " unit test";

      std::string color_end = "\033[0m";
      std::string red_start = "\033[1;31m";
      std::string green_start = "\033[1;32m";
      if (errors > 0)
	std::cout << red_start << ss.str() << color_end;
      else
	std::cout << green_start << ss.str() << color_end;
      show_time();
    }
  };
}

#endif

#define check_throw(expression) {			\
    try {						\
      expression;					\
      report_throw_failure();				\
    } catch(...) {}					\
  } 

#define begin_test_case(testname)  			\
  class testname : public test_case {			\
    using test_case::test_case;				\
    void test() {					\
      
#define end_test_case()				\
  }};						
