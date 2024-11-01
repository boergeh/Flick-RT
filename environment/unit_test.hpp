#ifndef flick_unit_test
#define flick_unit_test

#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>
#include <vector>
#include <chrono>
#include <iomanip>
//#include <format>

namespace flick {
  constexpr long double operator""_pct(long double p)
  {
    return p;
  }
  constexpr long double operator""_pct(unsigned long long p)
  {
    return p;
  }
  std::chrono::time_point<std::chrono::system_clock> time_1, time_2;
  void start_time_1() {
     time_1 = std::chrono::system_clock::now();
  }
  void start_time_2() {
     time_2 = std::chrono::system_clock::now();
  }
  void show_time_1() {
    auto time = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = time - time_1;
    //std::cout <<" - "<< std::format("{:.2}",duration.count()) << " s\n";
    std::cout <<" - "<< std::setprecision(2) << duration.count() << " s\n";
    std::cout << std::setprecision(6);
  }
  void show_time_2() {
    auto time = std::chrono::system_clock::now();
    std::chrono::duration<double> duration = time - time_2;
    //std::cout << std::format("{:.1}",1000*duration.count()) << " ms";
    std::cout << std::setprecision(1) << 1000*duration.count() << " ms";
    std::cout << std::setprecision(6);
  }
  double cpu_duration() {
    auto t0 = std::chrono::system_clock::now();
    for (size_t i=0; i<10000; i++)
      ;
    auto d = std::chrono::system_clock::now() - t0; 
    return 1e2*std::chrono::duration<double>(d).count();
  }
  
  class test_case {
  protected:
    std::string name_;
    bool do_printing_{true};
    int errors_{0};
    int checks_{0};
    virtual ~test_case() = default;
  public:
    void write_begin(const std::string& s) {
      std::cout << "\n\n" << " " << name_ << ", check number "
		<< std::to_string(checks_) <<":";
      if (not s.empty())
	std::cout << " \"" << s << "\".";
      std::cout << std::setprecision(3) << " ";
    }
    void write_end() {
      std::cout << std::setprecision(6) << std::endl << std::endl;
    }
    test_case(const std::string& name="") : name_{name}{}
    virtual void test() = 0;
    int errors() {return errors_;}
    void do_printing(bool b) {do_printing_ = b;}
    template<class T>
    void print(std::string s, T t) {
      if (do_printing_)
	std::cout << s << t << "," << std::flush;
    }
    void print_progress_begin() {
      if (do_printing_)
	std::cout << "[" << name_ << ", "<< std::flush;
    }
    void print_progress_end() {
      if (do_printing_) {
	std::cout << checks_ << " checks, ";
	show_time_2();
	std::cout << "]" << std::flush;
      }
    }
    void check(bool b, const std::string& s="") {
      checks_++;
      if (!b) {
	write_begin(s);
	std::cout << "Boolean test failed ";
	write_end();
	errors_++;
      }
    }
    void check_close(double value1, double value2, double percent=1e-12,
		     const std::string& s="") {
      checks_++;
      double diff = (value1/value2-1)*100;
      if (!std::isfinite(diff) | (fabs(diff) > percent)) {
	write_begin(s);
	std::cout << "The difference between "
		  << value1 << " and benchmark " << value2 << " is " << diff
		  << " %, which is larger than the accepted "
		  << percent << " %";
	write_end();
	errors_++;
      }
    }
    
    void check_small(double value, double accepted_distance=1e-14,
		     const std::string& s="") {
      checks_++;
      if (!std::isfinite(value) | (fabs(value) > accepted_distance)) {
	write_begin(s);
	std::cout << value
		  << " is further from zero than the accepted distance "
		  << accepted_distance;
	write_end();
	errors_++;
      }
    }
    void check_fast(double max_time, const std::string& s="") {
      checks_++;
      auto time = std::chrono::system_clock::now();
      std::chrono::duration<double> duration = time - time_2;
      if (duration > std::chrono::duration<double>(max_time)) {
	write_begin(s);
	std::cout << "Runtime of "
	  //<< std::format("{:.2}", duration.count())
		  << std::setprecision(2) << duration.count()
		  << " s is slower than the accepted " << max_time
		  << " s";
	std::cout << std::setprecision(6);
	write_end();
	errors_++;
      }
    }
    void report_throw_failure() {
      write_begin("");
      std::cout << "did not throw as expected";
      write_end();
      errors_++;
    }
  };

  class unit_test {
    std::string name_;
    std::vector<std::shared_ptr<test_case>> test_cases_;
  public:
    unit_test(const std::string& name="") : name_{name} {}
    template<class T>
    void include(std::string name="") {
      if (name.empty()) {
	name = T().class_name;
      }
      test_cases_.push_back(std::make_shared<T>(name));
    }
    void run_test_cases() {
      start_time_1();
      std::string s;
      if (test_cases_.size() > 1)
	s = "s";
      std::cout << "Running " << test_cases_.size() << " test case"
		<< s << " in " << name_ << " unit test: ";
      int errors = 0;
      for (auto& current_case: test_cases_) {
	start_time_2();
	current_case->print_progress_begin();
	try {
	  current_case->test();
	} catch (const std::exception& e) {
	  std::stringstream ss;
	  ss << "Standard exception: " << e.what();
	  current_case->write_begin(ss.str());
	  current_case->write_end();
	  errors++;
	}
	current_case->print_progress_end();
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
      ss << " error" << s << " in " << name_ << " unit test";
      std::string color_end = "\033[0m";
      std::string red_start = "\033[1;31m";
      std::string green_start = "\033[1;32m";
      if (errors > 0)
	std::cout << red_start << ss.str() << color_end;
      else
	std::cout << green_start << ss.str() << color_end;
      show_time_1();
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
  struct testname : public test_case {			\
    using test_case::test_case;				\
    const std::string class_name=#testname;	        \
    void test() {					\
    
#define end_test_case()		            		\
  }};						
