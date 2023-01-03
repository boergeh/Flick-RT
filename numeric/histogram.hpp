#ifndef flick_histogram
#define flick_histogram

#include <iostream>
#include <vector>
#include <iomanip>
#include <fstream>

namespace flick {
  class equal_bins {
  public:
    double min;
    double max;
    size_t n;
    double midpoint(size_t i) const {
      return min + (max-min)/n*(i+0.5);
    }
  };

  class histogram {
    // Histogram with bin ranges [min max>. Adding value at max will
    // throw out-of-range error.
  private:
    equal_bins x_bins;
    equal_bins y_bins;
    std::vector<std::vector<double>> h;
  public:
    histogram()
      : x_bins{0,1,1}, y_bins{0,1,1}, h{1,std::vector<double>(1,0.0)} {}
    histogram(equal_bins xb, equal_bins yb)
      : x_bins{xb}, y_bins{yb}, h{xb.n,std::vector<double>(yb.n,0.0)} {} 
    histogram(equal_bins xb)
      : x_bins{xb}, y_bins{equal_bins{-1,1,1}},
	h{xb.n,std::vector<double>(1,0.0)} {} 
    void add(double x, double y, double value) {
      size_t i = x_bins.n*(x-x_bins.min)/(x_bins.max-x_bins.min);
      size_t j = y_bins.n*(y-y_bins.min)/(y_bins.max-y_bins.min);
      h.at(i).at(j) += value;
    }
    void add(double x, double value) {
      add(x,0,value);
    }
    void add(const histogram& h_to_be_added) {
      for (size_t i = 0; i < x_bins.n; ++i)
	for (size_t j = 0; j < y_bins.n; ++j)
	  h.at(i).at(j) += h_to_be_added.bin_value(i,j);
    }
    double bin_value(size_t i, size_t j) const {
      return h.at(i).at(j);
    } 
    void write(std::ostream& os) const {
      if (y_bins.n > 1) {
	os << "nan" << " ";
	for (size_t j = 0; j < y_bins.n; ++j) {
	  os <<  y_bins.midpoint(j) << " ";
	}
	os << '\n';
      }
      for (size_t i = 0; i < x_bins.n; ++i) {
	os << x_bins.midpoint(i) << " ";
	for (size_t j = 0; j < y_bins.n; ++j) {
	  os << h[i][j] << " ";
	}
	os << '\n';
      }
    }
    friend std::ostream& operator<<(std::ostream &os, const histogram& h) {
      h.write(os);
      return os;
    }
  };
  void write(const histogram& h, const std::string& file_name,
	     size_t precision=4) {
    std::ofstream of(file_name);
    of << std::setprecision(precision);
    h.write(of);
    of.close();
  }
}

#endif

