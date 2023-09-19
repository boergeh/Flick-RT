#ifndef flick_table
#define flick_table
 
#include "function.hpp"

namespace flick {
  template<class Interpolation>
  class table
  {
    std::string header_;
    std::vector<double> row_values_;
    std::vector<double> col_values_;
    std::vector<function<Interpolation>> column_functions_;
  public:
    table() = default;
    std::string header() const {
      return header_;
    }
    function<Interpolation> row(size_t n) {
      function<Interpolation> f;
      for (size_t i=0; i < row_values_.size(); ++i)
	f.append({row_values_[i],column_functions_[i].y().at(n)});
      return f;
    }
    function<Interpolation> column(size_t n) {
      return function<Interpolation>{column_functions_.at(n).x(),
	column_functions_.at(n).y()};
    }
    double value(double row_value, double column_value) const {
      function<Interpolation> f;
      for (size_t i=0; i < col_values_.size(); ++i)
	f.append({col_values_[i], column_functions_[i].value(row_value)});
      return f.value(column_value);
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const table<Interpolation>& f) {
      os << f.header();
      os << f.row_values_.size() << std::endl;
      os << f.col_values_.size() << std::endl;
      for (size_t i=0; i<f.row_values_.size(); ++i) {
	os << f.row_values_[i] << " ";
      }
      os << std::endl;
      for (size_t i=0; i<f.col_values_.size(); ++i) {
	os << f.col_values_[i] << " ";
      }
      os << std::endl;
      for (size_t i=0; i<f.row_values_.size(); ++i) {
	for (size_t j=0; j<f.col_values_.size(); ++j) {
	  os << f.column_functions_[j].value(f.row_values_[i]) << " ";	  
	}
	os << std::endl;
      }
      return os;
    }
    friend std::istream& operator>>(std::istream &is,
				    table<Interpolation>& f) {
      f.header_ = read_header(is);
      size_t n_rows, n_cols;
      is >> n_rows >> n_cols;
      f.row_values_.resize(n_rows);
      for (size_t i=0; i<n_rows; ++i) {
	is >> f.row_values_[i];
      }
      f.col_values_.resize(n_cols);
      for (size_t i=0; i<n_cols; ++i) {
	is >> f.col_values_[i];
      }
      f.column_functions_.resize(n_cols);
      double value;
      for (size_t i=0; i<n_rows; ++i) {
	for (size_t j=0; j<n_cols; ++j) {
	  is >> value;
	  f.column_functions_[j].append({f.row_values_[i],value});	  
	}
      }
      return is;
    }
  };

  using pl_table = table<piecewise_linear>;
  using pp_table = table<piecewise_power>;
  using pe_table = table<piecewise_exponential>;
}

#endif
