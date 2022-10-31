#ifndef flick_sorted_vector
#define flick_sorted_vector

#include "../environment/exception.hpp"

namespace flick {
  enum class step_type{linear, exponential};
  class sorted_vector
  {
    mutable size_t current_lookup_{0};
    step_type step_type_{step_type::linear};
    std::vector<double> values_;
  public:

    class iterator {
    protected:
      size_t current_bin{0};
      int directed_step;
      sorted_vector *sv_;
    public:
      void move_to_next_bin() {
	current_bin += directed_step;
      }
      void move_to_bin_at(double value) {
	current_bin = sv_->find(value);
      }
      virtual size_t next_index() = 0;
      virtual size_t previous_index() = 0;
      virtual bool is_in_end_bin() = 0;
    };    
    class ascending_iterator : public iterator {
    public:
      ascending_iterator(sorted_vector& sv){
	sv_ = &sv;
	directed_step = 1;
      }
      size_t next_index() {
	return current_bin + 1;
      }
      size_t previous_index() {
	return current_bin;
      }
      bool is_in_end_bin() {
	if (current_bin == sv_->values_.size()-2)
	  return true;
	return false;
      }
    };
    class descending_iterator : public iterator {
    public:
      descending_iterator(sorted_vector& sv){
	sv_ = &sv;
	directed_step = -1;
      }
      size_t next_index() {
	return current_bin;
      }
      size_t previous_index() {
	return current_bin + 1;
      }
      bool is_in_end_bin() {
	if (current_bin == 0)
	  return true;
	return false;
      }
    };
    
    sorted_vector(){}
    sorted_vector(const std::vector<double>& v)
      : values_{v} {
      for (size_t i = 0; i < v.size()-1; ++i)
	ensure(v[i] < v[i+1]);
    }
    void set_step_type(step_type st) {
      if (st == step_type::exponential)
	ensure(values_.at(0) > 0);
      step_type_ = st;
    }
    void append(double v) {
      if (values_.size() > 0)
	ensure(v > values_.back());
      values_.push_back(v);
    } 
    void scale(double factor) {
      ensure(factor > 0.0);
      for (auto& v : values_)
	v *= factor;
    }
    void shift(double distance) {
      for (auto& v : values_)
	v += distance;
    }
    void log_transform() {
      ensure(values_.at(0) > 0);
      for (auto& v : values_)
	v = log(v);
    }    
    void exp_transform() {
      for (auto& v : values_)
	v = exp(v);
    }    
    size_t find(double value) const
    // find index of close element using last lookup and rate of
    // change
    {
      if (below_boundary(value))
	return 0;
      else if (above_boundary(value))
	return values_.size()-2;
      size_t lower_bound = 0;
      size_t upper_bound = values_.size()-2;
      size_t n = current_lookup_;   
      while (!inside_bin(value,n)) {
	double step = estimate_number_of_bins_to_skip(n,value);
	if (step > 0)
	  lower_bound = n+1;
	else if (step < 0)
	  upper_bound = n-1; 
	n += step;
	if (n < lower_bound)
	  n = lower_bound;
	else if (n > upper_bound)
	  n = upper_bound;
      }
      return n;
    }
    double operator[](size_t n) const {
      return values_[n];
    }
    size_t size() const {
      return values_.size();
    }
    double upper_value() const {
      return values_.back();
    }
    double lower_value() const {
      return values_.front();
    }
    std::vector<double> all_values() const {
      return values_;
    }
  private:
    void ensure(bool b) const {
      if (!b)
	throw exception("sorted_vector");
    }
    double estimate_number_of_bins_to_skip(size_t n, double value) const {
      if (step_type_ == step_type::linear)
	 return (value - values_[n]) / (values_[n+1]-values_[n]);
      else if (step_type_ == step_type::exponential && value > 0)
	return log(value/values_[n]) / log(values_[n+1]/values_[n]);
      else
	throw std::invalid_argument("sorted_vector");
    }   
    bool inside_bin(double value, size_t n) const {
      if (value >= values_[n] && value < values_[n+1]) {
	current_lookup_ = n;
	return true;
      }
      return false;
    }
    bool below_boundary(double value) const {
      if (value <= values_[0])
	return true;
      return false;
    }
    bool above_boundary(double value) const {
      if (value >= values_.end()[-2])
	return true;
      return false;
    }    
  };
  
  std::ostream& operator<<(std::ostream& os, const sorted_vector& v) {
    for (size_t i=0; i<v.size(); ++i)
      os << v[i] << " ";
    return os;
  }
}
#endif
