#ifndef flick_value_collection
#define flick_value_collection

namespace flick {
  class value_collection {
    std::vector<size_t> points_;
    std::vector<double> values_;
    double target_accuracy_;
    size_t total_points_ = 0;
    size_t initial_set_ = 4;
    double noise_floor_ = 0;
  public:
    value_collection(double accuracy)
      : target_accuracy_{accuracy} {
    }
    void add(double value, size_t points) {
      points_.push_back(points);
      values_.push_back(value);
      total_points_ += points;
    }
    void initial_set(size_t n) {
      initial_set_ = n;
    }
    void noise_floor(double value) {
      noise_floor_ = value;
    }
    bool accurate() const {      
      if (values_.size() > initial_set_ and std::isfinite(mean()) and std::isfinite(std()) and not std::isfinite(accuracy()))
	return true;
      if (values_.size() < initial_set_ or not std::isfinite(accuracy()))
	return false;
      return  accuracy() < target_accuracy_; 
    }
    double accuracy() const {
      return std()/(mean()+noise_floor_);
    }
    double mean() const {
      double m = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	m += points_[i]*values_[i];
      }
      return m/total_points_;
    }
    friend std::ostream& operator<<(std::ostream &os, const value_collection& d) {
      os << std::setprecision(5) << d.mean() << " \u00b1 "
	 << std::setprecision(2) << d.std()/d.mean()*100 << "%" 
	 << ", n = " << d.total_points_;
      return os;
    }    
    double std() const {
      double m = mean();
      double s = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	s += points_[i]*pow(values_[i]-m,2);
      }
      return sqrt(s/total_points_);
    }
  }; 
}

#endif
