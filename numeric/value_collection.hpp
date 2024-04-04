#ifndef flick_value_collection
#define flick_value_collection

namespace flick {
  class value_collection {
    std::vector<size_t> n_points_;
    std::vector<double> values_;
    double target_accuracy_;
    size_t total_n_points_ = 0;
    size_t initial_set_ = 4;
    double noise_floor_ = 0;
    bool log_values_ = false;
  public:
    value_collection(double accuracy)
      : target_accuracy_{accuracy} {
    }
    void set_log_values(bool b) {
      log_values_ = b;
    }
    void add(double value, size_t n_points) {
      n_points_.push_back(n_points);
      values_.push_back(value);
      total_n_points_ += n_points;
    }
    void initial_set(size_t n) {
      initial_set_ = n;
    }
    void noise_floor(double value) {
      noise_floor_ = value;
    }
    bool accurate() const {
      if (values_.size() < initial_set_)
	return false;
      if (std::isfinite(mean()) and std::isfinite(std()) and not std::isfinite(accuracy()))
	return true;
      if (not std::isfinite(accuracy()))
	return false;
      return  accuracy() < target_accuracy_; 
    }
    bool last_accurate() const {
      if (values_.size() < initial_set_)
	return false;
      if (not std::isfinite(accuracy()))
	return true;
      return  last_accuracy() < target_accuracy_; 
    }
    double last_accuracy() const {
      if (values_.size() > 1) {
	double v2 = values_.end()[-1];
	double v1 = values_.end()[-2];
	return fabs(v2-v1)/(v1+noise_floor_);
      }
      return std::numeric_limits<double>::max();
    }
    double accuracy() const {
      double mu = mean();
      double sigma = std();
      if (log_values_) {
	mu = exp(log_mean());
	sigma = mu*(exp(log_std())-exp(-log_std()));
      }
      return sigma/(mu+noise_floor_);
    }
    double last_value() const {
      return values_.back();
    }
    double mean() const {
      double m = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	m += n_points_[i]*values_[i];
      }
      return m/total_n_points_;
    }
    double log_mean() const {
      double m = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	m += n_points_[i]*log(values_[i]);
      }
      return m/total_n_points_;
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const value_collection& d) {
      os << std::setprecision(5) << d.last_value() << " \u00b1 "
	 << std::setprecision(2) << d.last_accuracy()*100 << "%" 
	 << " n = " << d.n_points_.back();
      os << std::setprecision(6);
      return os;
    }    
    double std() const {
      double m = mean();
      double s = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	s += n_points_[i]*pow(values_[i]-m,2);
      }
      return sqrt(s/total_n_points_);
    }
    double log_std() const {
      double m = log_mean();
      double s = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	s += n_points_[i]*pow(log(values_[i])-m,2);
      }
      return sqrt(s/total_n_points_);
    }
  }; 
}

#endif
