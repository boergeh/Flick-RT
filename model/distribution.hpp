#ifndef flick_distribution
#define flick_distribution

namespace flick {
  class distribution {
    std::vector<double> weights_;
    std::vector<double> values_;
    double target_accuracy_{0.1};
    size_t n_packages_{100};
  public:
    distribution(double target_accuracy)
      : target_accuracy_{target_accuracy} {
      n_packages_ = 0.5/pow(target_accuracy,2);
    }
    size_t n_packages() const {
      return n_packages_;
    }
    void add(double value) {
      weights_.push_back(n_packages_);
      values_.push_back(value);
      n_packages_ *= 2;
    }
    bool bad_accuracy() const {
      if (values_.size() < 3 or !(mean() > std::numeric_limits<double>::epsilon()))
	return true;
      return target_accuracy_ < accuracy(); 
    }
    double accuracy() const {
      return std()/mean();
    }
    double mean() const {
      double m=0;
      for (size_t i=0; i<values_.size(); ++i) {
	m += weights_[i]*values_[i];
      }
      return m/total_weight();
    }
    friend std::ostream& operator<<(std::ostream &os, const distribution& d) {
      os << std::setprecision(5) << d.mean() << " \u00b1 "
	 << std::setprecision(2) << d.std() 
	 << ", n = " << d.total_weight();
      return os;
    }    
  protected:
    double std() const {
      double m = mean();
      double s = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	s += weights_[i]*pow(values_[i]-m,2);
      }
      return sqrt(s/total_weight());
    }
    double total_weight() const {
      double w = 0;
      for (size_t i=0; i<weights_.size(); ++i) {
	w += weights_[i];
      }
      return w;
    }
  };
}

#endif
