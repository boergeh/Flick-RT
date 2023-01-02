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
      n_packages_ = 0.5*1/pow(target_accuracy,2);
    }
    size_t n_packages() {
      return n_packages_;
    }
    void add(double value) {
      weights_.push_back(n_packages_);
      values_.push_back(value);
      n_packages_ *= 2;
    }
    bool bad_accuracy() {
      if (values_.size() < 3)
	return true;
      else if (mean() < std::numeric_limits<double>::epsilon())
	return false;
      return target_accuracy_ < std()/mean(); 
    }
    double mean() {
      double m=0;
      for (size_t i=0; i<values_.size(); ++i) {
	m += weights_[i]*values_[i];
      }
      return m/total_weight();
    }
  protected:
    double std() {
      double m = mean();
      double s = 0;
      for (size_t i=0; i<values_.size(); ++i) {
	s += weights_[i]*pow(values_[i]-m,2);
      }
      return sqrt(s/total_weight());
    }
    double total_weight() {
      double w = 0;
      for (size_t i=0; i<weights_.size(); ++i) {
	w += weights_[i];
      }
      return w;
    }
  };
}

#endif
