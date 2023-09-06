#ifndef flick_distribution
#define flick_distribution

#include "constants.hpp"

namespace flick {
namespace distribution {

  template<class F>
  class newton_raphson {
    const F& f_;
    double percent_accuracy_;
    size_t max_iterations_ = 1e4;
    double epsilon_ = std::numeric_limits<double>::epsilon();
  public:
    newton_raphson(const F& f, double percent_accuracy = 1e-3)
      : f_{f}, percent_accuracy_{percent_accuracy} {}
    std::optional<double> hunt(double y_goal, double x_guess) const {
      double error = std::numeric_limits<double>::max();
      std::optional<double> x = x_guess;
      double previous_x = *x;
      double previous_error = error;
      size_t n=0;
      while(error > percent_accuracy_) {
	double y = f_.value(*x);
	x = *x - (y-y_goal) / f_.derivative(*x);
	error = 100*fabs(2*(*x-previous_x)/(*x+previous_x+epsilon_));
	if (not std::isfinite(error) or error > previous_error or n > max_iterations_)
	  return std::nullopt;
	previous_x = *x;
	previous_error = error;
	n++;
      }
      //std::cout << n-1 << std::endl;
      return x;
    }
  };

  class erf {
  public:
    double value(double x) const {
      return std::erf(x);
    }
    double derivative(double x) const {
      return 2/sqrt(constants::pi)*exp(-pow(x,2));
    }
  };

  class erf_inv {
    newton_raphson<erf> nr_;
    mutable double x_guess_ = 0;
  public:
    erf_inv(double percent_accuracy=1e-5) 
      : nr_(erf(),percent_accuracy) {}
    double operator()(double x) const {
      double y_goal = x;
      std::optional<double> xo = nr_.hunt(y_goal, x_guess_);
      if (not xo.has_value()) {
	x_guess_ = 0;
	xo = nr_.hunt(y_goal, x_guess_);
      }
      x_guess_ = *xo;
      return *xo;
    }
  };
  
  class basic_distribution {
  protected:
    const double pi_ = constants::pi;
  public:
    virtual double pdf(double x) const = 0;
    virtual double cdf(double x) const = 0;
    virtual double quantile(double x) const = 0;
  };
  
  class normal : public basic_distribution
  // https://en.wikipedia.org/wiki/Normal_distribution
  {
    double mu_;
    double sigma_;
    erf_inv erf_inv_;
  public:
    normal(double mu, double sigma)
      : mu_{mu}, sigma_{sigma} {}
    double pdf(double x) const {
      return 1/(sigma_*sqrt(2*pi_))*exp(-0.5*pow((x-mu_)/sigma_,2));
    }
    double cdf(double x) const {
      return 0.5*(1+std::erf((x-mu_)/(sigma_*sqrt(2))));
    }
    double quantile(double p) const {
      return mu_+sigma_*sqrt(2)*erf_inv_(2*p-1);
    }
  };

   class log_normal : public basic_distribution
  // https://en.wikipedia.org/wiki/Log-normal_distribution
  {
    double mu_;
    double sigma_;
    erf_inv erf_inv_;
  public:
    log_normal(double mu, double sigma)
      : mu_{mu}, sigma_{sigma} {}
    double pdf(double x) const {
      return 1/(x*sigma_*sqrt(2*pi_))*exp(-pow(log(x)-mu_,2)/(2*pow(sigma_,2)));
    }
    double cdf(double x) const {
      return 0.5*(1+std::erf((log(x)-mu_)/(sigma_*sqrt(2))));
    }
    double quantile(double p) const {
      return exp(mu_+sqrt(2*pow(sigma_,2))*erf_inv_(2*p-1));
    }
  };

  // henyey_greenstein etc here with mu as input x?
  // rieman zeta function for planck. or just planck. yes.
  
  stdvec quantiles(const basic_distribution& d, size_t n_points) {
    double dx = 0.1/(n_points+2);
    stdvec p = range(dx,1-dx,n_points).linspace();
    stdvec x(n_points);
    for (size_t i=0; i < n_points; ++i) {
      x[i] = d.quantile(p[i]);
    }
    return x;
  }
}
}

#endif
