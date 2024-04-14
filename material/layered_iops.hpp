#ifndef flick_layered_iops
#define flick_layered_iops

#include "material.hpp"
#include "ab_functions.hpp"

namespace flick {
  class layered_iops {
    std::shared_ptr<material::base> m_;
    stdvector boundaries_;
    size_t n_terms_;
    stdvector oda_;
    stdvector ods_;
    std::vector<std::vector<stdvector>> alpha_;
    std::vector<std::vector<stdvector>> beta_;
    stdvector refidx_;
  public:
    layered_iops(std::shared_ptr<material::base> m, const stdvector& boundaries, size_t n_terms)
      : m_{m}, boundaries_{boundaries},
	n_terms_{n_terms},
	oda_(n_layers()),
	ods_(n_layers()),
	alpha_(4, std::vector<stdvector>(n_layers(), stdvector(n_terms))),
	beta_(2, std::vector<stdvector>(n_layers(), stdvector(n_terms))),
	refidx_(n_layers())
    {
      if (boundaries.size()< 2 or boundaries[1] < boundaries[0] or
	  (not std::is_sorted(boundaries.begin(), boundaries.end())))
	throw std::runtime_error("boundary error in layered_iops");
      if (n_terms < 3)
	throw std::runtime_error("number of terms error in layered_iops");
    }
    void set_wavelength(double wl) {
      m_->set_wavelength(wl);
      update();
    }
    size_t n_layers() const {
      return boundaries_.size() - 1;
    }
    stdvector scattering_optical_depth() const {
      return ods_;
    }
    stdvector absorption_optical_depth() const {
      return oda_;
    }
    stdvector single_scattering_albedo() const {
      return (ods_ + oda_) / ods_;
    }
    std::vector<stdvector> alpha_terms(size_t n) const {
      return alpha_.at(n);
    }
    std::vector<stdvector> beta_terms(size_t n) const {
      return beta_.at(n);
    }
    stdvector refractive_index() const {
      return refidx_;
    }
    stdvector absorption_coefficient() const {
      return oda_ / layer_thicknesses();
    }
    stdvector scattering_coefficient() const {
      return ods_ / layer_thicknesses();
    }
    friend std::ostream& operator<<(std::ostream &os,
				    const layered_iops& iops) {
      os << "Scattering optical depth: " << iops.scattering_optical_depth()
	 << "\n";
      os << "Absorption optical depth: " << iops.absorption_optical_depth()
	 << "\n";
      os << "Phase function terms: " << iops.alpha_terms(0)[0] << "\n";
      return os;
    }  
  private:
    void update() {
      m_->set_direction({0,0});
      m_->set_position({0,0,boundaries_[0]});
      for (size_t i=0; i < n_layers(); i++) {
	double h = layer_thickness(i);
	double dh = average_scattering_height(boundaries_[i],boundaries_[i+1])-boundaries_[i];
	move(dh);
	set_alpha_beta(i);
	refidx_[i] = m_->real_refractive_index();
	move(-dh);
	set_optical_depth(i);
	move(h);
      }
    }
    double average_scattering_height(double h_low, double h_high) {
      stdvector h = range(h_low,h_high,100).linspace();
      stdvector s(h.size());
      stdvector sh(h.size());
      vector p0 = m_->pose().position();
      for (size_t i=0; i<sh.size(); i++) {
	m_->set_position({0,0,h[i]});
	s[i] = m_->scattering_coefficient();
	sh[i] = s[i]*h[i];
      }
      m_->set_position(p0);
      pl_function shf(h,sh);
      pl_function sf(h,s);
      double h_avg = shf.integral()/sf.integral();
      if (std::isfinite(h_avg))
	return h_avg;
      return h_low + (h_high-h_low)/2;
    }
    void set_alpha_beta(size_t i) {
      auto [alpha, beta] = material::fitted_mueller_alpha_beta(*m_,n_terms_);
      for (size_t n=0; n < alpha_.size(); n++)
	alpha_[n][i] = alpha[n];
      for (size_t n=0; n < beta_.size(); n++)
	beta_[n][i] = beta[n];
    }
    void set_optical_depth(size_t i) {
      oda_[i] = m_->absorption_optical_depth(layer_thickness(i));
      ods_[i] = m_->scattering_optical_depth(layer_thickness(i));
    }
    double layer_thickness(size_t i) const {
      return boundaries_.at(i+1)-boundaries_.at(i);
    }
    void move(double distance) const {
      vector new_position = m_->pose().position() +
	m_->pose().direction()*distance;
      m_->set_position(new_position);
    }
    stdvector layer_thicknesses() const {
      stdvector t(n_layers());
      for (size_t i=0; i < t.size(); i++)
	t[i] = layer_thickness(i);
      return t;
    }
  };
}

#endif
