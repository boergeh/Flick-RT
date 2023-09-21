#ifndef flick_layered_iops
#define flick_layered_iops

#include "material.hpp"
#include "ab_functions.hpp"

namespace flick {
  class layered_iops {
    material::base& m_;
    stdvector boundaries_;
    size_t n_terms_;
    stdvector oda_;
    stdvector ods_;
    std::vector<std::vector<stdvector>> alpha_;
    std::vector<std::vector<stdvector>> beta_;
    stdvector refidx_;
  public:
    layered_iops(material::base& m, const stdvector& boundaries, size_t n_terms)
      : m_{m}, boundaries_{boundaries},
	n_terms_{n_terms},
	oda_(n_layers()),
	ods_(n_layers()),
	alpha_(4, std::vector<stdvector>(n_layers(), stdvector(n_terms))),
	beta_(2, std::vector<stdvector>(n_layers(), stdvector(n_terms))),
	refidx_(n_layers())
    {
      m_.set_direction({0,0});
      if (n_terms <= 2)
	throw std::runtime_error("layered_iops");
      for (size_t i=0; i < n_layers(); i++) {
	move(layer_thickness(i)/2);
	set_alpha_beta(i);
	refidx_[i] = m_.real_refractive_index();
	move(-layer_thickness(i)/2);
	set_optical_depth(i);
	move(layer_thickness(i));
      }
    }
    void set_wavelength(double wl) {
      m_.set_wavelength(wl);
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
  private:
    friend std::ostream& operator<<(std::ostream &os,
				    const layered_iops& iops) {
      os << "Scattering optical depth: " << iops.scattering_optical_depth() << "\n";
      os << "Absorption optical depth: " << iops.absorption_optical_depth() << "\n";
      os << "Phase function terms: " << iops.alpha_terms(0)[0] << "\n";
      return os;
    }  
    void set_alpha_beta(size_t i) {
      auto [alpha, beta] = material::fitted_mueller_alpha_beta(m_,n_terms_);
      for (size_t n=0; n < alpha_.size(); n++)
	alpha_[n][i] = alpha[n];
      for (size_t n=0; n < beta_.size(); n++)
	beta_[n][i] = beta[n];
    }
    void set_optical_depth(size_t i) {
      m_.set_position({0,0,boundaries_[i]});
      oda_[i] = m_.absorption_optical_depth(layer_thickness(i));
      ods_[i] = m_.scattering_optical_depth(layer_thickness(i));
    }
    double layer_thickness(size_t i) const {
      return boundaries_.at(i+1)-boundaries_.at(i);
    }
    void move(double distance) const {
      vector new_position = m_.pose().position() +
	m_.pose().direction()*distance;
      m_.set_position(new_position);
    }
    stdvector layer_thicknesses() const {
      stdvector t(n_layers());
      for (size_t i=0; i < t.size(); i++)
	t[i] = layer_thickness(i);
      return t;
    }
  };

  class accurt_user_specified {
    layered_iops& iops_;
    stdvector wls_;
  public:
    accurt_user_specified(layered_iops& iops, const stdvector& wavelengths)
      : iops_{iops}, wls_{wavelengths} {}
  private:
    friend std::ostream& operator<<(std::ostream &os, const accurt_user_specified& u) {
      os << "# AccuRT configuration file for the user_specified material #\n"
	 << "PROFILE_LABEL = layer_numbering #\n"
	 << "MATERIAL_PROFILE = 1 #\n"
	 << "TURN_OFF_DELTA_FIT = true #\n\n"
	 << "WAVELENGTHS = ";
      for (auto& wl:u.wls_)
	os << wl*1e9 << " ";
      os << "#\n\n";
      std::vector<stdvector> n;
      std::vector<stdvector> a;
      std::vector<stdvector> s;
      std::vector<std::vector<stdvector>> p;
      for (auto& wl:u.wls_) {
	u.iops_.set_wavelength(wl);
	for (size_t j=0; j<u.iops_.n_layers(); j++) {
	  n.push_back(u.iops_.refractive_index());
	  a.push_back(u.iops_.absorption_coefficient());
	  s.push_back(u.iops_.scattering_coefficient());
	  p.push_back(u.iops_.alpha_terms(0));
	}
      }
      os << "REFRACTIVE_INDEX = ";
      for (size_t i=0; i<u.wls_.size(); i++) {
	os << *std::max_element(std::begin(n[i]), std::end(n[i])) << " ";
      }
      os << "#\n\n";
      size_t l = u.iops_.n_layers();
      for (int i=l-1; i>=0; i--) {
	for (size_t j=0; j<u.wls_.size(); j++) {
	  os << "A_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << a[i][j] << " #\n";
	  os << "S_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = " << s[i][j] << " #\n";
	  os << "P_" << std::to_string(l-i) << "_" << std::to_string(j+1)
	     << " = ";
	  for (size_t k=0; k<p[i][j].size(); k++) {
	    os << p[i][j][k]/p[i][j][0]/(2*k+1) << " ";
	  }
	  os << " #\n";
	  os << "\n";
	}
      }  
      return os;
    }
  };
}

#endif
