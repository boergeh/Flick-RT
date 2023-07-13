#ifndef flick_layered_iops
#define flick_layered_iops

#include "material.hpp"
#include "ab_functions.hpp"

namespace flick {
  class layered_iops {
    material::base& m_;
    size_t n_terms_;
    stdvector oda_;
    stdvector ods_;
    stdvector bottoms_;
    std::vector<std::vector<stdvector>> alpha_;
    std::vector<std::vector<stdvector>> beta_;
    stdvector refidx_;
  public:
    layered_iops(material::base& m, const stdvector& bottoms, size_t n_terms)
      : m_{m}, n_terms_{n_terms}, bottoms_{bottoms},
	oda_(bottoms.size()),
	ods_(bottoms.size()),
	alpha_(4, std::vector<stdvector>(bottoms.size(), stdvector(n_terms))),
	beta_(2, std::vector<stdvector>(bottoms.size(), stdvector(n_terms))),
	refidx_(bottoms.size())
    {
      assert(n_terms > 2);
      for (size_t i=0; i < bottoms_.size(); i++) {
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
    size_t n_layers() {
      return bottoms_.size();
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
    void set_alpha_beta(size_t i) {
      auto [alpha, beta] = material::fitted_mueller_alpha_beta(m_,n_terms_);
      for (size_t n=0; n < alpha_.size(); n++)
	alpha_[n][i] = alpha[n];
      for (size_t n=0; n < beta_.size(); n++)
	beta_[n][i] = beta[n];
    }
    void set_optical_depth(size_t i) {
      double scaling_factor = alpha_[0][i][0]*4*constants::pi;
      oda_[i] = m_.absorption_optical_depth(layer_thickness(i));
      ods_[i] = m_.scattering_optical_depth(layer_thickness(i))*scaling_factor;
    }
    double layer_thickness(size_t i) const {
      if (i > 0)
	return bottoms_[i]-bottoms_[i-1];
      return bottoms_[0];
    }
    void move(double distance) const {
      vector new_position = m_.pose().position() +
	m_.pose().direction()*distance;
      m_.set(m_.pose().move_to(new_position));
    }
    stdvector layer_thicknesses() const {
      stdvector t(bottoms_.size());
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
    friend std::ostream& operator<<(std::ostream &os, accurt_user_specified& u) {
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
      for (size_t i=0; i<u.iops_.n_layers(); i++) {
	for (size_t j=0; j<u.wls_.size(); j++) {
	  os << "A_" << std::to_string(i+1) << "_" << std::to_string(j+1)
	     << " = " << a[i][j] << " #\n";
	  os << "S_" << std::to_string(i+1) << "_" << std::to_string(j+1)
	     << " = " << s[i][j] << " #\n";
	  os << "P_" << std::to_string(i+1) << "_" << std::to_string(j+1)
	     << " = ";
	  for (size_t k=0; k<p[i][j].size(); k++) {
	    os << p[i][j][k]*4*constants::pi << " ";
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
