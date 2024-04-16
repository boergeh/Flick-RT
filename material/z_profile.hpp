#ifndef flick_material_z_profile
#define flick_material_z_profile

#include "material.hpp"
#include "iop_profile.hpp"
#include "../numeric/function.hpp"

namespace flick {
namespace material {
  template<class Function>
  class z_profile : public base {
  protected:
    iop_z_profile<Function> a_profile_;
    iop_z_profile<Function> s_profile_;
    double real_refractive_index_{1};
  public:
    z_profile() = default;
    const iop_z_profile<Function>& a_profile() const {
      return a_profile_;
    }
    const iop_z_profile<Function>& s_profile() const {
      return s_profile_;
    }
    const stdvector& height_grid() const {
      return a_profile_.height_grid();
    }
    virtual double real_refractive_index() const {
      return real_refractive_index_;
    } 
    double absorption_coefficient() const {
      return a_profile_.value(pose().position().z());
    }
    double scattering_coefficient() const {
      return s_profile_.value(pose().position().z());
    }
    double absorption_optical_depth(double distance) const {
      return a_profile_.optical_depth(pose(),distance);
    }
    double scattering_optical_depth(double distance) const {
      return s_profile_.optical_depth(pose(),distance);
    }
    double absorption_distance(double absorption_optical_depth) const {
      return s_profile_.distance(pose(),absorption_optical_depth);
    }
    double scattering_distance(double scattering_optical_depth) const {
      return s_profile_.distance(pose(),scattering_optical_depth);
    }
  private:
    friend std::ostream& operator<<(std::ostream &os, const z_profile& p) {
      os << "\n";
      os << "Absorption coefficient profile" << "\n";
      os << p.a_profile_;
      os << "Scattering coefficient profile" << "\n";
      os << p.s_profile_;
      return os;
    }
  };

  
  template<class Function>
  class scaled_z_profile : public z_profile<Function> {
  private:
    std::shared_ptr<base> m_;
    stdvector z_;
    stdvector scaling_factor_;
  public:
    scaled_z_profile(const std::shared_ptr<base>& m, const stdvector& z,
		     const stdvector& scaling_factor)
      : m_{m}, z_{z}, scaling_factor_{scaling_factor} {
      if (z_.size() != scaling_factor_.size())
	throw std::runtime_error("scaled_z_profile");
      make_iop_profile();
    } 
    void set_wavelength(double wl) override {
      m_->set_wavelength(wl);
      make_iop_profile();
    }
    void make_iop_profile() {
      z_profile<Function>::real_refractive_index_ = m_->real_refractive_index();
      stdvector a(z_.size());
      stdvector s(z_.size());
      vector p0 = m_->pose().position();
      for (size_t i=0; i<z_.size(); ++i) {
	m_->set_position({p0.x(), p0.y(), z_[i]});
	a[i] = m_->absorption_coefficient()*scaling_factor_[i];
	s[i] = m_->scattering_coefficient()*scaling_factor_[i];
      }
      m_->set_position(p0);
      z_profile<Function>::a_profile_ = iop_z_profile<Function>(Function(z_,a));
      z_profile<Function>::s_profile_ = iop_z_profile<Function>(Function(z_,s));
    }
    mueller mueller_matrix(const unit_vector& scattering_direction) const override {
      return m_->mueller_matrix(scattering_direction);
    }
  };
 
  template<class Function>
  std::shared_ptr<scaled_z_profile<Function>> make_scaled_z_profile(const std::shared_ptr<base>& m,
								    const stdvector& z,
								    const stdvector& scaling_factor) {
    return std::make_shared<scaled_z_profile<Function>>(scaled_z_profile<Function>(m,z,scaling_factor));   
  }
}
}

#endif
