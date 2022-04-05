#ifndef flick_emitter
#define flick_emitter

#include "radiation_package.hpp"
#include <initializer_list>

namespace flick {
  
  class wavelength_distribution {
  protected:
    uniform_random ur_;
  public:
    virtual double draw() = 0;
    virtual ~wavelength_distribution()=default;
  };
  class monocromatic : public wavelength_distribution {
    double wl_;
  public:
    monocromatic(double wl) : wl_{wl} {} 
    double draw() {
      return wl_;
    }
  };
  //class gaussian_solar_spectrum : public wavelength_distribution {
  //class gaussian_smooth_spectrum : public wavelength_distribution {
  
  class direction_distribution {
  protected:
    direction_generator dg_;
  public:
    virtual unit_vector draw() = 0;
    virtual ~direction_distribution()=default;
  };

  class unidirectional : public direction_distribution {
    unit_vector direction_;
  public:
    unidirectional(const unit_vector& d) : direction_{d}{}
    unit_vector draw() {
      return direction_;
    }
  };

  class conic : public direction_distribution {
    unit_vector cone_direction_;
    double solid_angle_;
  public:
    conic(double solid_angle, const unit_vector& cone_direction)
      : cone_direction_{cone_direction}, solid_angle_{solid_angle} {} 
    unit_vector draw() {
      return dg_.conic(solid_angle_, cone_direction_);
    }
  };
  
  class isotropic : public direction_distribution {
  public:
    isotropic() {} 
    unit_vector draw() {
      return dg_.isotropic();
    }
  };

  class emitter
  // Point radiometric device
  {
    pose placement_;
    size_t total_packages_;
    size_t packages_left_;
    stokes initial_stokes_;
    std::shared_ptr<wavelength_distribution> wld_;
    std::shared_ptr<direction_distribution> dd_;
  public:
    emitter(size_t packages, const stokes& initial_stokes)
      : total_packages_{packages}, packages_left_{packages},
	initial_stokes_{initial_stokes},
	wld_{std::make_shared<monocromatic>(500e-9)},
	dd_{std::make_shared<isotropic>()}
    {}
    emitter& move_by(const vector& v) {
      placement_.move_by(v);
      return *this;
    }   
    emitter& rotate_by(const quaternion& rotation) {
      placement_.rotate_by(rotation);
      return *this;
    }   
    template <class Wd, class... Args>
    emitter& wavelength(Args... a) {
      wld_ = std::make_shared<Wd>(a...);
      return *this;
    }
    template <class Dd, class... Args>
    emitter& direction(Args... a) {
      dd_ = std::make_shared<Dd>(a...);
      return *this;
    }   
    emitter& refill() {
      packages_left_ = total_packages_;
      return *this;
    }    
    radiation_package emit()
    {
      radiation_package rp(wld_->draw(), initial_stokes_);      
      rp.move_to(placement_.position());
      rp.rotate_to(dd_->draw());
      --packages_left_;
      return rp;
    }    
    size_t packages_left() const {
      return packages_left_;
    }
    friend std::ostream& operator<<(std::ostream &os, const emitter& em) {
      os << em.placement_ << " " << em.packages_left_ << " " << em.initial_stokes_;
      return os;
    }
  };
}
#endif
