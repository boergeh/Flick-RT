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
    vector position_{0,0,0};
    size_t total_packages_{0};
    size_t packages_left_{0};
    stokes initial_stokes_{1,0,0,0};
    std::shared_ptr<wavelength_distribution> wld_;
    std::shared_ptr<direction_distribution> dd_;
  public:
    emitter() = default;
    emitter(vector position, const stokes& initial_stokes, size_t n_packages)
      : position_{position},
	total_packages_{n_packages},
	packages_left_{n_packages},
	initial_stokes_{initial_stokes},
	wld_{std::make_shared<monocromatic>(500e-9)},
	dd_{std::make_shared<isotropic>()} {
    }
    emitter(size_t n_packages) : emitter({0,0,0,},{1,0,0,0},n_packages) {
    }
    template <class Wd, class... Args>
    void set_wavelength(Args... a) {
      wld_ = std::make_shared<Wd>(a...);
    }
    template <class Dd, class... Args>
    void set_direction(Args... a) {
      dd_ = std::make_shared<Dd>(a...);
    }   
    void add_packages(size_t n) {
      total_packages_ += n;
      packages_left_ += n;
    }    
    radiation_package emit()
    {
      pose p{position_, dd_->draw()};
      radiation_package rp(p, initial_stokes_);
      rp.wavelength(wld_->draw());
      packages_left_--;
      return rp;
    }    
    size_t packages_left() const {
      return packages_left_;
    }
    bool is_empty() const {
      return (packages_left_ == 0);
    }
    friend std::ostream& operator<<(std::ostream &os, const emitter& em) {
      os << em.position_ << " " << em.initial_stokes_  << " " << em.packages_left_;
      return os;
    }
  };
}

#endif
