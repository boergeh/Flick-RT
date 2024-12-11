#ifndef flick_vector
#define flick_vector

#include <iostream>
#include <valarray>
#include <limits>
#include <fstream>
#include "constants.hpp"

namespace flick {
  class limits {
    double l;
    double u;
  public:
    limits(double lower, double upper) : l(lower), u(upper) {}
    double lower() const {return l;}
    double upper() const {return u;}
  };

  class vector {
  protected:
    double x_;
    double y_;
    double z_;
  public:
    vector() : x_{0}, y_{0}, z_{0} {}
    vector(double x, double y, double z)
      : x_{x}, y_{y}, z_{z} {} 
    static vector from_spherical(double r, double theta, double phi) {
      return vector{r*sin(theta)*cos(phi),
	  r*sin(theta)*sin(phi), r*cos(theta)};
    }
    static vector from_cylindrical(double s, double phi, double z) {
      return vector{s*cos(phi),s*sin(phi),z};
    }
    double x() const {return x_;}
    double y() const {return y_;}
    double z() const {return z_;}
    double r() const {return sqrt(x_*x_ + y_*y_ + z_*z_);}
    double s() const {return sqrt(x_*x_ + y_*y_);}
    double u() const {return z_/r();}
    double theta() const
    // [0, pi]
    {
      return acos(u());
    }
    double phi() const
    // [0, 2*pi]
    {
      double hyp = s();
      if (hyp < std::numeric_limits<double>::epsilon())
      	return 0;
      if (y_ > 0)
    	return acos(x_/hyp);
      return 2*constants::pi-acos(x_/hyp);
    }
    vector& operator+=(const vector& v) {
      x_ += v.x();
      y_ += v.y();
      z_ += v.z();
      return *this;
    }
    vector& operator-=(const vector& v) {
      x_ -= v.x();
      y_ -= v.y();
      z_ -= v.z();
      return *this;
    }
    vector& operator*=(double k) {
      x_ *= k;
      y_ *= k;
      z_ *= k;
      return *this;
    }
    vector& operator/=(double k) {
      x_ /= k;
      y_ /= k;
      z_ /= k;
      return *this;
    }  
  };
  vector operator*(double k, const vector& v) {
    return vector(k*v.x(), k*v.y(), k*v.z());
  }
  vector operator*(const vector& v, double k) {
    return k*v;
  }
  vector operator/(const vector& v, double k) {
    return 1/k*v;
  }
  vector operator+(const vector& v1, const vector& v2) {
    return vector(v1.x()+v2.x(), v1.y()+v2.y(), v1.z()+v2.z());
  }
  vector operator-(const vector& v1, const vector& v2) {
    return v1+(-1*v2);
  }
  vector operator-(const vector& v) {
    return -1*v;
  }
  std::ostream& operator<<(std::ostream& ostr, const vector& v) {
      ostr << v.x() << " "<< v.y() << " "<< v.z();
    return ostr;
  }
  double norm(const vector& v) {
    return v.r();
  }
  double dot(const vector& a, const vector& b) {
    return a.x()*b.x() + a.y()*b.y() + a.z()*b.z();
  }
  vector cross(const vector& a, const vector& b) {
    return vector(a.y()*b.z()-a.z()*b.y(),
		  a.z()*b.x()-a.x()*b.z(),
		  a.x()*b.y()-a.y()*b.x());
  }
  double rms(const vector& a, const vector& b) {
    return sqrt(((a.x()-b.x())*(a.x()-b.x())+
		 (a.y()-b.y())*(a.y()-b.y())+
		 (a.z()-b.z())*(a.z()-b.z()))/3);
  }

  
  class unit_vector : public vector {
  public:
    unit_vector() : vector{0,0,1} {}
    unit_vector(double theta, double phi)
      : vector{vector::from_spherical(1,theta,phi)} {}
    unit_vector(double x, double y, double z)
      : vector{x,y,z} {
      double r = norm(*this);
      x_ /= r;
      y_ /= r;
      z_ /= r;
    }
    // unit_vector(double mu) {
    //  vector{sqrt(fabs(1-mu*mu)), 0, mu};
    // }

    double mu() const {
      return z_/r();
    }
  };

  unit_vector operator-(const unit_vector& uv) {
    return unit_vector(-uv.x(), -uv.y(), -uv.z());
  }
  unit_vector normalize(const vector& v) {
    return unit_vector{v.x(),v.y(),v.z()};
  }

  class rotate {
  private:
    vector v;
    double angle;
  public:
    rotate(const vector& v, double angle) : v(v), angle(angle) {}
    vector about(const unit_vector& k)
    // Rodrigues' rotation formula
    {
      return v*cos(angle)+cross(k,v)*sin(angle)+k*dot(k,v)*(1-cos(angle));
    }
  };
  
  std::valarray<double> cast_to_valarray(const vector& v) {
    return std::valarray<double>{v.x(),v.y(),v.z()};
  }
}
#endif
