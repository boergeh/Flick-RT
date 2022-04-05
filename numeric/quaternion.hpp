#ifndef flick_quaternion
#define flick_quaternion

#include <iostream>

namespace flick {
  // see Wikipedia quaternion
  struct quaternion {    
    double a{0};
    double b{0};
    double c{0};
    double d{0};
    quaternion(double a_, double b_, double c_, double d_)
      : a(a_), b(b_), c(c_), d(d_){}
    quaternion& operator+=(const quaternion& q) {
      a += q.a;
      b += q.b;
      c += q.c;
      d += q.d;
      return *this;
    }
    quaternion& operator-=(const quaternion& q) {
      a -= q.a;
      b -= q.b;
      c -= q.c;
      d -= q.d;
      return *this;
    }
  };
  quaternion operator*(double l, const quaternion& q) {
    return quaternion{l*q.a, l*q.b, l*q.c,l*q.d};
  }
  quaternion operator*(const quaternion& q, double l) {
    return l*q;
  }
  quaternion operator*(const quaternion& q1, const quaternion& q2) {
    return quaternion{q1.a*q2.a-q1.b*q2.b-q1.c*q2.c-q1.d*q2.d,
	q1.a*q2.b+q1.b*q2.a+q1.c*q2.d-q1.d*q2.c,
	q1.a*q2.c-q1.b*q2.d+q1.c*q2.a+q1.d*q2.b,
	q1.a*q2.d+q1.b*q2.c-q1.c*q2.b+q1.d*q2.a};
  }
  quaternion operator+(const quaternion& q1, const quaternion& q2) {
    quaternion q = q1;
    return q += q2;
  }
  quaternion operator-(const quaternion& q1, const quaternion& q2) {
    return q1+(-1*q2);
  }
  quaternion operator-(const quaternion& q) {
    return -1*q;
  }
  quaternion inv(const quaternion& q) {
    return 1/(q.a*q.a+q.b*q.b+q.c*q.c+q.d*q.d)*quaternion{q.a,-q.b,-q.c,-q.d};
  }
  double norm(const quaternion& q) {
    return sqrt(q.a*q.a + q.b*q.b + q.c*q.c + q.d*q.d);
  }
  quaternion conj(const quaternion& q) {
    return quaternion{q.a,-q.b,-q.c,-q.d};
  }
  quaternion identity() {
    return quaternion{1,0,0,0};
  }
  std::ostream& operator<<(std::ostream& ostr, const quaternion& q) {
      ostr << q.a << " "<< q.b << " "<< q.c << " "<< q.d;
    return ostr;
  }

}
#endif
