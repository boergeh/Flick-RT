#ifndef flick_std_operators
#define flick_std_operators
#include <vector>

namespace flick {
  using stdvector = std::vector<double>;

  stdvector& operator+=(stdvector& v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] += v2[i];
    return v1;
  }
  stdvector operator+(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] += v2[i];
    return v1;
  }
  stdvector operator+(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] += k;
    return v;
  }
  stdvector operator+(double k, stdvector v) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] += k;
    return v;
  }
  stdvector operator-(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] -= v2[i];
    return v1;
  }
  stdvector operator-(stdvector v1, double v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] -= v2;
    return v1;
  }
  stdvector operator*(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] *= v2[i];
    return v1;
  }
  stdvector operator*(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] *= k;
    return v;
  }
  stdvector operator*(double k, stdvector v) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] *= k;
    return v;
  }
  stdvector operator/(stdvector v1, const stdvector& v2) {  
    for (size_t i=0; i<v1.size(); ++i)
      v1[i] /= v2[i];
    return v1;
  }
  stdvector operator/(stdvector v, double k) {  
    for (size_t i=0; i<v.size(); ++i)
      v[i] /= k;
    return v;
  }
  double sum(const stdvector& v) {
    double s=0;
    for (size_t i=0; i<v.size(); ++i)
      s += v[i];
    return s;
  }
  stdvector operator^(stdvector v, double k) {
    double s=0;
    for (size_t i=0; i<v.size(); ++i)
      v[i] = std::pow(v[i],k);
    return v;
  }
  double rms(const stdvector& v) {
    return sqrt(pow(sum(v^2),1./2));
  }

}
#endif
