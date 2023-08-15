#ifndef flick_numeric_units
#define flick_numeric_units

namespace flick {
  namespace units {
    using ld = long double;
    using ll = unsigned long long;

    // length
    constexpr ld operator""_nm(ld x) {return x*1e-9;}
    constexpr ld operator""_nm(ll x) {return x*1e-9;}
    constexpr ld operator""_mm(ld x) {return x*1e-3;}
    constexpr ld operator""_mm(ll x) {return x*1e-3;}
    constexpr ld operator""_m(ld x) {return x;}
    constexpr ld operator""_m(ll x) {return x;}
    constexpr ld operator""_km(ld x) {return x*1e3;}
    constexpr ld operator""_km(ll x) {return x*1e3;}

    // temperature
    constexpr ld operator""_K(ld x) {return x;}
    constexpr ld operator""_K(ll x) {return x;}

    // pressure
    constexpr ld operator""_Pa(ld x) {return x;}
    constexpr ld operator""_Pa(ll x) {return x;}
    constexpr ld operator""_hPa(ld x) {return 100*x;}
    constexpr ld operator""_hPa(ll x) {return 100*x;}

    // other
    constexpr ld operator""_psu(ld x) {return x;}
    constexpr ld operator""_psu(ll x) {return x;}
  }
}

#endif
