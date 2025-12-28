#pragma once

namespace vfep {

constexpr double R_universal = 8.314462618;      // J/(mol*K)
constexpr double sigmaSB     = 5.670374419e-8;   // W/(m^2*K^4)
constexpr double P_atm       = 101325.0;         // Pa

inline double clamp(double x, double lo, double hi) {
    return (x < lo) ? lo : (x > hi) ? hi : x;
}

} // namespace vfep
