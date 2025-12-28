#include "LiIonRunaway.h"

#include <algorithm>
#include <cmath>

namespace vfep {

namespace {
static double clamp01(double x) {
    if (x < 0.0) return 0.0;
    if (x > 1.0) return 1.0;
    return x;
}

static bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}
} // namespace

void LiIonRunaway::step(double dt, double T_K, double& heat_W, double& ventGas_kgps) {
    heat_W = 0.0;
    ventGas_kgps = 0.0;

    if (!st_.enabled) return;
    if (!isFinitePositive(dt)) return;
    if (!std::isfinite(T_K)) return;

    if (st_.energyRemaining_J <= 1.0 && st_.ventGasRemaining_kg <= 1e-9) return;

    if (!st_.triggered && T_K >= cfg_.T_onset_K) st_.triggered = true;
    if (!st_.triggered) return;

    const double denom = std::max(1e-9, (cfg_.T_full_K - cfg_.T_onset_K));
    const double prog = clamp01((T_K - cfg_.T_onset_K) / denom);

    // Rate ramps from 0.15/tau to 1.0/tau as temperature approaches full runaway
    const double tau = std::max(1e-6, cfg_.tau_s);
    const double rate = (0.15 + 0.85 * prog) / tau; // 1/s

    // dt-independent first-order decay: dX = X*(1 - exp(-rate*dt))
    const double oneMinus = 1.0 - std::exp(-rate * dt);

    const double dE = std::min(st_.energyRemaining_J, st_.energyRemaining_J * oneMinus);
    st_.energyRemaining_J -= dE;
    heat_W = dE / dt;

    const double dM = std::min(st_.ventGasRemaining_kg, st_.ventGasRemaining_kg * oneMinus);
    st_.ventGasRemaining_kg -= dM;
    ventGas_kgps = dM / dt;
}

} // namespace vfep
