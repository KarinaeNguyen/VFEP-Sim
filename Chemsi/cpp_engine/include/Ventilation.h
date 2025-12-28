#pragma once

#include "Reactor.h"

namespace vfep {

struct VentilationConfig {
    // Air changes per hour
    double ACH = 3.0;

    // Supply temperature
    double T_supply_K = 295.15;
};

class Ventilation {
public:
    explicit Ventilation(const ChemistryIndex& idx) : idx_(idx) {}

    void setConfig(const VentilationConfig& cfg) { cfg_ = cfg; }
    const VentilationConfig& config() const { return cfg_; }

    void apply(double dt, Reactor& r);

private:
    VentilationConfig cfg_;
    ChemistryIndex idx_;
};

} // namespace vfep
