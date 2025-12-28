#pragma once

#include <algorithm>

namespace vfep {

struct LiIonConfig {
    double mass_pack_kg       = 8.0;
    double T_onset_K          = 423.15;
    double T_full_K           = 773.15;
    double E_total_J_per_kg   = 6.0e6;
    double tau_s              = 35.0;
    double ventGas_kg_per_kg  = 0.12;
};

struct LiIonState {
    bool enabled = false;
    bool triggered = false;
    double energyRemaining_J = 0.0;
    double ventGasRemaining_kg = 0.0;
};

class LiIonRunaway {
public:
    void configure(const LiIonConfig& cfg) { cfg_ = cfg; reset(); }

    // Resets energy/vent inventories and trigger state.
    // Does NOT change enabled flag.
    void reset() {
        st_.triggered = false;
        st_.energyRemaining_J = cfg_.mass_pack_kg * cfg_.E_total_J_per_kg;
        st_.ventGasRemaining_kg = cfg_.mass_pack_kg * cfg_.ventGas_kg_per_kg;
    }

    void setEnabled(bool e) { st_.enabled = e; }
    const LiIonState& state() const { return st_; }

    void step(double dt, double T_K, double& heat_W, double& ventGas_kgps);

private:
    LiIonConfig cfg_;
    LiIonState st_;
};

} // namespace vfep
