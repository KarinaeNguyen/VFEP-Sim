#pragma once

#include "Reactor.h"

namespace vfep {

struct SuppressionConfig {
    bool enabled = false;
    double mdot_total_kgps = 0.30;

    // Mass fractions of total discharge. Recommended: sum to 1.0.
    double frac_inhibitor = 0.25;
    double frac_inert     = 0.75;

    // Cooling effect modeled as a lumped removal term (W per kg/s of discharge)
    double cooling_W_per_kgps = 250000.0;

    // Minimal actuator model (truth telemetry)
    // - rpm_target: commanded steady-state RPM when enabled and tank > 0
    // - rpm_ramp_rate_rpmps: first-order ramp rate (RPM per second)
    double rpm_target = 3200.0;
    double rpm_ramp_rate_rpmps = 4000.0;
};

struct SuppressionState {
    double inhibitor_kg = 0.0;
    double inert_kg     = 0.0;
    double tank_kg      = 200.0;

    // Actuator truth state
    double vfep_rpm = 0.0;
};

class Suppression {
public:
    explicit Suppression(const ChemistryIndex& idx) : idx_(idx) {}

    void setConfig(const SuppressionConfig& cfg) { cfg_ = cfg; }
    const SuppressionConfig& config() const { return cfg_; }

    void resetTank(double tank_kg) {
        st_.tank_kg = tank_kg;
        st_.inhibitor_kg = 0.0;
        st_.inert_kg = 0.0;
        st_.vfep_rpm = 0.0;
    }

    // Returns cooling power in Watts (positive = cooling).
    double apply(double dt, Reactor& r,
                 double& inhib_kgm3,
                 double& inert_kgm3,
                 double& agent_mdot_kgps);

    // Actuator truth telemetry (do not infer from mdot in the visualizer).
    double vfepRPM() const noexcept { return st_.vfep_rpm; }

private:
    SuppressionConfig cfg_;
    SuppressionState st_;
    ChemistryIndex idx_;
};

} // namespace vfep
