// Suppression.cpp (patched: inert_kgm3 derived from reactor state to avoid divergence)
#include "Suppression.h"

#include <algorithm>
#include <cmath>

namespace vfep {

namespace {
constexpr double kTiny = 1e-12;

static bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}
} // namespace

double Suppression::apply(double dt,
                          Reactor& r,
                          double& inhib_kgm3,
                          double& inert_kgm3,
                          double& agent_mdot_kgps) {
    const double V = r.config().volume_m3;

    // Default outputs
    // - inhibitor_kgm3 is modeled as a scalar "loading" in the space (not part of reactor moles).
    // - inert_kgm3 is derived from reactor composition to remain consistent with ventilation removal.
    inhib_kgm3 = (V > 0.0) ? (st_.inhibitor_kg / V) : 0.0;
    inert_kgm3 = 0.0;
    agent_mdot_kgps = 0.0;

    // If volume is invalid, we cannot compute concentrations reliably.
    if (V <= 0.0) {
        return 0.0;
    }

    // Derive current inert concentration from reactor state (single source of truth).
    // This remains consistent even if ventilation later removes INERT moles.
    if (idx_.iINERT >= 0 && idx_.iINERT < static_cast<int>(r.moles().size())
        && idx_.iINERT < static_cast<int>(r.species().size())) {
        const auto& sp = r.species();
        const double M = sp[idx_.iINERT].molarMass_kg_per_mol;
        const double nInert = r.moles()[idx_.iINERT];
        if (M > kTiny && std::isfinite(nInert) && nInert > 0.0) {
            inert_kgm3 = (nInert * M) / V;
        } else {
            inert_kgm3 = 0.0;
        }
    }

    // Reject non-physical dt
    if (!isFinitePositive(dt)) return 0.0;

    // --------------------
    // Minimal VFEP actuator truth (deterministic)
    // --------------------
    // Ramp RPM toward target when enabled and tank has inventory; otherwise ramp down to zero.
    const double rpm_target = (cfg_.enabled && st_.tank_kg > kTiny) ? std::max(0.0, cfg_.rpm_target) : 0.0;
    const double rpm_rate   = std::max(0.0, cfg_.rpm_ramp_rate_rpmps);

    if (std::isfinite(rpm_target) && std::isfinite(rpm_rate)) {
        const double maxStep = rpm_rate * dt;
        const double err = rpm_target - st_.vfep_rpm;
        const double step = std::clamp(err, -maxStep, maxStep);
        st_.vfep_rpm = std::max(0.0, st_.vfep_rpm + step);
    } else {
        st_.vfep_rpm = 0.0;
    }


    // No agent available or disabled
    if (!cfg_.enabled || st_.tank_kg <= kTiny) return 0.0;

    // Limit mdot by remaining tank inventory
    const double mdot_cap = st_.tank_kg / dt; // kg/s to empty in exactly dt
    const double mdot = std::max(0.0, std::min(cfg_.mdot_total_kgps, mdot_cap));
    if (mdot <= kTiny) return 0.0;

    // Normalize fractions safely (so they sum to 1 when possible)
    double fracInhib = std::clamp(cfg_.frac_inhibitor, 0.0, 1.0);
    double fracInert = std::clamp(cfg_.frac_inert,     0.0, 1.0);

    const double fracSum = fracInhib + fracInert;
    if (fracSum > kTiny) {
        fracInhib /= fracSum;
        fracInert /= fracSum;
    } else {
        // If user set both to 0, default to all inert to avoid division-by-zero
        fracInhib = 0.0;
        fracInert = 1.0;
    }

    const double mdot_inhib = mdot * fracInhib;
    const double mdot_inert = mdot * fracInert;

    // Update tank and stored inhibitor mass (inhibitor modeled as space loading)
    const double dm = mdot * dt;
    st_.tank_kg = std::max(0.0, st_.tank_kg - dm);

    st_.inhibitor_kg += mdot_inhib * dt;

    // Note: st_.inert_kg is no longer used for inert_kgm3 reporting.
    // You may keep it as "delivered inert mass" for bookkeeping; it does not affect the reactor state.
    st_.inert_kg += mdot_inert * dt;

    inhib_kgm3 = st_.inhibitor_kg / V;
    agent_mdot_kgps = mdot;

    // Add inert gas to reactor mixture (O2 displacement)
    if (idx_.iINERT >= 0 && idx_.iINERT < static_cast<int>(r.moles().size())
        && idx_.iINERT < static_cast<int>(r.species().size())) {
        const auto& sp = r.species();
        const double M = sp[idx_.iINERT].molarMass_kg_per_mol;
        if (M > kTiny) {
            const double dn = (mdot_inert * dt) / M;
            r.addMoles(idx_.iINERT, dn);
        }
    }

    // Update inert_kgm3 again after injection to reflect current reactor state.
    if (idx_.iINERT >= 0 && idx_.iINERT < static_cast<int>(r.moles().size())
        && idx_.iINERT < static_cast<int>(r.species().size())) {
        const auto& sp = r.species();
        const double M = sp[idx_.iINERT].molarMass_kg_per_mol;
        const double nInert = r.moles()[idx_.iINERT];
        if (M > kTiny && std::isfinite(nInert) && nInert > 0.0) {
            inert_kgm3 = (nInert * M) / V;
        } else {
            inert_kgm3 = 0.0;
        }
    }

    // Cooling power returned as positive Watts for cooling
    return cfg_.cooling_W_per_kgps * mdot;
}

} // namespace vfep
