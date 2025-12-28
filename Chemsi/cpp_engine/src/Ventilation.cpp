// Ventilation.cpp (patched: conserve exchanged moles + energy-consistent temperature update)
// Key changes vs prior version:
// - Keeps your robust dnOut-based mole exchange (no long-run drift).
// - Replaces linear temperature interpolation with a Cp-weighted energy balance
//   consistent with Reactor's Cp model (Cp_mix [J/K] = Σ n_i * cp_i).
//
// Energy model used (consistent with constant molar cp per species, no reference enthalpy):
//   E := Σ (n_i * cp_i) * T  = Cp_mix * T
// Outflow removes species-specific Cp contributions at reactor temperature.
// Inflow adds Cp contributions at supply temperature.
// Then T_new = E_new / Cp_new.
//
// This eliminates "free" heating/cooling from ventilation and aligns with Reactor.cpp contracts.

#include "Ventilation.h"

#include <algorithm>
#include <cmath>

namespace vfep {

namespace {
constexpr double kTiny = 1e-15;

static bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}

// Clamp/renormalize a 3-component spec (O2, CO2, H2O) and return N2 remainder.
static void normalizeAmbientFractions(double& yO2, double& yCO2, double& yH2O, double& yN2) {
    yO2  = std::clamp(yO2,  0.0, 1.0);
    yCO2 = std::clamp(yCO2, 0.0, 1.0);
    yH2O = std::clamp(yH2O, 0.0, 1.0);

    const double sum = yO2 + yCO2 + yH2O;
    if (sum > 1.0) {
        const double inv = 1.0 / sum;
        yO2  *= inv;
        yCO2 *= inv;
        yH2O *= inv;
        yN2 = 0.0;
    } else {
        yN2 = 1.0 - sum;
    }
}

} // namespace

void Ventilation::apply(double dt, Reactor& r) {
    if (!isFinitePositive(dt)) return;

    // ACH [1/hr] -> lambda [1/s]
    const double lambda = std::max(0.0, cfg_.ACH) / 3600.0;
    if (lambda <= 0.0) return;

    // Exchange fraction over dt (clamped for stability; prefer substepping upstream for large dt)
    const double frac = std::clamp(lambda * dt, 0.0, 1.0);
    if (frac <= 0.0) return;

    // Supply composition (mole fractions)
    double yO2  = 0.2095;
    double yCO2 = 0.00042;
    double yH2O = 0.0100;
    double yN2  = 0.0;
    normalizeAmbientFractions(yO2, yCO2, yH2O, yN2);

    auto& n = r.moles();
    const auto& sp = r.species();

    const double T  = r.temperatureK();
    const double Ts = cfg_.T_supply_K;

    // If temperatures are invalid, we still do mole exchange but skip temperature update.
    const bool tempsOK = std::isfinite(T) && std::isfinite(Ts);

    // Energy proxy consistent with Reactor Cp model: E = (Σ n_i cp_i) * T
    double E_old = 0.0;
    if (tempsOK) {
        const double Cp_old = r.mixtureCp_J_per_K();
        if (std::isfinite(Cp_old) && Cp_old > kTiny) {
            E_old = Cp_old * T;
        } else {
            // If Cp is degenerate, skip energy update.
            // Mole exchange will still proceed (Reactor guards against Cp near zero).
            E_old = 0.0;
        }
    }

    // Remove a fraction of existing gas moles (outflow), and compute the actual removed moles.
    // Also compute Cp removed so we can remove energy consistently.
    double dnOut = 0.0;
    double Cp_out_J_per_K = 0.0; // Σ(dn_i * cp_i) for removed gas

    auto removeFrac = [&](int idx) {
        if (idx < 0 || idx >= static_cast<int>(n.size())) return;
        if (idx >= static_cast<int>(sp.size())) return;
        if (sp[idx].phase != Phase::Gas) return;

        const double ni = n[idx];
        if (!(ni > kTiny) || !std::isfinite(ni)) {
            // Sanitize non-finite / non-positive inventories
            n[idx] = 0.0;
            return;
        }

        const double dn = ni * frac;  // moles removed for this species
        dnOut += dn;

        if (tempsOK) {
            const double cpi = sp[idx].cp_J_per_molK;
            if (std::isfinite(cpi) && cpi > 0.0) {
                Cp_out_J_per_K += dn * cpi;
            }
        }

        n[idx] = std::max(0.0, ni - dn);
    };

    removeFrac(idx_.iN2);
    removeFrac(idx_.iO2);
    removeFrac(idx_.iCO2);
    removeFrac(idx_.iH2O);
    removeFrac(idx_.iFUEL);
    removeFrac(idx_.iINERT);

    if (!(dnOut > kTiny) || !std::isfinite(dnOut)) {
        return;
    }

    // Add back inflow moles at supply composition (exactly matching actual outflow)
    // Also compute Cp_in so we can add energy consistently.
    double Cp_in_J_per_K = 0.0; // Σ(dnIn_i * cp_i) for incoming gas

    auto addIn = [&](int idx, double y) {
        if (idx < 0 || idx >= static_cast<int>(n.size())) return;
        if (idx >= static_cast<int>(sp.size())) return;
        if (sp[idx].phase != Phase::Gas) return;

        const double yy = std::max(0.0, y);
        if (yy <= 0.0) return;

        const double dnIn = dnOut * yy;
        n[idx] += dnIn;

        if (tempsOK) {
            const double cpi = sp[idx].cp_J_per_molK;
            if (std::isfinite(cpi) && cpi > 0.0) {
                Cp_in_J_per_K += dnIn * cpi;
            }
        }
    };

    addIn(idx_.iN2,  yN2);
    addIn(idx_.iO2,  yO2);
    addIn(idx_.iCO2, yCO2);
    addIn(idx_.iH2O, yH2O);

    // In supply air, assume no fuel and no inert agent
    // (they are removed by outflow only)
    // idx_.iFUEL and idx_.iINERT intentionally not added.

    // Energy-consistent temperature update (Cp-weighted mixing)
    if (tempsOK) {
        const double Cp_new = r.mixtureCp_J_per_K();
        if (std::isfinite(Cp_new) && Cp_new > kTiny) {
            // Remove outflow energy at reactor temperature; add inflow energy at supply temperature.
            // E_old already represents Cp_old*T (if Cp_old was valid); if Cp_old was invalid, E_old=0
            // and we skip meaningful energy conservation—still avoids NaNs.
            const double E_out = Cp_out_J_per_K * T;
            const double E_in  = Cp_in_J_per_K  * Ts;

            double E_new = (E_old - E_out) + E_in;

            if (!std::isfinite(E_new)) {
                // If energy fails numerically, do not modify temperature.
                return;
            }

            const double T_new = E_new / Cp_new;
            if (std::isfinite(T_new)) {
                r.setTemperatureK(T_new);
            }
        }
    }
}

} // namespace vfep
