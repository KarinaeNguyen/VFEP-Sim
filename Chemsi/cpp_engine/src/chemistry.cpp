#include "Chemistry.h"
#include "Constants.h"

#include <algorithm>
#include <cmath>

namespace vfep {

namespace {
constexpr double kTiny = 1e-15;
constexpr double kMinTemp_K = 250.0;

// inhibitor effect: exp(-k_inhib * conc)
constexpr double kInhibCoeff = 5.0; // 1/(kg/m^3) tune later

static bool isFinitePositive(double x) {
    return std::isfinite(x) && x > 0.0;
}
} // namespace

Chemistry::Chemistry(const std::vector<Species>& sp, ChemistryIndex idx, CombustionModel model)
: sp_(sp), idx_(idx), model_(model) {}

ReactionResult Chemistry::react(
    double dt,
    double T_K,
    double V_m3,
    std::vector<double>& n_mol,
    double inhibitor_kg_per_m3
) {
    ReactionResult rr;

    // Basic guards
    if (!isFinitePositive(dt)) return rr;
    if (!isFinitePositive(V_m3)) return rr;
    if (!std::isfinite(T_K)) return rr;
    if (n_mol.empty()) return rr;

    const int iF   = idx_.iFUEL;
    const int iO2  = idx_.iO2;
    const int iCO2 = idx_.iCO2;
    const int iH2O = idx_.iH2O;

    if (iF < 0 || iO2 < 0 || iCO2 < 0 || iH2O < 0) return rr;
    if (iF >= static_cast<int>(n_mol.size())
        || iO2 >= static_cast<int>(n_mol.size())
        || iCO2 >= static_cast<int>(n_mol.size())
        || iH2O >= static_cast<int>(n_mol.size())) {
        return rr;
    }

    const double nFuel = std::max(0.0, n_mol[iF]);
    const double nO2   = std::max(0.0, n_mol[iO2]);
    if (nFuel <= kTiny || nO2 <= kTiny) return rr;

    // Concentrations (mol/m^3)
    const double cFuel = nFuel / V_m3;
    const double cO2   = nO2   / V_m3;

    // Arrhenius temperature clamp
    const double Tuse = std::max(kMinTemp_K, T_K);
    const double kT = model_.A * std::exp(-model_.Ea / (R_universal * Tuse));
    if (!std::isfinite(kT) || kT <= 0.0) return rr;

    const double inhib = std::max(0.0, inhibitor_kg_per_m3);
    const double inhibFactor = std::exp(-kInhibCoeff * inhib);

    // Rate of fuel consumption (mol/m^3/s)
    double rFuel = kT
        * std::pow(std::max(0.0, cFuel), model_.orderFuel)
        * std::pow(std::max(0.0, cO2),   model_.orderO2)
        * inhibFactor;

    if (!std::isfinite(rFuel) || rFuel <= 0.0) return rr;

    // Stoichiometry
    const double nuO2  = std::max(0.0, model_.nuO2());
    const double nuCO2 = std::max(0.0, model_.nuCO2());
    const double nuH2O = std::max(0.0, model_.nuH2O());

    const double maxFuelByFuel = nFuel;
    const double maxFuelByO2   = (nuO2 > kTiny) ? (nO2 / nuO2) : 0.0;

    const double fuelToConsume_kin = rFuel * V_m3 * dt; // mol
    const double fuelToConsume =
        std::max(0.0, std::min({fuelToConsume_kin, maxFuelByFuel, maxFuelByO2}));

    if (fuelToConsume <= kTiny) return rr;

    const double o2Consumed = nuO2  * fuelToConsume;
    const double co2Formed  = nuCO2 * fuelToConsume;
    const double h2oFormed  = nuH2O * fuelToConsume;

    // Update moles with non-negativity guards
    n_mol[iF]   = std::max(0.0, n_mol[iF]  - fuelToConsume);
    n_mol[iO2]  = std::max(0.0, n_mol[iO2] - o2Consumed);
    n_mol[iCO2] = std::max(0.0, n_mol[iCO2] + co2Formed);
    n_mol[iH2O] = std::max(0.0, n_mol[iH2O] + h2oFormed);

    rr.dMolFuel = -fuelToConsume;
    rr.dMolO2   = -o2Consumed;
    rr.dMolCO2  = +co2Formed;
    rr.dMolH2O  = +h2oFormed;

    // Heat release rate (W). Positive means heat generation.
    const double Q_J = model_.heatRelease_J_per_molFuel * fuelToConsume;
    rr.heat_W = std::isfinite(Q_J) ? (Q_J / dt) : 0.0;

    return rr;
}

} // namespace vfep
