#pragma once

#include <vector>

#include "Species.h"
#include "Chemistry.h"

namespace vfep {

struct ReactorConfig {
    double volume_m3   = 120.0;
    double area_m2     = 180.0;
    double T_amb_K     = 295.15;
    double h_W_m2K     = 10.0;
    double emissivity  = 0.85;
};

class Reactor {
public:
    Reactor(std::vector<Species> sp, ChemistryIndex idx, CombustionModel model);

    // IMPORTANT:
    // Chemistry stores references to the species vector; default moves can invalidate them.
    // For safety and determinism (especially in UI/visualizer contexts), Reactor is non-copyable and non-movable.
    Reactor(const Reactor&) = delete;
    Reactor& operator=(const Reactor&) = delete;
    Reactor(Reactor&&) = delete;
    Reactor& operator=(Reactor&&) = delete;

    void setConfig(const ReactorConfig& cfg) noexcept { cfg_ = cfg; }
    const ReactorConfig& config() const noexcept { return cfg_; }

    // Species list is immutable after construction to preserve indexing and cached gasIdx_.
    const std::vector<Species>& species() const noexcept { return sp_; }

    // Moles vector: must remain the same size as species() for the lifetime of Reactor.
    // Invariant expected by all subsystems:
    //   n_mol_.size() == sp_.size()
    std::vector<double>& moles() noexcept { return n_mol_; }
    const std::vector<double>& moles() const noexcept { return n_mol_; }

    double temperatureK() const noexcept { return T_K_; }
    void setTemperatureK(double T) noexcept { T_K_ = T; }

    // Totals/mixture properties operate over gas species only.
    [[nodiscard]] double totalGasMoles() const noexcept;

    // Mixture heat capacity [J/K] computed from molar cp [J/mol-K] and moles [mol].
    // Expected to be finite and >= 0.0. Reactor::step() must guard against Cp near zero.
    [[nodiscard]] double mixtureCp_J_per_K() const noexcept;

    // Returns gas-phase mole fraction in [0,1]. NaN-safe:
    // if totals are zero/non-finite or i is invalid/non-gas, returns 0.0.
    [[nodiscard]] double gasMoleFraction(int i) const noexcept;

    // Adds dn moles to species i. dn may be negative; inventory is clamped at >= 0.
    void addMoles(int i, double dn) noexcept;

    // externalCooling_W convention:
    //   > 0 removes heat (cooling)
    //   < 0 adds heat (heating)
    //
    // outCombustionHRR_W: Positive = heat release rate from combustion only (W).
    void step(double dt,
              double inhibitor_kg_per_m3,
              double externalCooling_W,
              double& outCombustionHRR_W) noexcept;

private:
    // Positive = heat leaving reactor to ambient (W)
    [[nodiscard]] double heatLoss_W() const noexcept;

    ReactorConfig cfg_;
    std::vector<Species> sp_;
    ChemistryIndex idx_;
    CombustionModel model_;
    Chemistry chemistry_;

    std::vector<double> n_mol_;
    std::vector<int> gasIdx_; // cached indices of gas species
    double T_K_ = 295.15;
};

} // namespace vfep
