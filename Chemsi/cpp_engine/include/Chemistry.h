#pragma once

#include <vector>
#include "Species.h"

namespace vfep {

struct CombustionModel {
    // Fuel pseudo-formula CxHyOz
    double C = 1.0;
    double H = 2.0;
    double O = 0.0;

    // Kinetics (tunable)
    double A = 2.0e6;    // effective rate coefficient
    double Ea = 8.0e4;   // J/mol
    double orderFuel = 1.0;
    double orderO2   = 1.0;

    // Heat release (J/mol fuel reacted)
    double heatRelease_J_per_molFuel = 6.0e5;

    double nuO2()  const { return C + H / 4.0 - O / 2.0; }
    double nuCO2() const { return C; }
    double nuH2O() const { return H / 2.0; }
};

struct ChemistryIndex {
    int iN2    = -1;
    int iO2    = -1;
    int iCO2   = -1;
    int iH2O   = -1;
    int iFUEL  = -1;
    int iINERT = -1;
    int iINHIB = -1;
};

struct ReactionResult {
    double dMolFuel = 0.0;
    double dMolO2   = 0.0;
    double dMolCO2  = 0.0;
    double dMolH2O  = 0.0;
    double heat_W   = 0.0; // Positive = heat release rate
};

class Chemistry {
public:
    Chemistry(const std::vector<Species>& sp, ChemistryIndex idx, CombustionModel model);

    ReactionResult react(
        double dt,
        double T_K,
        double V_m3,
        std::vector<double>& n_mol,
        double inhibitor_kg_per_m3
    );

private:
    const std::vector<Species>& sp_;
    ChemistryIndex idx_;
    CombustionModel model_;
};

} // namespace vfep
