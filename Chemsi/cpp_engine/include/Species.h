#pragma once
#include <string>

namespace vfep {

enum class Phase { Gas, Aerosol, Solid };

struct Species {
    std::string name;
    double molarMass_kg_per_mol; // kg/mol
    double cp_J_per_molK;        // J/(mol*K) constant cp approximation
    Phase phase;
};

} // namespace vfep
