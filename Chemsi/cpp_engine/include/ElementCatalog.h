#ifndef ELEMENT_CATALOG_H
#define ELEMENT_CATALOG_H

#include "Atom.h"

// Factory class to spawn specific elements with correct physics data
class ElementCatalog {
public:
    // --- PURPLE K ELEMENTS ---
    
    // Potassium (K) - The Radical Scavenger
    static Atom CreatePotassium(double x, double y, double z) {
        // Mass: 39.098 amu, Radius: 275 pm, Charge: +1
        return Atom("K", 19, 39.098, 275.0, 1.0, 0.82, x, y, z);
    }

    // Hydrogen (H)
    static Atom CreateHydrogen(double x, double y, double z) {
        return Atom("H", 1, 1.008, 120.0, 1.0, 2.20, x, y, z);
    }

    // Carbon (C)
    static Atom CreateCarbon(double x, double y, double z) {
        return Atom("C", 6, 12.011, 170.0, 0.0, 2.55, x, y, z);
    }

    // Oxygen (O)
    static Atom CreateOxygen(double x, double y, double z) {
        return Atom("O", 8, 15.999, 152.0, -2.0, 3.44, x, y, z);
    }

    // --- RACK & INSULATION ELEMENTS (The Fuel/Obstacle) ---
    
    // Chlorine (Cl) - Found in PVC Insulation (C2H3Cl)
    static Atom CreateChlorine(double x, double y, double z) {
        return Atom("Cl", 17, 35.45, 175.0, -1.0, 3.16, x, y, z);
    }
    
    // Copper (Cu) - The Wiring
    static Atom CreateCopper(double x, double y, double z) {
        return Atom("Cu", 29, 63.546, 140.0, 0.0, 1.90, x, y, z);
    }
    
    // Iron (Fe) - The Rack Structure
    static Atom CreateIron(double x, double y, double z) {
        return Atom("Fe", 26, 55.845, 126.0, 0.0, 1.83, x, y, z);
    }
};

#endif