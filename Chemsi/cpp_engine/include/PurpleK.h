#ifndef PURPLEK_H
#define PURPLEK_H

#include "ChemicalAgent.h"
#include "ElementCatalog.h"

class PurpleK : public ChemicalAgent {
public:
    PurpleK() : ChemicalAgent("Purple K (KHCO3)") {
        // Build the KHCO3 Molecule using our Catalog
        addAtom(ElementCatalog::CreatePotassium(0,0,0));
        addAtom(ElementCatalog::CreateHydrogen(50,0,0));
        addAtom(ElementCatalog::CreateCarbon(100,0,0));
        addAtom(ElementCatalog::CreateOxygen(150,50,0));
        addAtom(ElementCatalog::CreateOxygen(150,-50,0));
        addAtom(ElementCatalog::CreateOxygen(200,0,0));
    }

    // --- THE REACTION LOGIC ---
    // Returns a list of NEW chemical agents (The Products)
    // 2 KHCO3 -> K2CO3 + CO2 + H2O
    std::vector<ChemicalAgent> decompose(double temperature) {
        std::vector<ChemicalAgent> products;

        if (temperature > 190.0) { // Celsius
            // 1. Create Water Mist (H2O)
            ChemicalAgent water("Water Mist (H2O)");
            water.addAtom(ElementCatalog::CreateHydrogen(0,0,0));
            water.addAtom(ElementCatalog::CreateHydrogen(20,0,0));
            water.addAtom(ElementCatalog::CreateOxygen(10,10,0));
            products.push_back(water);

            // 2. Create CO2 Gas
            ChemicalAgent co2("Carbon Dioxide (CO2)");
            co2.addAtom(ElementCatalog::CreateCarbon(0,0,0));
            co2.addAtom(ElementCatalog::CreateOxygen(30,0,0));
            co2.addAtom(ElementCatalog::CreateOxygen(-30,0,0));
            products.push_back(co2);

            // 3. Create Potassium Carbonate (K2CO3) Residue
            ChemicalAgent residue("Potassium Carbonate (K2CO3)");
            residue.addAtom(ElementCatalog::CreatePotassium(0,0,0));
            residue.addAtom(ElementCatalog::CreatePotassium(50,0,0));
            residue.addAtom(ElementCatalog::CreateCarbon(25,25,0));
            residue.addAtom(ElementCatalog::CreateOxygen(25,50,0));
            residue.addAtom(ElementCatalog::CreateOxygen(25,0,0));
            residue.addAtom(ElementCatalog::CreateOxygen(50,25,0));
            products.push_back(residue);
        }
        
        return products;
    }
};

#endif