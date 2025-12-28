#ifndef CHEMICAL_AGENT_H
#define CHEMICAL_AGENT_H

#include <vector>
#include <string>
#include "Atom.h"

// The Middle Layer: Represents one unit of the Chemical Compound
class ChemicalAgent {
public:
    std::string name;
    std::vector<Atom> molecularStructure;
    double molecularWeight;

    ChemicalAgent(std::string n) : name(n), molecularWeight(0) {}

    void addAtom(Atom a) {
        molecularStructure.push_back(a);
        molecularWeight += a.mass;
    }

    // --- NEW: VIRTUAL DECOMPOSE FUNCTION ---
    // The "virtual" keyword tells C++: "Check if the child class (PurpleK) has a better version of this function."
    virtual std::vector<ChemicalAgent> decompose(double temperature) {
        // Default behavior: Do nothing (return empty list)
        return {}; 
    }
};

#endif