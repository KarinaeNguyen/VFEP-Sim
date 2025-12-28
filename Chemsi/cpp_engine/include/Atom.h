#ifndef ATOM_H
#define ATOM_H

#include <string>
#include <cmath>
#include <utility>

// Safer scoped enum
enum class AtomState { SOLID, LIQUID, GAS, PLASMA };

class Atom {
public:
    // ------------------------------------------------------------------
    // Backward-compatible public fields (used by tests)
    // ------------------------------------------------------------------
    int id = 0;
    std::string type;

    // ------------------------------------------------------------------
    // Identity
    // ------------------------------------------------------------------
    std::string symbol;
    int atomicNumber = 0;

    // ------------------------------------------------------------------
    // Physical / chemical properties
    // ------------------------------------------------------------------
    double mass = 0.0;
    double radius = 0.0;
    AtomState state = AtomState::SOLID;
    double charge = 0.0;
    double electronegativity = 0.0;

    // ------------------------------------------------------------------
    // Position / velocity
    // ------------------------------------------------------------------
    double x = 0.0, y = 0.0, z = 0.0;
    double vx = 0.0, vy = 0.0, vz = 0.0;

    // ------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------

    // (id, type, mass, x, y, z)
    Atom(int id_,
         const char* type_,
         double mass_,
         double px,
         double py,
         double pz)
        : id(id_),
          type(type_ ? type_ : ""),
          symbol(type_ ? type_ : ""),
          atomicNumber(0),
          mass(mass_),
          radius(0.0),
          state(AtomState::SOLID),
          charge(0.0),
          electronegativity(0.0),
          x(px), y(py), z(pz),
          vx(0.0), vy(0.0), vz(0.0)
    {}

    Atom(std::string s,
         int z_num,
         double m,
         double r,
         double c,
         double en,
         double px,
         double py,
         double pz)
        : id(0),
          type(s),
          symbol(std::move(s)),
          atomicNumber(z_num),
          mass(m),
          radius(r),
          state(AtomState::SOLID),
          charge(c),
          electronegativity(en),
          x(px), y(py), z(pz),
          vx(0.0), vy(0.0), vz(0.0)
    {}

    // ------------------------------------------------------------------
    // Helpers
    // ------------------------------------------------------------------
    double getVolume() const {
        constexpr double pi = 3.14159265358979323846;
        return (4.0 / 3.0) * pi * std::pow(radius, 3);
    }
};

#endif // ATOM_H
