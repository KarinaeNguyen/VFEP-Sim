#ifndef PROJECTILE_H
#define PROJECTILE_H

#include "PurpleK.h" // Essential for the specific reaction

class Projectile {
public:
    PurpleK payload; // Storing the specific object, not the generic parent

    double x, y, z;
    double vx, vy, vz;
    double mass;

    // Constructor: Accepts PurpleK specifically
    Projectile(PurpleK agent, double startX, double startY, double startZ) 
        : payload(agent), x(startX), y(startY), z(startZ), vx(0), vy(0), vz(0) {
        
        this->mass = payload.molecularWeight;
    }

    void updatePhysics(double dt) {
        x += vx * dt;
        y += vy * dt;
        z += vz * dt;
    }
};

#endif