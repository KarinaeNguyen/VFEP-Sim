#ifndef RACK_H
#define RACK_H

#include <string>
#include "ElementCatalog.h" // To define materials like Iron/Copper

class ServerRack {
public:
    int id;
    std::string label; // e.g., "Rack-A01"
    
    // Position (Bottom-Left-Back Corner)
    double x, y, z;
    
    // Dimensions
    double width, depth, height;

    // State
    bool isOnFire;
    double surfaceTemp;

    // Constructor
    ServerRack(int id, std::string label, double x, double y, double z) 
        : id(id), label(label), x(x), y(y), z(z), 
          width(0.6), depth(1.0), height(2.2), // Standard 42U Rack Size (m)
          isOnFire(false), surfaceTemp(30.0) {}

    // 📐 COLLISION DETECTION (The "Hit Box")
    // Checks if a point (px, py, pz) is inside this rack
    bool checkCollision(double px, double py, double pz) {
        return (px >= x && px <= x + width) &&
               (py >= y && py <= y + depth) &&
               (pz >= z && pz <= z + height);
    }

    // 🔥 IGNITION METHOD
    // Turns this rack into a heat source
    void ignite() {
        isOnFire = true;
        surfaceTemp = 800.0; // Cable fire temp
    }
};

#endif