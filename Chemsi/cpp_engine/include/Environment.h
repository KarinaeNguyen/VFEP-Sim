#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <vector>

struct RoomConfig {
    // 1. Dimensions (The "Box")
    double length, width, height; // Meters

    // 2. Atmospheric Properties
    double gravity;       // -9.81 m/s^2
    double airDensity;    // kg/m^3 (Standard: 1.225)
    double temperature;   // Celsius (Global Ambient)
    double humidity;      // % (0.0 to 1.0)

    // 3. Airflow (AC Units)
    // A simple vector representing the constant wind from CRAC units
    double windVx, windVy, windVz;

    // Constructor with Default Data Center Values
    RoomConfig() 
        : length(10.0), width(8.0), height(4.0), // 10x8x4m Room
          gravity(-9.81),
          airDensity(1.225),
          temperature(24.0), // Standard Server Room Temp
          humidity(0.45),    // 45% Humidity
          windVx(0.0), windVy(0.0), windVz(0.0) {} // Still air initially
};

#endif